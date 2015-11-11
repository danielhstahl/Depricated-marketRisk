#define _USE_MATH_DEFINES
#include <iostream>
//#include "Swap.h"
#include <ctime>
#include <vector>
#include <chrono> //for accurate multithreading time using std::chrono
//#include "EulerSimulation.h"
//#include "Matrix.h"
//#include "Double.h"
#include <map>
#include "SimulNorm.h"
#include <unordered_map>
#include "Date.h"
#include "BlackScholes.h"
#include "Vasicek.h"
#include "MarketData.h"
#include "Newton.h"
#include "MC.h"
#include "YieldSpline.h"

int main(){
  MC mc(1000);//100000 run initially
  /*create portfolio...note that this is a "fake" portfolio for demonstration*/
  int sizeOfPortfolio=3000;//3000 distinct assets
  Date currDate;
  std::vector<AssetFeatures> portfolio;
  for(int i=0; i<sizeOfPortfolio/4;i++){
    //Date dt=currDate+(i+1)*.005;
    portfolio.push_back(AssetFeatures(currDate+(i+1)*.005, "bond"));
  }
  for(int i=0; i<sizeOfPortfolio/4;i++){
    portfolio.push_back(AssetFeatures(currDate+(i+1)*.005, .03, .25, "caplet"));
  }
  for(int i=0; i<sizeOfPortfolio/4;i++){
    portfolio.push_back(AssetFeatures(currDate+(i+1)*.005, .03, .25, currDate+((i+1)*.001+2), "swaption"));
    //Strike, Tenor, Tenor, SwapMaturity, OptionMaturity
  }
  for(int i=0; i<sizeOfPortfolio/4;i++){
    portfolio.push_back(AssetFeatures(currDate+(i+1)*.005, .99, currDate+((i+1)*.001+.25), "call"));
  }
  //std::cout<<"test1: "<<portfolio[0].Maturity<<std::endl;
  //std::cout<<"test2: "<<portfolio[300].Maturity<<std::endl;
  currDate.setScale("day");
  Date simulateToDate=currDate+10;//ten day VaR
  double future=simulateToDate-currDate; //convert date to double
  std::map<double, double> holdPossibleDates; //using "map" to ensure sorted keys
  double diff=0;
  for(int i=0; i<sizeOfPortfolio;i++){
    if(portfolio[i].Maturity-simulateToDate<=0){//if any maturities are less than the simulated date
      diff=portfolio[i].Maturity-currDate;
      holdPossibleDates[diff]=diff;
    }
  }
  diff=simulateToDate-currDate;
  holdPossibleDates[diff]=diff;
 /*end portfolio*/

  double a=.4; //speed
  double b=.05; //long run average
  double sig=.01; //volatility of the process
  double t=1;
  double r0=.03;
  double delta=.25;

  /*Newton Test */
  /*This demonstrates the workflow for estimating and using a one dimensional Vasicek model for VaR purposes in a fixed income portfolio.  The risk neutral parameter is assumed to be "mu".  Hence the "speed" and "volatility" are constant between the risk neutral and real world Vasicek processes. */

  //First, need to estimate volatility (by Girsonav's theorem, this doesn't change under measure change for a diffusion)
  int m=10;//number of volatility dates...calibrate to volatility smile

  //Create "fake" volatility surface.  Note that this is created from the assumption that the actual generating process is indeed a Vasicek process.  In a real world implementation, this data would be given to us by the market and we wouldn't "create" it
  std::vector<ForwardValue> volatilitySurface; //careful...currently this isnt the "Black" volatility surface using the "Black" model.  This is the volatility surface of implied volatility on bond options.  This can be converted to caplets easily; but it may be better to have a function that switches between Black volatility and Bond volatility
  Date volDate;
  volDate.setScale("year");
  double volDateD=0;
  for(int i=0;i<m;i++){
    volDateD=(i+1)*.25;
    volatilitySurface.push_back(ForwardValue(volDate+volDateD, volDate+(volDateD+delta), Vasicek_Volatility(a, sig,volDateD+delta, volDateD)));
  }

  //create "fake" yield curve.  Note that this is created from the assumption that the yield curve is generated by bonds that are "derivatives" of a Vasicek process.
  std::vector<SpotValue> yield;
  Date dt;
  dt.setScale("year");
  double delt=0;
  for(int i=0; i<20; i++){
    delt=(i+1)*.25;

    yield.push_back(SpotValue(dt+delt, -log(Vasicek_Price(r0, a, b, sig, delt))/delt));
  }
  //HandleYield<NelsonSiegel> hy(yield); //uses nelsonsiegel method to find yield
  //Now we actually fit the data.  Note that the constructor automatically fits "Speed" and "Volatility" to the volatility surface using newton's method.  The constructor also estimates the "Theta" that fits the yield curve.
  YieldSpline yld(yield);
  Vasicek<YieldSpline> vs(yld, volatilitySurface, r0); //this construtor currently prints the estimate of "theta".  In this example, these estimates are constant (since theta(t)=mu for all t)
  /*for(int i=0; i<20; i++){
    delt=(i+1)*.25;

    std::cout<<dt+delt<<": ";

    std::cout<<-log(Vasicek_Price(r0, a, b, sig, delt))/delt<<" ";
    std::cout<<yld.Yield(delt)/delt<<" ";
    std::cout<<yld.Forward(delt)<<std::endl;
  }
*/
  /*some tests to ensure correct computations */
  std::cout<<vs.Bond_Price(t)<<std::endl; //This should be the same as the below
  std::cout<<Vasicek_Price(r0, a, b, sig, t)<<std::endl;  //This should be the same as above.
  std::cout<<Vasicek_Caplet(r0, a, b, sig, .04, delta, t)<<std::endl; //option on simple interest rate
  std::cout<<vs.Caplet(r0, .04, .04, delta, t)<<std::endl;

  std::cout<<vs.Caplet(r0, 0, .04, delta, t)<<std::endl;
  std::cout<<vs.Bond_Call(.02, 0, .92, 3, 2)<<std::endl;
  std::cout<<vs.Bond_Put(.01, 0, .99, 3, 2)<<std::endl;
  std::cout<<vs.Swap_Rate(.25, 4.5)<<std::endl;
  std::cout<<"Swaption: "<<vs.Swaption(r0, 0, .04, .25, 4.5, .5)<<std::endl;//Rate, FutureTime, Strike, Tenor, SwapMaturity, OptionMaturity
  std::vector<double> cashFl;
  for(int i=0; i<15; i++){
    cashFl.push_back(.5+.25*(i+1));
  }
  std::cout<<"Swaption Analytic: "<<vs.Swaption(r0, 0, .04, cashFl, .5)<<std::endl;
  std::cout<<"Future: "<<vs.EuroDollarFuture(r0, .04, .25, 1)<<std::endl;
  std::cout<<"Forward: "<<vs.Forward(r0, .04, .25, 1)<<std::endl;
  std::cout<<"Coupon Bond Option: "<<vs.Bond_Put(.03, .04, 1, .01, cashFl, .5)<<std::endl;
  std::vector<double> cashFlowBond;
  for(int i=0; i<15; i++){
    cashFlowBond.push_back(.25*(i+1));
  }
  std::cout<<"Coupon Bond: "<<vs.Bond_Price(.05, .5, .01, cashFl)<<std::endl;

  /*end tests */
  SimulNorm simul;
  vs.setFutureTimes(holdPossibleDates);
  /*test current portfolio */

  double portVal=0;
  for(int i=0; i<sizeOfPortfolio; ++i){
    double diffDate=portfolio[i].Maturity-currDate;
    double r=r0;
    double t=0;
    if(portfolio[i].type=="bond"){ //concerned here...I want this to be more automated (eg, simply call bond price automatically, but how to deal with parameters?)

      portVal+=vs.Bond_Price(r, t, portfolio[i].Maturity-currDate);//zero coupons remember...but are we assuming reinvestment??  Good grief that could get complicated...
    }
    else if(portfolio[i].type=="caplet"){
      portVal+=vs.Caplet(r, t, portfolio[i].Strike, portfolio[i].Tenor, portfolio[i].Maturity-currDate);
    }
    else if(portfolio[i].type=="swaption"){
      portVal+=vs.Swaption(r, t, portfolio[i].Strike, portfolio[i].Tenor, portfolio[i].UnderlyingMaturity-currDate, portfolio[i].Maturity-currDate);
    }
    else if(portfolio[i].type=="call"){
      portVal+=vs.Bond_Call(r, t, portfolio[i].Strike, portfolio[i].UnderlyingMaturity-currDate, portfolio[i].Maturity-currDate);
    }
  }
  std::cout<<"Current value of portfolio: "<<portVal<<std::endl;
  auto start = std::chrono::system_clock::now();
  mc.simulateDistribution([&](){
    std::unordered_map<double, double> ratePath=vs.simulate(simul); //simulates path
    double valueOfPortfolio=0;
    //bond, caplet, call
    double r=0;
    double t=0;
    for(int i=0; i<sizeOfPortfolio; ++i){
      double diffDate=portfolio[i].Maturity-currDate;
      r=0;
      t=0;
      if(ratePath.find(diffDate)!=ratePath.end()){ //if maturity is less than simulation
        r=ratePath.find(diffDate)->second;
        t=diffDate;
      }
      else{
        r=ratePath.find(future)->second;
        t=future;
      }
    //  std::cout<<r<<std::endl;
      if(portfolio[i].type=="bond"){ //concerned here...I want this to be more automated (eg, simply call bond price automatically, but how to deal with parameters?)
        valueOfPortfolio+=vs.Bond_Price(r, t, portfolio[i].Maturity-currDate);//zero coupons remember...but are we assuming reinvestment??  Good grief that could get complicated...
      }
      else if(portfolio[i].type=="caplet"){
        valueOfPortfolio+=vs.Caplet(r, t, portfolio[i].Strike, portfolio[i].Tenor, portfolio[i].Maturity-currDate);
      }
      else if(portfolio[i].type=="swaption"){
        valueOfPortfolio+=vs.Swaption(r, t, portfolio[i].Strike, portfolio[i].Tenor, portfolio[i].UnderlyingMaturity-currDate, portfolio[i].Maturity-currDate);
      }
      else if(portfolio[i].type=="call"){
        valueOfPortfolio+=vs.Bond_Call(r, t, portfolio[i].Strike, portfolio[i].UnderlyingMaturity-currDate, portfolio[i].Maturity-currDate);
      }
    }
    return valueOfPortfolio;
  });
  auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
  std::cout<<mc.getVaR(.99)<<std::endl;
  std::cout<<mc.getEstimate()<<std::endl;
  std::cout<<mc.getError()<<std::endl;

}
