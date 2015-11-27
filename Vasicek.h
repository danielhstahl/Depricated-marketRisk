#ifndef __VASICEK_H_INCLUDED__
#define __VASICEK_H_INCLUDED__
#include <vector>
#include <cmath>
#include "BlackScholes.h"
#include "MarketData.h"
#include "Newton.h"
#include "SimulNorm.h"
#include "Swap.h"
#include "BondUtilities.h"
#include <map>
#include <unordered_map>
#include <algorithm>

typedef double Tenor;
typedef double Price;
typedef double Mu;
typedef double Speed;//mean reversion
typedef double ShortRateSigma;
typedef double BondMaturity;
typedef double SwapMaturity;
typedef double Rate;
typedef double SwapRate;
typedef double FutureTime;
typedef double Coupon;
typedef double Maturity;
typedef std::vector<double> Times;

typedef std::vector<SpotValue> YieldCurve; //SpotValue is defined in "MarketData" as {Date, double}
typedef std::vector<ForwardValue> VolatilityCurve; //ForwardValue is defined in "MarketData" as {Date, Date, double}

//vanilla implementations
Discount Vasicek_Price(Rate, Speed, Mu, ShortRateSigma, BondMaturity);
//Discount Vasicek_Price(Rate, Speed, YieldCurve&, ShortRateSigma, BondMaturity);

double Vasicek_Volatility(Speed a, Mu sigma, BondMaturity tM, OptionMaturity t);
double Vasicek_Put(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Put(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
//double Vasicek_Swaption(Rate, Speed, Mu, ShortRateSigma, Strike, std::vector<double>&);
double Vasicek_Caplet(Rate, Speed, Mu, ShortRateSigma, Strike, Tenor, OptionMaturity);
double Vasicek_Caplet(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, Tenor, OptionMaturity); //underlying here is B(t, T+\delta), discount is B(t, T)



/*BSSigma Vasicek_Volatility (Speed a, ShortRateSigma sigma, BondMaturity tM, OptionMaturity t){ //can plug this into Black Scholes for option on bond
  return (1-exp(-a*(tM-t)))*(sigma/a)*sqrt((1-exp(-2*a*t))/(2*a*t));
}*/

template <typename T>
class Vasicek{
  private:
    YieldCurve yield; //requires continuously compounded yield
    VolatilityCurve vCurve;
    Rate r0;
    Speed a;
    Mu b;
    ShortRateSigma sigma;
    T* yieldClass;
    Times times;
    std::unordered_map<int, std::vector<double> > storeParameters;
    double minDiffT;
    int n;

    Date initial_t;
    double maxMaturity=0; //this is the "largest" maturity available

    void estimateSpeedVolatility(){
      estimateSpeedVolatility(.3, .1);
    }

    void estimateSpeedVolatility(double guessMu, double guessSigma){
      Newton nt;
      std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
      int n=vCurve.size();
      std::vector<double> dataToMinimizeOver(n); //the volatility from above
      std::vector<std::vector<double> > additionalParameters(n, std::vector<double>(2)); //a and sigma
      for(int i=0; i<n; i++){ //populate our functions to minimize over.
        vCurve[i].beginDate.setScale("year");
        vCurve[i].endDate.setScale("year");
        dataToMinimizeOver[i]=vCurve[i].value;
        Date currDate;
        additionalParameters[i][0]=vCurve[i].endDate-currDate;
        additionalParameters[i][1]=vCurve[i].beginDate-currDate;
        meanSquare.push_back([](std::vector<double> &guess, std::vector<double> &additionalParameters){
          return Vasicek_Volatility(guess[0], guess[1], additionalParameters[0], additionalParameters[1]);
        });
      }
      std::vector<double> guess(2);
      guess[0]=guessMu; //these seem decent guesses
      guess[1]=guessSigma;//these seem decent guesses
      nt.optimize(meanSquare, dataToMinimizeOver, guess, additionalParameters);
      a=guess[0];
      sigma=guess[1];
    }

    //double CtT(double, double);
  //std::unordered_map<int, std::unordered_map<int, double> > ctTMap;
    //std::unordered_map<int, double> diffusionMap;
  public:
    Vasicek(Speed a_, ShortRateSigma sigma_, Rate r0_){ //note...this is used for getting the "futures adjustment" section.
      a=a_;
      sigma=sigma_;
      r0=r0_;
    }
    Vasicek(Speed a_, ShortRateSigma sigma_){ //note...this is used for getting the "futures adjustment" section.
      a=a_;
      sigma=sigma_;
    }
    Vasicek(VolatilityCurve &vCurve_){ //note...this is used for getting the "futures adjustment" section.
      vCurve=vCurve_;
      estimateSpeedVolatility();
    }
    Vasicek(VolatilityCurve &vCurve_, Rate r0_){ //note...this is used for getting the "futures adjustment" section.
      vCurve=vCurve_;
      r0=r0_;
      estimateSpeedVolatility();
    }
    Vasicek(T &yield_, VolatilityCurve &vCurve_, Rate r0_){
      vCurve=vCurve_;
      r0=r0_;
      estimateSpeedVolatility();
      yieldClass=yield_;

    }

    Vasicek(T *yield_, Speed a_, ShortRateSigma sigma_, Rate r0_){
      a=a_;
      //yield=yield_;
      sigma=sigma_;
      r0=r0_;
      yieldClass=yield_;
      //createContinuousYield();
      //estimateTheta();
    }
    void setYield(T *yield_){
      yieldClass=yield_;
    }
    void setR(double r0_){
      r0=r0_;
    }
    Discount Bond_Price(Rate r, FutureTime t1, BondMaturity t2) {
      if(t2<=t1){
        return 1;
      }
      double ctT=(1-exp(-a*(t2-t1)))/a;
      Discount atT=-ctT*r;
      double yieldFunc1=yieldClass->Yield(t1);

      double yieldFunc2=yieldClass->Yield(t2);
      ctT=ctT*yieldClass->Forward(t1)-(sigma*sigma/(4*a))*ctT*ctT*(1-exp(-2*a*t1));
      //std::cout<<"yield1: "<<yieldFunc1;
      //std::cout<<" yield2: "<<yieldFunc2;
      //std::cout<<" forward: "<<yieldClass->Forward(t1)<<std::endl;
      return exp(atT+ctT-yieldFunc2+yieldFunc1); //see hull and white's seminal paper
    }
    Discount Bond_Price(BondMaturity t2){
      //return exp(-get_yield(t2));
      return exp(-yieldClass.Yield(t2));
    }
    Price Bond_Put(Rate r, FutureTime t0, Strike k, BondMaturity tM, OptionMaturity t){
      return Vasicek_Put(Bond_Price(r, t0, tM), Bond_Price(r, t0, t), r, a, sigma, k, tM-t0, t-t0);
    }
    Price Bond_Put(Strike k , BondMaturity tM, OptionMaturity t){
      return Vasicek_Put(Bond_Price(tM), Bond_Price(t), r0, a, sigma, k, tM, t);
    }
    Price Bond_Call(Rate r, FutureTime t0, Strike k, BondMaturity tM, OptionMaturity t){
      return Vasicek_Call(Bond_Price(r, t0, tM), Bond_Price(r, t0, t), r, a, sigma, k, tM-t0, t-t0);
    }
    Price Bond_Call(Strike k , BondMaturity tM, OptionMaturity t){
      return Vasicek_Call(Bond_Price(tM), Bond_Price(t), r0, a, sigma, k, tM, t);
    }
    Price Caplet(Rate r, FutureTime t0, Strike k , Tenor delta, OptionMaturity t){
      //std::cout<<"max t: "<<delta+t<<std::endl;
      return Vasicek_Caplet(Bond_Price(r, t0, t+delta), Bond_Price(r, t0, t), r, a, sigma, k, delta, t-t0);
    }

    Price Caplet(Strike k , Tenor delta, OptionMaturity t){
      return Vasicek_Caplet(Bond_Price(t+delta), Bond_Price(t), r0, a, sigma, k, delta, t);
    }
    Price Forward(Rate r, FutureTime t0, Tenor delta, Maturity t){
      return (Bond_Price(r, t0, t)/Bond_Price(r, t0, t+delta)-1)/delta;
    }
    Price Bond_Price(Rate r, FutureTime t0, Coupon cp, std::vector<double>& cashFlows){
      int n=cashFlows.size();
      double retVal=0;
      int l=0;
      while(cashFlows[l]<t0){
        l++;
      }
      for(int i=l; i<(n-1); i++){
        retVal+=cp*Bond_Price(r, t0, cashFlows[i]);
      }
      retVal+=(1+cp)*Bond_Price(r, t0, cashFlows[n-1]);
      return retVal;
    }
    void findHistoricalB(std::vector<SpotValue>& historicalRates){//assumes that real world rate process follows Vasicek
      int n=historicalRates.size();
      double dt=0;
      b=0;//global b
      double shortRateTime=7.0/360.0;
      historicalRates[0].value=convertLiborToContinuous(historicalRates[0].value, shortRateTime);
      for(int i=0; i<(n-1); ++i){
        historicalRates[i+1].date.setScale("year");
        dt=historicalRates[i+1].date-historicalRates[i].date;
        historicalRates[i+1].value=convertLiborToContinuous(historicalRates[i+1].value, shortRateTime); //convert to continuous time
        b+=(historicalRates[i+1].value-exp(-dt)*historicalRates[i].value)/(1.0-exp(-dt));
      }
      b=b/(n-1);
      std::cout<<"This is b: "<<b<<std::endl;
    }
    Price Bond_Call(Rate r, FutureTime t0, Strike k, Coupon cp, std::vector<double>& cashFlowTimes, OptionMaturity optMat){
      double xStrike=0;
      Newton nt;
      //std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
      int n=cashFlowTimes.size();
      double guess=r;
      nt.zeros([&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlowTimes[i]);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlowTimes[n-1]);
        return retVal-k;
      },[&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlowTimes[i])*((exp(-a*(cashFlowTimes[i]-optMat))-1)/a);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlowTimes[n-1])*((exp(-a*(cashFlowTimes[n-1]-optMat))-1)/a);
        return retVal;
      },guess);
      double retVal=0;
      for(int i=0; i<(n-1); i++){
        retVal+=cp*Bond_Call(r, t0, Bond_Price(guess, optMat, cashFlowTimes[i]), cashFlowTimes[i], optMat);
        //std::cout<<"retVal: "<<retVal<<" i: "<<i<<std::endl;
      }
      retVal+=(1+cp)*Bond_Call(r, t0, Bond_Price(guess, optMat, cashFlowTimes[n-1]), cashFlowTimes[n-1], optMat);
      return retVal;
    }
    Price Bond_Put(Rate r, FutureTime t0, Strike k, Coupon cp, std::vector<double>& cashFlowTimes, OptionMaturity optMat){
      double xStrike=0;
      Newton nt;
      //std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
      int n=cashFlowTimes.size();
      double guess=r;
      nt.zeros([&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlowTimes[i]);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlowTimes[n-1]);
        return retVal-k;
      },[&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlowTimes[i])*((exp(-a*(cashFlowTimes[i]-optMat))-1)/a);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlowTimes[n-1])*((exp(-a*(cashFlowTimes[n-1]-optMat))-1)/a);
        return retVal;
      },guess);
      double retVal=0;
      for(int i=0; i<(n-1); i++){
        retVal+=cp*Bond_Put(r, t0, Bond_Price(guess, optMat, cashFlowTimes[i]), cashFlowTimes[i], optMat);
        //std::cout<<"retVal: "<<retVal<<" i: "<<i<<std::endl;
      }
      retVal+=(1+cp)*Bond_Put(r, t0, Bond_Price(guess, optMat, cashFlowTimes[n-1]), cashFlowTimes[n-1], optMat);
      return retVal;
    }
    Price EuroDollarFuture(Rate r, FutureTime t0, Tenor delta, Maturity t){
      //double dT=1-exp(-a*t);
      //double dt=1-exp(-a*t0);
      double dtT=exp(-a*(t-t0));

      double ATTDelta=(1-exp(-a*delta))/a;
      double expT=exp(-a*t0);
      double forwardT=yieldClass.Forward(t);
      double mutT=r*dtT+forwardT-dtT*yieldClass.Forward(t0)+(sigma*sigma*.5/(a*a))*(1-dtT+exp(-2*a*t)-exp(-a*(t+t0)));
      double sigT=(sigma*sigma*.5/a)*(1-exp(-a*2*(t-t0)));
      double yieldFunc1=yieldClass.Yield(t);
      double yieldFunc2=yieldClass.Yield(t+delta);
      return (exp(yieldFunc2-yieldFunc1+ATTDelta*mutT+.5*ATTDelta*ATTDelta*sigT-ATTDelta*forwardT+(sigma*sigma*.25/a)*ATTDelta*ATTDelta*(1-exp(-2*t)))-1)/delta;

      //ctT=yieldFunc1-yieldFunc2+atT*yieldClass.Forward(t)-(sigma*sigma*.25/(a*a*a))*ctT*ctT*(exp(2*a*t)-1);
      //return (exp(atT*muT+.5*sigT*atT*atT-ctT)-1)/delta;
    }
    /*Price EuroDollarFutureAdjustment(Rate r, FutureTime t0, Tenor delta, Maturity t){
      double delT=t-t0;
      double dt=1-exp(-a*delT);
      double dtT=exp(-a*delta);
      double valToReturn=(1-dtT)/(a);
      valToReturn=(valToReturn/delta)*(valToReturn*(1-exp(-2*a*delT))+(2*dt/a))*(sigma*sigma/(4*a));
      return valToReturn;
    }*/
    Price EuroDollarFutureAdjustment(Rate r, FutureTime t0, Tenor delta, Maturity t){ //this is based on continuously compouned rates! this is an approximation; even assuming a hull-white model
      return (sigma*sigma*.25/(delta*a*a*a))*(1-exp(-2*a*(t-t0))-exp(-2*a*delta)+exp(-2*a*(t+delta-t0)));
    }

    /*Price EuroDollarFutureAdjustment(Rate r, FutureTime t0, Tenor delta, Maturity t){
      //0<t<t0<t1<t2
      double t1=t;
      double t2=t+delta;
      t=t0;
      t0=t1;


      double coef1=sigma*sigma/(2*a*a*a);
      return (Bond_Price(r, 0, t1)/Bond_Price(r, 0, t2))*(exp(coef1*(exp(-a*t1)-exp(-a*t2))*(exp(a*t0)-1)*(2-exp(-a*(t2-t0)-exp(-a*t2))))-1)/delta;
    }*/

    void setFutureTimes(std::map<int, double> &endingTimes){
      if(!storeParameters.empty()){
        storeParameters.erase(storeParameters.begin()); //time->{drift->dr, vol->v}
      }
      double t=0;
      //int i=0;
      //double r=r0;
      auto explicitPath=[&](double t1, double t2, int key){ //
        //double nextR=0;
        double expT=exp(-a*(t2-t1));
      //  double expT1=1-exp(-a*t1);//only used under risk neutral drift
        //double expT2=1-exp(-a*t2); //only used under risk neutral drift

        std::vector<double> driftVol(3);
        //driftVol[0]=yieldClass.Forward(t2)-yieldClass.Forward(t1)*expT+((sigma*sigma)/(2*a*a))*(expT2*expT2-expT*expT1*expT1);; //drift under risk neutral
        driftVol[0]=b*(1-expT); //drift
        driftVol[1]=expT; //damp
        driftVol[2]=sigma*sqrt((1-exp(-2*a*(t2-t1)))/(2*a)); //volatility
        storeParameters[key]=driftVol;
      };
      for(auto& iter : endingTimes) { //
        explicitPath(t,  iter.second, iter.first);
        t=iter.second;

      }

    }

    std::unordered_map<int, double> simulate(SimulNorm &nextNorm){ //pasthrough of some random generating thing
      std::unordered_map<int, double> simul;
      double r=r0;
      for(auto& iter : storeParameters) { //
        r=storeParameters[iter.first][0]+storeParameters[iter.first][1]*r+storeParameters[iter.first][2]*nextNorm.getNorm();//I hope SimulNorm is thread safe...
        simul[iter.first]=r;
      }
      return simul;
    }
    Price Swap_Price(Rate r, FutureTime t, SwapRate k, Tenor delta, Tenor freq, SwapMaturity mat){ //freq is the frequency of the bond payments...typically this is the same as delta, but not always.
      int numDates=(int)((mat-t)/freq)+1;//
      double firstDate=mat-(numDates-1)*freq;
      std::vector<double> bonds1(numDates);
      std::vector<double> bonds2(numDates);
      for(int i=0; i<numDates; i++){
        bonds1[i]=Bond_Price(r, t, firstDate+i*freq);
        bonds2[i]=Bond_Price(r, t, firstDate+i*freq+delta);
      }
      return getSwapPrice(bonds1, bonds2, k, delta);
    }
    Price Swap_Price(Rate r, FutureTime t, SwapRate k, Tenor delta, SwapMaturity mat){
      return Swap_Price(r, t, k, delta, delta, mat);
    }
    Price Swap_Price(SwapRate k, Tenor delta, SwapMaturity mat){
      return Swap_Price(r0, 0, k, delta, delta, mat);
    }
    Price Swap_Price(SwapRate k, Tenor delta, Tenor freq, SwapMaturity mat){
      return Swap_Price(r0, 0, k, delta, freq, mat);
    }
    SwapRate Swap_Rate(Tenor delta, SwapMaturity mat){
      return Swap_Rate(r0, 0, delta, delta, mat);
    }
    SwapRate Swap_Rate(Tenor delta, Tenor freq, SwapMaturity mat){
      return Swap_Rate(r0, 0, delta, freq, mat);
    }
    SwapRate Swap_Rate(Rate r, FutureTime t, Tenor delta, SwapMaturity mat){
      return Swap_Rate(r, t, delta, delta, mat);
    }
    SwapRate Swap_Rate(Rate r, FutureTime t, Tenor delta, Tenor freq, SwapMaturity mat){ //see http://www.columbia.edu/~mh2078/market_models.pdf
      int numDates=(int)((mat-t)/freq)+1; //this has been tested and works as it should
      double firstDate=mat-(numDates-1)*freq; //this should always be next payment time..so if t is a payment date, this should be t; else this should be the next payment.  This has been tested and works as it should
      std::vector<double> bonds1(numDates);
      std::vector<double> bonds2(numDates);
      for(int i=0; i<numDates; i++){
        bonds1[i]=Bond_Price(r, t, firstDate+i*freq); //this algorithm may be getting expensive...
        bonds2[i]=Bond_Price(r, t, firstDate+i*freq+delta);
      }
      return getSwapRate(bonds1, bonds2, delta);
    }


    //Price Swaption(Rate r, FutureTime tFut, Strike k, Tenor delta, SwapMaturity mat, OptionMaturity optMat){
      //return Swaption(r, tFut, k, delta, delta, mat, optMat);
    //}
    Price Swaption(Strike k, Tenor delta, SwapMaturity mat, OptionMaturity optMat){
      return Swaption(r0, 0, k, delta, mat, optMat);
    }
    Price Swaption(Rate r, FutureTime tFut, Strike k, Tenor delta, SwapMaturity mat, OptionMaturity optMat){
        int n=(int)(mat/delta); //swap maturity is the length of the swap at option expiration
        std::vector<double> cashFlowFimes(n);
        for(int i=0; i<n; ++i){
          cashFlowFimes[i]=optMat+(i+1)*delta;
        }
        return Swaption(r, tFut, k, cashFlowFimes, optMat);
    }
    Price Swaption(Rate r, FutureTime tFut, Strike k, std::vector<double>& cashFlowFimes, OptionMaturity optMat){
      Newton nt;
      int n=cashFlowFimes.size();
      std::vector<double> c(n);
      c[0]=k*(cashFlowFimes[0]-optMat);
      for(int i=1; i<n; i++){
        c[i]=k*(cashFlowFimes[i]-cashFlowFimes[i-1]);
      }
      c[n-1]=c[n-1]+1;
      double guess=r;
      nt.zeros([&](double r){
        double retVal=0;
        for(int i=0; i<n; i++){
          retVal+=c[i]*Bond_Price(r, optMat, cashFlowFimes[i]);
        }
        return retVal-1;
      },[&](double r){
        double retVal=0;
        for(int i=0; i<n; i++){
          retVal+=c[i]*Bond_Price(r, optMat, cashFlowFimes[i])*((exp(-a*(cashFlowFimes[i]-optMat))-1)/a);
        }

        return retVal;
      },guess);
      double retVal=0;
      //std::cout<<"tfut: "<<tFut<<" optMat: "<<optMat<<std::endl;
    //  std::cout<<cashFlows[0]<<std::endl;
      //std::cout<<Bond_Put(r, tFut, 1, 1, optMat)<<std::endl;
      for(int i=0; i<n; i++){
        /*if(isnan(c[i]*Bond_Put(r, tFut, Bond_Price(guess, optMat, cashFlows[i]), cashFlows[i], optMat))){
          std::cout<<"i: "<<i<<" tfut: "<<tFut<<" guess: "<<guess<<" optMat: "<<optMat<<"cashFlowi: "<<cashFlows[i]<<std::endl;
        }*/
        retVal+=c[i]*Bond_Put(r, tFut, Bond_Price(guess, optMat, cashFlowFimes[i]), cashFlowFimes[i], optMat);
        //std::cout<<"retVal: "<<retVal<<" i: "<<i<<std::endl;
      }
      //std::cout<<retVal<<std::endl;
      return retVal;
    }



};



#endif
