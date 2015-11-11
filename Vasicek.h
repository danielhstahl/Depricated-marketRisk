#ifndef __VASICEK_H_INCLUDED__
#define __VASICEK_H_INCLUDED__
#include <vector>
#include <cmath>
#include "BlackScholes.h"
#include "MarketData.h"
#include "Newton.h"
#include "SimulNorm.h"
#include "Swap.h"
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
typedef std::vector<double> Times;

typedef std::vector<SpotValue> YieldCurve; //SpotValue is defined in "MarketData" as {Date, double}
typedef std::vector<ForwardValue> VolatilityCurve; //ForwardValue is defined in "MarketData" as {Date, Date, double}

//vanilla implementations
Discount Vasicek_Price(Rate, Speed, Mu, ShortRateSigma, BondMaturity);
//Discount Vasicek_Price(Rate, Speed, YieldCurve&, ShortRateSigma, BondMaturity);

BSSigma Vasicek_Volatility(Speed, ShortRateSigma, BondMaturity, OptionMaturity); //bond volatility, can be used in Black Scholes option pricing...note that this is a seperate function since we may need to estimate Speed and ShortRateSigma from this.
double Vasicek_Put(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Put(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
//double Vasicek_Swaption(Rate, Speed, Mu, ShortRateSigma, Strike, std::vector<double>&);
double Vasicek_Caplet(Rate, Speed, Mu, ShortRateSigma, Strike, Tenor, OptionMaturity);
double Vasicek_Caplet(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, Tenor, OptionMaturity); //underlying here is B(t, T+\delta), discount is B(t, T)

template <typename T>
class Vasicek{
  private:
    YieldCurve yield;
    VolatilityCurve vCurve;
    Rate r0;
    Speed a;
    Mu b;
    ShortRateSigma sigma;
    //int maxInt=0;
    //BondMaturity t;
    //HandleYield<yieldType> yieldClass;
    T yieldClass;
    Times times;
    std::unordered_map<int, std::vector<double> > storeParameters;
    //Theta theta;
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
      nt.optimize(meanSquare, additionalParameters, dataToMinimizeOver, guess);
      a=guess[0];
      sigma=guess[1];
    }

    //double CtT(double, double);
  //std::unordered_map<int, std::unordered_map<int, double> > ctTMap;
    //std::unordered_map<int, double> diffusionMap;
  public:
    Vasicek(){

    }
    Vasicek(T &yield_, VolatilityCurve &vCurve_, Rate r0_){
      vCurve=vCurve_;
      r0=r0_;
      estimateSpeedVolatility();
      yieldClass=yield_;

    }

    Vasicek(T &yield_, Speed a_, ShortRateSigma sigma_, Rate r0_){
      a=a_;
      //yield=yield_;
      sigma=sigma_;
      r0=r0_;
      yieldClass=yield_;
      //createContinuousYield();
      //estimateTheta();
    }

    Discount Bond_Price(Rate r, FutureTime t1, BondMaturity t2) {
      if(t2<t1){
        return 1;
      }
      double ctT=(1-exp(-a*(t2-t1)))/a;
      double atT=-ctT*r;
      double yieldFunc1=yieldClass.Yield(t1);
      double yieldFunc2=yieldClass.Yield(t2);
      ctT=ctT*yieldClass.Forward(t1)-(sigma*sigma/(4*a))*ctT*ctT*(1-exp(-2*a*t1));
      return exp(atT+ctT)*exp(-yieldFunc2+yieldFunc1); //see hull and white's seminal paper
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
      for(int i=0; i<(n-1); ++i){
        historicalRates[i+1].date.setScale("year");
        dt=historicalRates[i+1].date-historicalRates[i].date;
        b+=(historicalRates[i+1].value-exp(-dt)*historicalRates[i].value)/(1.0-exp(-dt));
      }
      b=b/(n-1);
      std::cout<<"This is b: "<<b<<std::endl;
    }
    Price Bond_Call(Rate r, FutureTime t0, Strike k, Coupon cp, std::vector<double>& cashFlows, OptionMaturity optMat){
      double xStrike=0;
      Newton nt;
      //std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
      int n=cashFlows.size();
      double guess=r;
      nt.zeros([&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlows[i]);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlows[n-1]);
        return retVal-k;
      },[&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlows[i])*((exp(-a*(cashFlows[i]-optMat))-1)/a);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlows[n-1])*((exp(-a*(cashFlows[n-1]-optMat))-1)/a);
        return retVal;
      },guess);
      double retVal=0;
      for(int i=0; i<(n-1); i++){
        retVal+=cp*Bond_Call(r, t0, Bond_Price(guess, optMat, cashFlows[i]), cashFlows[i], optMat);
        //std::cout<<"retVal: "<<retVal<<" i: "<<i<<std::endl;
      }
      retVal+=(1+cp)*Bond_Call(r, t0, Bond_Price(guess, optMat, cashFlows[n-1]), cashFlows[n-1], optMat);
      return retVal;
    }

    Price Bond_Put(Rate r, FutureTime t0, Strike k, Coupon cp, std::vector<double>& cashFlows, OptionMaturity optMat){
      double xStrike=0;
      Newton nt;
      //std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
      int n=cashFlows.size();
      double guess=r;
      nt.zeros([&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlows[i]);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlows[n-1]);
        return retVal-k;
      },[&](double r){
        double retVal=0;
        for(int i=0; i<(n-1); i++){
          retVal+=cp*Bond_Price(r, optMat, cashFlows[i])*((exp(-a*(cashFlows[i]-optMat))-1)/a);
        }
        retVal+=(1+cp)*Bond_Price(r, optMat, cashFlows[n-1])*((exp(-a*(cashFlows[n-1]-optMat))-1)/a);
        return retVal;
      },guess);
      double retVal=0;
      for(int i=0; i<(n-1); i++){
        retVal+=cp*Bond_Put(r, t0, Bond_Price(guess, optMat, cashFlows[i]), cashFlows[i], optMat);
        //std::cout<<"retVal: "<<retVal<<" i: "<<i<<std::endl;
      }
      retVal+=(1+cp)*Bond_Put(r, t0, Bond_Price(guess, optMat, cashFlows[n-1]), cashFlows[n-1], optMat);
      return retVal;
    }

    Price EuroDollarFuture(Rate r, FutureTime t0, Tenor delta, Maturity t){
      double dT=1-exp(-a*t);
      double dt=1-exp(-a*t0);
      double dtT=exp(-a*(t-t0));
      double muT=r*dtT+yieldClass.Forward(t)+(sigma*sigma*.5)*dT*dT-yieldClass.Forward(t0)*dtT-dtT*sigma*sigma*.5*dt*dt;
      double sigT=(sigma*sigma*.5/a)*(1-exp(-a*2*(t-t0)));
      dtT=1-exp(-a*delta);
      double atT=dtT/a;
      double ctT=exp(-a*(t+delta))+dT-1.0;
      double yieldFunc1=yieldClass.Yield(t);
      double yieldFunc2=yieldClass.Yield(t+delta);
      ctT=yieldFunc1-yieldFunc2+atT*yieldClass.Forward(t)-(sigma*sigma*.25/(a*a*a))*ctT*ctT*(exp(2*a*t)-1);
      return (exp(atT*muT+.5*sigT*atT*atT-ctT)-1)/delta;
    }

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
      int numDates=(int)((mat-t-.00001)/freq);//makes sure rounds down if (mat-t)/freq is an integer
      double firstDate=mat-numDates*freq;
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

    SwapRate Swap_Rate(Rate r, FutureTime t, Tenor delta, Tenor freq, SwapMaturity mat){
      int numDates=(int)((mat-t)/freq);
      double firstDate=mat-numDates*freq; //this should always be t+delta
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
        int n=(int)(mat/delta)-1; //swap maturity is the length of the swap at option expiration
        std::vector<double> cashFlows(n);
        for(int i=0; i<n; ++i){
          cashFlows[i]=optMat+(i+1)*delta;
        }
        return Swaption(r, tFut, k, cashFlows, optMat);
    }

    Price Swaption(Rate r, FutureTime tFut, Strike k, std::vector<double>& cashFlows, OptionMaturity optMat){
      //double xStrike=0;
      Newton nt;
      //std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
      int n=cashFlows.size();
      std::vector<double> c(n);
      c[0]=k*(cashFlows[0]-optMat);
      for(int i=1; i<n; i++){
        c[i]=k*(cashFlows[i]-cashFlows[i-1]);
      }
      c[n-1]=c[n-1]+1;
      double guess=r;
      nt.zeros([&](double r){
        double retVal=0;
        for(int i=0; i<n; i++){
          retVal+=c[i]*Bond_Price(r, optMat, cashFlows[i]);
        }
        return retVal-1;
      },[&](double r){
        double retVal=0;
        for(int i=0; i<n; i++){
          retVal+=c[i]*Bond_Price(r, optMat, cashFlows[i])*((exp(-a*(cashFlows[i]-optMat))-1)/a);
        }
        return retVal;
      },guess);
      double retVal=0;
      for(int i=0; i<n; i++){
        retVal+=c[i]*Bond_Put(r, tFut, Bond_Price(guess, optMat, cashFlows[i]), cashFlows[i], optMat);
        //std::cout<<"retVal: "<<retVal<<" i: "<<i<<std::endl;
      }
      return retVal;
    }



};



#endif
