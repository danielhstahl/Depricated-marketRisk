#ifndef __COMPUTEPORTFOLIO_H_INCLUDED__
#define __COMPUTEPORTFOLIO_H_INCLUDED__
#include <vector>
#include <string>
#include "MarketData.h"
#include "Date.h"
#include <iostream>
#include <map>
typedef std::vector<AssetFeatures> Portfolio;

template <typename T>
class ComputePortfolio{
private:
  Portfolio portfolio;
  int n;
  T* PricingEngine;
  SimulNorm simul;
  double future;
  int intFuture;
public:
  ComputePortfolio(){

  }
  ComputePortfolio(T *engine, Portfolio &prt){
    PricingEngine=engine;
    portfolio=prt;
    n=portfolio.size();
  }
  ComputePortfolio(T *engine, Portfolio &prt, int daysToSimulate){
    PricingEngine=engine;
    portfolio=prt;
    n=portfolio.size();
    setDates(daysToSimulate);
  }
  /*~ComputePortfolio() {
      delete PricingEngine;
  }*/
  void setDates(int daysToSimulate){
    std::map<int, double> holdPossibleDates; //using "map" to ensure sorted keys
    double diff=0;
    int key=0;
    Date currDate;
    currDate.setScale("day");
    Date simulateToDate=currDate+daysToSimulate;
    future=simulateToDate-currDate;
    intFuture=simulateToDate.msDiff(currDate);
    for(int i=0; i<n;i++){
      if(portfolio[i].Maturity-simulateToDate<=0){//if any maturities are less than the simulated date
        diff=portfolio[i].Maturity-currDate;
        key=portfolio[i].Maturity.msDiff(currDate);
        holdPossibleDates[key]=diff;
      }
    }
    holdPossibleDates[intFuture]=future;
    PricingEngine->setFutureTimes(holdPossibleDates);
  }
  double execute(double r, double t){
    Date currDate;
    double portVal=0;
    double valueOfPortfolio=0;
    for(int i=0; i<n; ++i){
      if(portfolio[i].type=="bond"){ //concerned here...I want this to be more automated (eg, simply call bond price automatically, but how to deal with parameters?)
        portVal+=PricingEngine->Bond_Price(r, t, portfolio[i].Maturity-currDate);//zero coupons remember...but are we assuming reinvestment??  Good grief that could get complicated...
      }
      else if(portfolio[i].type=="caplet"){
        portVal+=PricingEngine->Caplet(r, t, portfolio[i].Strike, portfolio[i].Tenor, portfolio[i].Maturity-currDate);
      }
      else if(portfolio[i].type=="swaption"){
        portVal+=PricingEngine->Swaption(r, t, portfolio[i].Strike, portfolio[i].Tenor, portfolio[i].UnderlyingMaturity-currDate, portfolio[i].Maturity-currDate);
      }
      else if(portfolio[i].type=="call"){
        portVal+=PricingEngine->Bond_Call(r, t, portfolio[i].Strike, portfolio[i].UnderlyingMaturity-currDate, portfolio[i].Maturity-currDate);
      }
    }
    return portVal;
  }
  double execute(const std::unordered_map<int, double> &ratePath){
    Date currDate;
    double portVal=0;
    double r=0;
    double t=0;
    double diffDate=0;
    int key=0;
    for(int i=0; i<n; ++i){
      r=0;
      t=0;
      diffDate=portfolio[i].Maturity-currDate;
      key=portfolio[i].Maturity.msDiff(currDate);
      if(ratePath.find(key)!=ratePath.end()){ //if maturity is less than simulation
        r=ratePath.find(key)->second;
        t=diffDate;
      }
      else{
        r=ratePath.find(intFuture)->second;
        t=future;
      }
      //if(diffDate<t){
      //  std::cout<<"Something is wrong!: t:"<<t<<" diffDate: "<<diffDate<<std::endl;
      //}
      if(portfolio[i].type=="bond"){ //concerned here...I want this to be more automated (eg, simply call bond price automatically, but how to deal with parameters?)
        portVal+=PricingEngine->Bond_Price(r, t, diffDate);//zero coupons remember...but are we assuming reinvestment??  Good grief that could get complicated...
      }
      else if(portfolio[i].type=="caplet"){
        portVal+=PricingEngine->Caplet(r, t, portfolio[i].Strike, portfolio[i].Tenor, diffDate);
      }
      else if(portfolio[i].type=="swaption"){
        portVal+=PricingEngine->Swaption(r, t, portfolio[i].Strike, portfolio[i].Tenor, (portfolio[i].UnderlyingMaturity-currDate)-diffDate, diffDate);
      }
      else if(portfolio[i].type=="call"){
        portVal+=PricingEngine->Bond_Call(r, t, portfolio[i].Strike, portfolio[i].UnderlyingMaturity-currDate, diffDate);
      }
      /*if(std::isnan(portVal)){
        std::cout<<"r: "<<r<<" t: "<<t<<" type: "<<portfolio[i].type<<" diffdate: "<<diffDate<<" underlying: "<<portfolio[i].UnderlyingMaturity-currDate<<" date: "<<portfolio[i].UnderlyingMaturity<<std::endl;
        break;
      }*/
    }
  //  std::cout<<"test for bond: "<<PricingEngine->Swaption(r, t, .03, .25, 2.030556, .030556)<<" r: "<<r<<" t: "<<t<<std::endl;
    //std::cout<<"diffdate: "<<diffDate<<std::endl;
    //std::cout<<"portfolio value: "<<portVal<<"r: "<<r<<"t: "<<t<<std::endl;
    //std::cout<<"future: "<<intFuture<<" r: "<<std::endl;
    //std::cout<<"minR: "<<minR<<std::endl;
    return portVal;
  }


};
#endif
