#ifndef __HandleYield_H_INCLUDED__
#define __HandleYield_H_INCLUDED__
#include <cmath>
#include "Newton.h"
#include <Eigen/Dense>
#include "MarketData.h"
#include <type_traits>
typedef std::vector<SpotValue> YieldCurve;

template <int yieldType=0>
class HandleYield {
private:
  YieldCurve yield;
  std::vector<double> theta;
  double maxMaturity;
  //char type[];
public:
  HandleYield(YieldCurve& yield_){
    yield=yield_;
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T),yieldType==0), void>::type
  computeYieldFunction(){//NelsonSiegel
    Newton nt;
    std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
    int n=yield.size();
    std::vector<double> dataToMinimizeOver(n); //the yield from above
    std::vector<std::vector<double> > additionalParameters(n, std::vector<double>(1)); //t
    for(int i=0; i<n; i++){ //populate our functions to minimize over.
      yield[i].date.setScale("year");
      dataToMinimizeOver[i]=yield[i].value;
      Date currDate;
      additionalParameters[i][0]=yield[i].date-currDate;
      std::cout<<"data to minimize: "<<dataToMinimizeOver[i]<<std::endl;
      meanSquare.push_back([](std::vector<double> &guess, std::vector<double> &additionalParameters){
        double tLambda=additionalParameters[0]*guess[3];
        double expLambda=exp(-tLambda);
        return guess[0]+guess[1]*(1-expLambda)/tLambda+guess[2]*((1-expLambda)/tLambda-expLambda);
      });
    }
    theta=std::vector<double>(4);
    theta[0]=.05; //these seem decent guesses
    theta[1]=.05;//these seem decent guesses
    theta[2]=.05;//these seem decent guesses
    theta[3]=.05;//these seem decent guesses
    nt.optimize(meanSquare, additionalParameters, dataToMinimizeOver, theta);
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T), yieldType==1), void>::type
  computeYieldFunction(){ //polynomial
    computeYieldFunction(0);
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T), yieldType==1), void>::type  //polynomial
  computeYieldFunction(int n){ //solves for polynomial that perfectly fits yield curve.  Note that this method may have computational difficulties
    //int n=yield.size();
    int m=yield.size();
    if(m<n || n<1){
      n=m;
    }
    Eigen::MatrixXd HoldParameters(n, n);
    Eigen::VectorXd ThetaEigen(n);
    Eigen::VectorXd yieldValues(n);
    Date currDate;
    int k=0;
    double tt=0;
    for(int i=0; i<n; ++i){ //n needs to be greater than 2..
      k=i*(m/n);
      HoldParameters(i, 0)=1;
      HoldParameters(i, 1)=yield[k].date-currDate;
      yieldValues(i)=yield[k].value*HoldParameters(i, 1);
      if(HoldParameters(i, 1)>maxMaturity){
        maxMaturity=HoldParameters(i, 1);
      }
      tt=HoldParameters(i, 1)*HoldParameters(i, 1);
      for(int j=2; j<n; ++j){
        HoldParameters(i, j)=tt;//pow(HoldParameters(i, 1), j);
        tt*=HoldParameters(i, 1);
      }
    }
    ThetaEigen=HoldParameters.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(yieldValues);
    theta=Theta(n);
    for(int i=0; i<n; ++i){
      theta[i]=ThetaEigen(i);
      std::cout<<"Theta: "<<theta[i]<<std::endl;
    }
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T), yieldType==0), double>::type //NelsonSiegel
  Yield(double t){
    double expLambda=exp(-t*theta[3]);
    return theta[0]*t+theta[1]*(1-expLambda)/theta[3]+theta[2]*((1-expLambda)/theta[3]-expLambda*t);
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T), yieldType==0), double>::type //NelsonSiegel
  Forward(double t){ //f(0, t)
    double tLambda=t*theta[3];
    double expLambda=exp(-tLambda);
    return theta[0]+theta[1]*expLambda+theta[2]*tLambda*expLambda;
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T), yieldType==1), double>::type  //Polynomial
  Yield(double t){ //y(0, t)*t
    int n=theta.size();
    double thetVal=0;
    thetVal+=theta[0]+theta[1]*t;
    double getPow=t*t;
    for(int i=2; i<n; ++i){
      thetVal+=theta[i]*getPow;
      getPow=getPow*t;
    }
    return thetVal;///t;
  }
  template<class T> //this feels super hacky
  typename std::enable_if<(sizeof(T), yieldType)==1, double>::type //Polynomial
  Forward(double t){ //f(0, t)
    int n=theta.size();
    double thetVal=0;
    thetVal+=theta[1];
    double getPow=t;
    //double yield=0;
    //yield+=theta[0]+theta[1]*t;

    for(int i=2; i<n; ++i){
      thetVal+=theta[i]*getPow*i;
      getPow=getPow*t;
      //yield+=theta[i]*getPow;
    }
    return thetVal;//*t+yield;
  }


};
#endif
