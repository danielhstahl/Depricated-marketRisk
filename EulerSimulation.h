#ifndef __EULERSIMULATION_H_INCLUDED__
#define __EULERSIMULATION_H_INCLUDED__
#include <iostream>
#include <cmath>
#include "SimulNorm.h"
template<typename numberType>
class EulerSimulation{ //simulates SDE of type dX=alpha(X, t)dt+sigma(X, t)dW
  private:
    int steps;
    int numSimulations;
    double finalResult;
    double t;
    double sqrtT;
    double dt;
    double distVol=0;
  public:
    EulerSimulation(int steps_, int numSimulations_, double t_){
      steps=steps_;
      t=t_;
      numSimulations=numSimulations_;
      dt=t/(double)steps;
      sqrtT=sqrt(dt);
      //norm=std::normal_distribution<double>(0.0, sqrtT);
    }
    double getResult(){
      return finalResult;
    }
    double getMCVol(){
      return distVol;
    }
    template< typename ALPHA, typename SIGMA, typename CALLBACK, typename CALLBACKTH>
    void simulate(numberType initialValue, ALPHA&& alpha,  SIGMA&& sigma,  CALLBACKTH&& callbackTH, CALLBACK&& callback){ //optional callback takes the current "t" and current value as an argument...for time homogenous solutions
      numberType currVal=initialValue;//might be a double or a matrix;
      SimulNorm<numberType> y(initialValue.size(), 0.0, sqrtT);
      //y.setN(initialValue);
      double currT;
      double valTH=0;
      finalResult=0;
      distVol=0;
      #pragma omp parallel//multithread using openmp
      {
        #pragma omp for //multithread using openmp
        for(int j=0; j<numSimulations; j++){
          currVal=initialValue;
          valTH=0;

          for(int i=0; i<steps; i++){
            currT=(double)i*dt;
            callbackTH(currT, currVal, valTH); //manipulates "valTH"...pass by reference!
            //std::cout<<"got here"<<std::endl;
            currVal=currVal+alpha(currT, currVal)*dt+sigma(currT, currVal)*y.getNorm(); //Euler discretization
          }
          valTH=callback(currVal, valTH);
          finalResult+=valTH;
          distVol+=valTH*valTH; //distVol measures the MC error
        }
      }
      valTH=(double)numSimulations;
      finalResult=finalResult/valTH;
      distVol=sqrt((distVol/valTH-finalResult*finalResult)/valTH);//distVol measures the MC error
    }
    template< typename ALPHA, typename SIGMA>
    void simulate(numberType initialValue, ALPHA&& alpha, SIGMA&& sigma){ //if no callbacks are defined, simply return the expected value of the SDE
      simulate(initialValue, alpha, sigma, [&](double a, numberType b, double c){}, [&](numberType a, double b){return a.l1norm();});
    }
    template< typename ALPHA, typename SIGMA, typename CALLBACK>
    void simulate(numberType initialValue, ALPHA&& alpha, SIGMA&& sigma, CALLBACK&& callback){ //if only one callback defined, ignore path dependent parts
      simulate(initialValue, alpha, sigma, [&](double a, numberType b, double c){}, callback);
    }
};
#endif
