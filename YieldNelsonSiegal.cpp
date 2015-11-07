#include "YieldNelsonSiegal.h"
YieldNelsonSiegal::YieldNelsonSiegal(YieldCurve& yield_){
  yield=yield_;
  computeYieldFunction();
}
YieldNelsonSiegal::YieldNelsonSiegal(){
}

void YieldNelsonSiegal::computeYieldFunction(){//NelsonSiegel
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

double YieldNelsonSiegal::Yield(double t){
  double expLambda=exp(-t*theta[3]);
  return theta[0]*t+theta[1]*(1-expLambda)/theta[3]+theta[2]*((1-expLambda)/theta[3]-expLambda*t);
}
double YieldNelsonSiegal::Forward(double t){ //f(0, t)
  double tLambda=t*theta[3];
  double expLambda=exp(-tLambda);
  return theta[0]+theta[1]*expLambda+theta[2]*tLambda*expLambda;
}
