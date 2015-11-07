#include "YieldPolynomial.h"
YieldPolynomial::YieldPolynomial(YieldCurve& yield_){
  yield=yield_;
  computeYieldFunction();
}
YieldPolynomial::YieldPolynomial(){
}
void YieldPolynomial::computeYieldFunction(){
  computeYieldFunction(yield.size());
}
void YieldPolynomial::computeYieldFunction(int n){ //solves for polynomial that perfectly fits yield curve.  Note that this method may have computational difficulties
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

double YieldPolynomial::Yield(double t){ //y(0, t)*t
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
double YieldPolynomial::Forward(double t){ //f(0, t)
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
