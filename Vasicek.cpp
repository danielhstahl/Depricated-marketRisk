#include "Vasicek.h"

/*Vasicek::Vasicek(YieldCurve& yield, Rate r0){


}*/
Vasicek::Vasicek(YieldCurve &yield_, Speed a_, ShortRateSigma sigma_, Rate r0_){
  a=a_;
  yield=yield_;
  sigma=sigma_;
  r0=r0_;
  estimateTheta();
}
void Vasicek::estimateTheta(){ //can only be run after Speed, ShortRateSigma are found
  n=yield.size();
  theta=Theta(n);
  times=Times(n+1);
  double at=0;
  double coef1=0;
  double coef2=0;
  double cf=sigma*sigma/(2*a*a);
  double plac=0;
  times[0]=0;
  for(int i=0; i<n; i++){
    plac=0;
    yield[i].date.setScale("year");
    times[i+1]=yield[i].date-initial_t;
    at=(1-exp(-a*times[i+1]))/a;
    coef1=(exp(-a*times[i+1])/a)*(exp(a*times[i+1])-exp(a*times[i]))-(times[i+1]-times[i]);
    coef2=cf*(a*.5*at*at-times[i+1]+at);
    for(int j=0; j<i; j++){
      plac+=theta[j]*((exp(a*times[j+1])-exp(a*times[j]))*(exp(-a*times[i+1])/a)-(times[j+1]-times[j]));
    }
    theta[i]=(at*r0-yield[i].value*times[i+1]+coef2-plac)/coef1;
    std::cout<<theta[i]<<std::endl;
  }
}
double Vasicek::CtT(double t1, double t2){ //computes C(t, T) based off given t and T
  int j=0;
  int i=0;
  double at=(1-exp(-a*t2))/a;
  double bprice=-(sigma*sigma/(2*a*a))*(a*.5*at*at-t2+t1+at);
  while(t1<times[i]){ //linear search is still efficient for small n;
      i++;
  }
  bprice+=theta[i-1]*((exp(a*times[i])-exp(a*t1))*exp(-a*(t2-t1))/(a)-(times[i]-t1));
  while(t2>times[i]){
    bprice+=theta[i]*((exp(a*times[i+1])-exp(a*times[i]))*exp(-a*(t2-t1))/(a)-(times[i+1]-times[i]));
    i++;
  }
  bprice+=theta[i]*((exp(a*t2)-exp(a*times[i]))*exp(-a*(t2-t1))/(a)-(t2-times[i]));
  return bprice;
}
Discount Vasicek::Bond_Price(Rate r, FutureTime t1, BondMaturity t2){
  return exp(-r*(1-exp(-a*(t2-t1)))/a+CtT(t1, t2));
}
Discount Vasicek::Bond_Price(BondMaturity t2){
  return exp(-r0*(1-exp(-a*t2))/a+CtT(0, t2));
}
Discount Vasicek_Price(Rate r, Speed a, Mu b, ShortRateSigma sigma, BondMaturity t){
  double at=(1-exp(-a*t))/a;
  double ct=sigma*sigma;
  ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
  return exp(-at*r+ct);
}
/*Discount Vasicek_Price(Rate r, Speed a, YieldCurve b, ShortRateSigma sigma, BondMaturity t){
  double at=(1-exp(-a*t))/a;
  double ct=sigma*sigma;
  ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
  return exp(-at*r+ct);
}*/
BSSigma Vasicek_Volatility (Speed a, ShortRateSigma sigma, BondMaturity tM, OptionMaturity t){ //can plug this into Black Scholes for option on bond
  return (1-exp(-a*(tM-t)))*(sigma/a)*sqrt((1-exp(-2*a*t))/(2*a*t));
}
double Vasicek_Put(Rate r, Speed a, Mu b, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  Underlying S0=Vasicek_Price(r, a, b, sigma, tM);
  Discount discount=Vasicek_Price(r, a, b, sigma, t);
  return BSPut(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}

double Vasicek_Put(Underlying S0, Discount discount, Rate r, Speed a, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  return BSPut(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}

double Vasicek_Call(Rate r, Speed a, Mu b, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  Underlying S0=Vasicek_Price(r, a, b, sigma, tM);
  Discount discount=Vasicek_Price(r, a, b, sigma, t);
  return BSCall(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}
double Vasicek_Caplet(Rate r, Speed a, Mu b, ShortRateSigma sigma, Strike k, Tenor delta, OptionMaturity t){
  return (k*delta+1)*Vasicek_Put(r, a, b, sigma, 1/(delta*k+1), t+delta, t);
}

double Vasicek_Caplet(Underlying S0, Discount discount, Rate r, Speed a, ShortRateSigma sigma, Strike k, Tenor delta, OptionMaturity t){
  return (k*delta+1)*Vasicek_Put(S0, discount, r, a, sigma, 1/(delta*k+1), t+delta, t);
}
