#include "Vasicek.h"

Discount Vasicek_Price(Rate r, Speed a, Mu b, ShortRateSigma sigma, BondMaturity t){
  double at=(1-exp(-a*t))/a;
  double ct=sigma*sigma;
  ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
  return exp(-at*r+ct);
}

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
double Vasicek_Call(Underlying S0, Discount discount, Rate r, Speed a, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  return BSCall(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
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
