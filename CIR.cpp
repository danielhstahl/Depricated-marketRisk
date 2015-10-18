#include "CIR.h"
Discount CIR_Price(Rate r, Speed a, Mu b, ShortRateSigma sigma, BondMaturity t){
  double h=sqrt(a*a+2*sigma*sigma);
  double xp=(2*a*b)/(sigma*sigma);
  double den=2*h+(a+h)*(exp(t*h)-1);
  return pow(2*h*exp((a+h)*t*.5)/den, xp)*exp(-(2*exp(t*h)-1)*r/den);
}
