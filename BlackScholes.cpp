#include "BlackScholes.h"

double BSCall(Underlying S0, Strike k, Discount discount, BSSigma sigma, OptionMaturity t){
  double s=sqrt(2);
  double sqrtT=sqrt(t);
  double d1=log(S0/(discount*k))/(sigma*sqrtT)+sigma*sqrtT*.5;
  return S0*(.5+.5*erf(d1/s))-k*discount*(.5+.5*(erf((d1-sigma*sqrtT)/s)));
}

double BSPut(Underlying S0, Strike k, Discount discount, BSSigma sigma, OptionMaturity t){
  double s=sqrt(2);
  double sqrtT=sqrt(t);
  double d1=log(S0/(discount*k))/(sigma*sqrtT)+sigma*sqrtT*.5;
  return S0*(.5*erf(d1/s)-.5)+k*discount*(.5-.5*(erf((d1-sigma*sqrtT)/s)));
}
