#ifndef __BlackScholes_H_INCLUDED__
#define __BlackScholes_H_INCLUDED__
#include <cmath>

typedef double Underlying;
typedef double BSSigma;
typedef double Discount;
typedef double Strike;
typedef double OptionMaturity;
//typedef double Tenor;

double BSPut(Underlying, Strike, Discount, BSSigma, OptionMaturity);
double BSCall(Underlying, Strike, Discount, BSSigma, OptionMaturity);

#endif
