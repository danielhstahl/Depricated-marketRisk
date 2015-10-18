#ifndef __CIR_H_INCLUDED__
#define __CIR_H_INCLUDED__
#include <vector>
#include <cmath>
#include <string>
#include "BlackScholes.h"
#include "Date.h"

typedef double Tenor;
typedef double Mu;
struct Yield {
  Date date;
  double yield;
};
typedef std::vector<Yield> YieldCurve;
typedef double Speed;//mean reversion
typedef double ShortRateSigma;
typedef double BondMaturity;
typedef double Rate;

Discount CIR_Price(Rate, Speed, Mu, ShortRateSigma, BondMaturity);



#endif
