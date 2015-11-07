#ifndef __YIELDNELSONSIEGAL_H_INCLUDED__
#define __YIELDNELSONSIEGAL_H_INCLUDED__
#include <cmath>
#include "Newton.h"
#include "MarketData.h"

typedef std::vector<SpotValue> YieldCurve;
typedef std::vector<double> Theta;
class YieldNelsonSiegal {
private:
  YieldCurve yield;
  std::vector<double> theta;
  double maxMaturity;

public:
  YieldNelsonSiegal(YieldCurve&);
  YieldNelsonSiegal();

  void computeYieldFunction();

  double Yield(double);
  double Forward(double);

};
#endif
