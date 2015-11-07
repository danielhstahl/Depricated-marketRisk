#ifndef __YIELDPOLYNOMIAL_H_INCLUDED__
#define __YIELDPOLYNOMIAL_H_INCLUDED__
#include <cmath>
#include <Eigen/Dense>
#include "MarketData.h"
typedef std::vector<SpotValue> YieldCurve;
typedef std::vector<double> Theta;
class YieldPolynomial {
private:
  YieldCurve yield;
  std::vector<double> theta;
  double maxMaturity;

  //char type[];
public:
  YieldPolynomial(YieldCurve&);
  YieldPolynomial();
  void computeYieldFunction(int);
  void computeYieldFunction();
  double Yield(double);
  double Forward(double);
};
#endif
