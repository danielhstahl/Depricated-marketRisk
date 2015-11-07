#ifndef __YIELDSPLINE_H_INCLUDED__
#define __YIELDSPLINE_H_INCLUDED__
#include <cmath>
#include "MarketData.h"
//#include <type_traits>
#include "Spline.h"
typedef std::vector<SpotValue> YieldCurve;
class YieldSpline {
private:
  YieldCurve yield;
  std::vector<double> theta;
  double maxMaturity;
  std::vector<double> splineX;
  std::vector<double> splineY;
  std::vector<double> splineZ;
  //char type[];
public:
  YieldSpline(YieldCurve&);
  YieldSpline();

  void computeYieldFunction();

  double Yield(double);//Cubic Spline

  double Forward(double t);//Cubic Spline

};
#endif
