#ifndef __YIELDSPLINE_H_INCLUDED__
#define __YIELDSPLINE_H_INCLUDED__
#include <cmath>
#include "MarketData.h"
#include "Spline.h"
#include "Vasicek.h"
#include "BondUtilities.h"
typedef std::vector<SpotValue> YieldCurve;
typedef std::vector<SpotValue> FutureCurve;
typedef std::vector<SpotValue> LiborCurve;
typedef std::vector<SpotValue> SwapCurve;
class YieldSpline {
private:
  YieldCurve yield;
  FutureCurve forwardCurve;
  std::vector<double> theta;
  double maxMaturity;
  std::vector<double> splineX;
  std::vector<double> splineY;
  std::vector<double> splineZ;
  double r0;
  //char type[];
public:
  YieldSpline(YieldCurve&);
  YieldSpline();
  void computeSimpleSwapSpline(LiborCurve&, SwapCurve&);
  void computeYieldFunction();
  void computeSimpleFutureSpline(double, FutureCurve&, Vasicek<YieldSpline>&);

  double Yield(double);//Cubic Spline

  double Forward(double);//Cubic Spline

};

#endif
