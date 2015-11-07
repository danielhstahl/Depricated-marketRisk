#include "YieldSpline.h"

YieldSpline::YieldSpline(YieldCurve& yield_){
  yield=yield_;
  computeYieldFunction();
}
YieldSpline::YieldSpline(){
}

void YieldSpline::computeYieldFunction(){//Cubic Spline
  int n=yield.size();
  Date currDate;
  splineX=std::vector<double>(n+1);
  splineY=std::vector<double>(n+1);
  double dt=0;
  splineX[0]=0;
  splineY[0]=0;
  for(int i=0;i<n;++i){
    yield[i].date.setScale("year");
    dt=yield[i].date-currDate;
    splineX[i+1]=dt;
    splineY[i+1]=yield[i].value*dt;
  }
  splineZ=spline(splineX, splineY);
}

double YieldSpline::Yield(double t){//Cubic Spline
  return splint(splineX, splineY, splineZ, t);
}

double YieldSpline::Forward(double t){//Cubic Spline
  return splintD(splineX, splineY, splineZ, t);
}
