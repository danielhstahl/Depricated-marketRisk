#include <vector>
#include "Spline.h"
std::vector<double> spline(std::vector<double> &x, std::vector<double> &y)//, double yp1, double ypn)
/*Given arrays x[0..n] and y[0..n] containing a tabulated function, i.e., yi = f(xi), with
x1 < x2 <...< xN , and given values yp1 and ypn for the first derivative of the interpolating
function at points 0 and n, respectively, this routine returns an array y2[0..n] that contains
the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary.*/
  {
  int i,k,n;
  n=x.size();
  double p,qn,sig,un;
  std::vector<double> u(n-1);
  std::vector<double> y2(n);
  //u=vector(1,n-1);
  //if (yp1 > 0.99e30) {
    y2[0]=u[0]=0.0;
  //} //The lower boundary condition is set either to be “nat  ural”
  //else { //or else to have a specified first derivative.
  //  y2[0] = -0.5;
  //  u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
//  }
  for (i=1;i<n-1;i++) {/* This is the decomposition loop of the tridiagonal algorithm.
  y2 and u are used for temporary
  storage of the decomposed
  factors.*/
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  //if (ypn > 0.99e30) {//The upper boundary condition is set either to be“natural”
    qn=un=0.0;
  //}
  //else { //or else to have a specified first derivative.
  //  qn=0.5;
  //  un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  //}
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

  for (k=n-2;k>-1;k--){ //This is the backsubstitution loop of the tridiagonal algorithm.
    y2[k]=y2[k]*y2[k+1]+u[k];
  }
  return y2;
}
double splint(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, double x)
/*Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function (with the xai’s in order),
and given the array y2a[0..n-1], which is the output from spline above, and given a value of
x, this routine returns a cubic-spline interpolated value y.*/
{
  //void nrerror(char error_text[]);
  int klo,khi,k,n;
  n=xa.size();
  double h,b,a,y;
  klo=0; /*We will find the right place in the table by means of
  bisection. This is optimal if sequential calls to this
  routine are at random values of x. If sequential calls
  are in order, and closely spaced, one would do better
  to store previous values of klo and khi and test if
  they remain appropriate on the next call.*/
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) {
      khi=k;
    }
    else {
      klo=k;
    }
  } //klo and khi now bracket the input value of x.
  h=xa[khi]-xa[klo];
  //if (h == 0.0) {
  //  std::cout<<"Bad xa input to routine splint"<<std::endl; //The xa’s must be dis  tinct.
  //}
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h; //Cubic spline polynomial is now evaluated.
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return y;
}
double splintD(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, double x)
/*Given the arrays xa[0..n-1] and ya[0..n-1], which tabulate a function (with the xai’s in order),
and given the array y2a[0..n-1], which is the output from spline above, and given a value of
x, this routine returns a cubic-spline interpolated value y.*/
{
  //void nrerror(char error_text[]);
  int klo,khi,k,n;
  n=xa.size();
  double h,b,a,y;
  klo=0; /*We will find the right place in the table by means of
  bisection. This is optimal if sequential calls to this
  routine are at random values of x. If sequential calls
  are in order, and closely spaced, one would do better
  to store previous values of klo and khi and test if
  they remain appropriate on the next call.*/
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) {
      khi=k;
    }
    else {
      klo=k;
    }
  } //klo and khi now bracket the input value of x.
  h=xa[khi]-xa[klo];
  //if (h == 0.0) {
  //  std::cout<<"Bad xa input to routine splint"<<std::endl; //The xa’s must be dis  tinct.
  //}
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h; //Cubic spline polynomial is now evaluated.
  //y=a*ya[klo]+b*ya[khi]+(a*a*a*y2a[klo]-a*y2a[klo]+b*b*b*y2a[khi]-b*y2a[khi])*(h*h)/6.0;
  y=-ya[klo]/h+ya[khi]/h+((1.0-3*a*a)*y2a[klo]+(3*b*b-1.0)*y2a[khi])*h/6.0;//derivative
  return y;
}
