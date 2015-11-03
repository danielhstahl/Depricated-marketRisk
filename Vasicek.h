#ifndef __VASICEK_H_INCLUDED__
#define __VASICEK_H_INCLUDED__
#include <vector>
#include <cmath>
#include "BlackScholes.h"
#include "MarketData.h"
#include "Newton.h"
#include "SimulNorm.h"
#include "Swap.h"
#include <map>
#include <unordered_map>
#include <algorithm>
#include <Eigen/Dense>
#include "datatable.h"
#include "bsplineapproximant.h"

typedef double Tenor;
typedef double Price;
typedef double Mu;
typedef double Speed;//mean reversion
typedef double ShortRateSigma;
typedef double BondMaturity;
typedef double SwapMaturity;
typedef double Rate;
typedef double SwapRate;
typedef double FutureTime;
typedef double Coupon;
typedef std::vector<double> Times;
typedef std::vector<double> Theta;
typedef std::vector<SpotValue> YieldCurve; //SpotValue is defined in "MarketData" as {Date, double}
typedef std::vector<ForwardValue> VolatilityCurve; //ForwardValue is defined in "MarketData" as {Date, Date, double}
class Vasicek{
  private:
    YieldCurve yield;
    VolatilityCurve vCurve;
    Rate r0;
    Speed a;
    ShortRateSigma sigma;
    //int maxInt=0;
    //BondMaturity t;
    Times times;
    std::unordered_map<double, std::vector<double> > storeParameters;
    Theta theta;
    double minDiffT;
    int n;
    SPLINTER::BSplineApproximant* bspline;
    double thetaCalcPolynomial(double);
    //void estimateTheta();
    void estimateSpeedVolatility();
    void createContinuousYield();
    void createNSS();
    void BSpline();

    void estimateSpeedVolatility(double, double);
    Date initial_t;
    double maxMaturity=0; //this is the "largest" maturity available
    double thetaCalc(double);
    //double CtT(double, double);
  //std::unordered_map<int, std::unordered_map<int, double> > ctTMap;
    //std::unordered_map<int, double> diffusionMap;
  public:
  //  Vasicek(YieldCurve&); //all other parameters free
    Vasicek(YieldCurve&, Speed, ShortRateSigma, Rate); //fully specified
    Vasicek(YieldCurve&, VolatilityCurve&, Rate); //fully specified
    //Vasicek(YieldCurve&, Speed, ShortRateSigma); //in case want to simulate bonds of various maturity and rate
    double get_yield(double);
    double get_forward_rate(double);
    double get_forward_rate_polynomial(double);
    double get_yield_polynomial(double);
    double get_yield_spline(double);
    double get_forward_rate_spline(double);
    void deletePointers();
    //open question: should these "double" times be changed to "Date" times?  and let the class handle switching?
    Discount Bond_Price(Rate, FutureTime, BondMaturity);
    Discount Bond_Price(BondMaturity);
    Discount Bond_Price(Rate, FutureTime, Coupon, std::vector<double>&);

    Price Bond_Put(Rate, FutureTime, Strike, BondMaturity, OptionMaturity) ;
    Price Bond_Put(Strike, BondMaturity, OptionMaturity) ;
    Price EuroDollarFuture(Rate, FutureTime, Tenor, Maturity);
    Price Forward(Rate, FutureTime, Tenor, Maturity);
    Price Bond_Call(Rate, FutureTime, Strike, BondMaturity, OptionMaturity) ;
    Price Bond_Call(Strike, BondMaturity, OptionMaturity) ;
    Price Bond_Call(Rate, FutureTime, Strike, Coupon, std::vector<double>&, OptionMaturity);
    Price Bond_Put(Rate, FutureTime, Strike, Coupon, std::vector<double>&, OptionMaturity);
    Price Caplet(Rate, FutureTime, Strike, Tenor, OptionMaturity) ;
    Price Caplet(Strike, Tenor, OptionMaturity);
    std::unordered_map<double, double> simulate(SimulNorm&);
    void setFutureTimes(std::map<double, double>&);
    Price Swap_Price(Rate, FutureTime, SwapRate, Tenor, Tenor, SwapMaturity);
    Price Swap_Price(Rate, FutureTime, SwapRate, Tenor, SwapMaturity);
    Price Swap_Price(SwapRate, Tenor, SwapMaturity);
    Price Swap_Price(SwapRate, Tenor, Tenor, SwapMaturity);
    SwapRate Swap_Rate(Tenor, SwapMaturity);
    SwapRate Swap_Rate(Rate, FutureTime, Tenor, SwapMaturity);
    SwapRate Swap_Rate(Rate, FutureTime, Tenor, Tenor, SwapMaturity);
    SwapRate Swap_Rate(Tenor, Tenor, SwapMaturity);
    //Price Swaption(Rate, FutureTime, Strike, Tenor, Tenor, SwapMaturity, OptionMaturity);
    Price Swaption(Rate, FutureTime, Strike, Tenor, SwapMaturity, OptionMaturity);
    Price Swaption(Strike, Tenor, SwapMaturity, OptionMaturity);
    Price Swaption(Rate, FutureTime, Strike, std::vector<double>&, OptionMaturity);
};
//vanilla implementations
Discount Vasicek_Price(Rate, Speed, Mu, ShortRateSigma, BondMaturity);
//Discount Vasicek_Price(Rate, Speed, YieldCurve&, ShortRateSigma, BondMaturity);

BSSigma Vasicek_Volatility(Speed, ShortRateSigma, BondMaturity, OptionMaturity); //bond volatility, can be used in Black Scholes option pricing...note that this is a seperate function since we may need to estimate Speed and ShortRateSigma from this.
double Vasicek_Put(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Put(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
//double Vasicek_Swaption(Rate, Speed, Mu, ShortRateSigma, Strike, std::vector<double>&);
double Vasicek_Caplet(Rate, Speed, Mu, ShortRateSigma, Strike, Tenor, OptionMaturity);
double Vasicek_Caplet(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, Tenor, OptionMaturity); //underlying here is B(t, T+\delta), discount is B(t, T)


#endif
