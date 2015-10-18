#ifndef __VASICEK_H_INCLUDED__
#define __VASICEK_H_INCLUDED__
#include <vector>
#include <cmath>
#include "BlackScholes.h"
#include "MarketData.h"

typedef double Tenor;
typedef double Mu;
typedef double Speed;//mean reversion
typedef double ShortRateSigma;
typedef double BondMaturity;
typedef double Rate;
typedef double FutureTime;
typedef std::vector<double> Times;
typedef std::vector<double> Theta;
typedef std::vector<SpotValue> YieldCurve; //SpotValue is defined in "MarketData" as {Date, double}
class Vasicek{
  private:
    YieldCurve yield;
    Rate r0;
    Speed a;
    ShortRateSigma sigma;
    //BondMaturity t;
    Times times;
    Theta theta;
    int n;
    void estimateTheta();
    Date initial_t;
    double CtT(double, double);
  public:
  //  Vasicek(YieldCurve&); //all other parameters free
    Vasicek(YieldCurve&, Speed, ShortRateSigma, Rate); //fully specified
    //Vasicek(YieldCurve&, Speed, ShortRateSigma); //in case want to simulate bonds of various maturity and rate
    Discount Bond_Price(Rate, FutureTime, BondMaturity);
    Discount Bond_Price(BondMaturity);
    //Discount Bond_Price();
    //Discount Bond_Price();

};
//vanilla implementations
Discount Vasicek_Price(Rate, Speed, Mu, ShortRateSigma, BondMaturity);
//Discount Vasicek_Price(Rate, Speed, YieldCurve&, ShortRateSigma, BondMaturity);

BSSigma Vasicek_Volatility(Speed, ShortRateSigma, BondMaturity, OptionMaturity); //bond volatility, can be used in Black Scholes option pricing...note that this is a seperate function since we may need to estimate Speed and ShortRateSigma from this.
double Vasicek_Put(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Put(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Call(Rate, Speed, Mu, ShortRateSigma, Strike, BondMaturity, OptionMaturity);
double Vasicek_Caplet(Rate, Speed, Mu, ShortRateSigma, Strike, Tenor, OptionMaturity);
double Vasicek_Caplet(Underlying, Discount, Rate, Speed, ShortRateSigma, Strike, Tenor, OptionMaturity); //underlying here is B(t, T+\delta), discount is B(t, T)

/*class Bond {
 private:
   double a;
   double b;
   double sigma;
   double t;
   std::string type;
   double vol=0;
 public:
   Bond(double, double, double); //default to vasicek
   Bond(double, double, double, double); //default to vasicek, with t
   Bond(double, double, double, char const[]); //type "CIR" or "VASICEK"
   Bond(double, double, double, double, char const[]); //type "CIR" or "VASICEK"
   double getPrice(double);
   double getPrice(double, double);
   double getCIRPrice(double);
   double getCIRPrice(double, double);
   double getVasicekPrice(double);
   double getVasicekPrice(double, double);
   double getVasicekVolatility(double);
};*/

#endif
