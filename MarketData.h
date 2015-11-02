#ifndef __MARKETDATA_H_INCLUDED__
#define __MARKETDATA_H_INCLUDED__
#include "Date.h"
#define NelsonSiegel 0
#define Polynomial 1

 struct AssetFeatures{
   Date Maturity;
   Date UnderlyingMaturity;
   double Strike;
   double Tenor;
   std::string type;
   AssetFeatures(const Date &maturity, double strike, const char type_[]){
     Maturity=maturity;
     Strike=strike;
     type=std::string(type_);
   }
   AssetFeatures(const Date &maturity, const char type_[]){
     Maturity=maturity;
     type=std::string(type_);
   }
   AssetFeatures(const Date &maturity, double strike, double tenor, const char type_[]){
     Maturity=maturity;
     Strike=strike;
     Tenor=tenor;
     type=std::string(type_);
   }
   AssetFeatures(const Date &maturity, double strike, double tenor, const Date &underlyingMaturity, const char type_[]){
     Maturity=maturity;
     Strike=strike;
     UnderlyingMaturity=underlyingMaturity;
     Tenor=tenor;
     type=std::string(type_);
   }
   AssetFeatures(const Date &maturity, double strike, const Date &underlyingMaturity, const char type_[]){
     Maturity=maturity;
     UnderlyingMaturity=underlyingMaturity;
     Strike=strike;
     type=std::string(type_);
   }
 };
 struct ForwardValue{
   Date beginDate;
   Date endDate;
   double value;
   ForwardValue(const Date &dt1, const Date &dt2, double val){
     beginDate=dt1;
     endDate=dt2;
     value=val;
   }
 };

 struct SpotValue{
   Date date;
   double value;
   SpotValue(const Date &dt, double val){
     date=dt;
     value=val;
   }
 };
#endif
