#ifndef __MARKETDATA_H_INCLUDED__
#define __MARKETDATA_H_INCLUDED__
#include "Date.h"


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
