#include "Swap.h"

/*Swap::Swap(const std::vector<double> &libor_, double tenor_){
  //maxT=horizon; //eg 30 (years)
  tenor=tenor_; //eg .5 (half a year)
  libor=libor_; //array of libor rates given as (b_i-b_{i+1})/(b_{i+1}*tenor)
  n=libor.size()+1;
  bonds=std::vector<double>(n);//b(t, t+tenor*i)
  bonds[0]=1;//b(t, t)
  for(int i=1; i<n; i++){
    bonds[i]=bonds[i-1]/(tenor*libor[i-1]+1.0);
    //std::cout<<bonds[i]<<" "<<(std::pow(bonds[i], -1.0/(double)i)-1.0)/(tenor*i)<<std::endl;
  }
}*/
/*Swap::Swap(const std::vector<SpotValue> &yield, double tenor_){
    tenor=tenor_; //eg .5 (half a year)
    n=yield.size()+1;
    //bonds=std::vector<ForwardValue>;//b(t, t+tenor*i)

    //bonds[0].value=1;//b(t, t)
    Date currDate;
    bonds.push_back(ForwardValue(currDate, currDate, 1));

    for(int i=1; i<n; i++){
      bonds.push_back(ForwardValue(yield[i-1].date, yield[i-1].date, exp(-yield[i-1].value*(yield[i-1].date-currDate)))); //I dont understand why this errors but the other constructor doesn't.....

    }
}*/
Swap::Swap(double tenor_){
    tenor=tenor_; //eg .5 (half a year)
}
Swap::Swap(const std::vector<ForwardValue> &libor_){
  //maxT=horizon; //eg 30 (years)
  //tenor=tenor_; //eg .5 (half a year)
  libor=libor_; //array of libor rates given as (b_i-b_{i+1})/(b_{i+1}*tenor)
  n=libor.size();
  //bonds=std::vector<ForwardValue>;//b(t, t+tenor*i)
  bonds.push_back(ForwardValue(libor[0].beginDate, libor[0].beginDate, 1));
//  bonds[0].value=1;//b(t, t)
//  bonds[0].beginDate=libor[0].beginDate;//b(t, t)
  //bonds[0].endDate=libor[0].beginDate;//b(t, t)
  for(int i=1; i<(n+1); i++){
    libor[i-1].endDate.setScale("year");
    bonds.push_back(ForwardValue(libor[i-1].beginDate, libor[i-1].endDate, bonds[i-1].value/((libor[i-1].endDate-libor[i-1].beginDate)*libor[i-1].value+1.0)));
    //bonds[i].value=bonds[i-1].value/((libor[i-1].endDate-libor[i-1].beginDate)*libor[i-1].value+1.0);
    //bonds[i].beginDate=libor[i-1].beginDate;//b(t, t)
  //  bonds[i].endDate=libor[i-1].endDate;//b(t, t)
    std::cout<<bonds[i].value<<" "<<(std::pow(bonds[i].value, -1.0/(double)i)-1.0)/(libor[i-1].endDate-libor[i-1].beginDate)<<std::endl;
  }
}
double Swap::getPrice(double timeToMaturity, double k){
  return getPrice(timeToMaturity, 1, k);
}
double Swap::getPrice(std::vector<double> &BondPrices, std::vector<double> &BondPricesDelta, double k){
  double sPrice=0;
  k=1.0/tenor+k;
  for(int i : BondPrices){
    sPrice+=BondPrices[i]/tenor-k*BondPricesDelta[i];
  }
  return sPrice;
}
double Swap::getRate(std::vector<double> &BondPrices, std::vector<double> &BondPricesDelta){ //gets swap rate
  double den=0;
  double num=0;
  for(int i : BondPrices){
    num+=BondPrices[i];
    den+=BondPricesDelta[i];
  }
  return (num/den-1.0)/tenor;
}
double Swap::getPrice(double timeToMaturity, double discountToNextPayment, double k){ //k is the swap rate (should be a constant for each swap).  discountToNextPayment is the (unannualized) discount from now to the first payment
  double khat=(1.0/tenor+k)*discountToNextPayment; //in case less than exactly "tenor" to next payment
  double tenorHat=tenor*discountToNextPayment;
  double price=0;
  int numToMaturity=(int)(timeToMaturity/tenor);
  for(int i=1; i<numToMaturity;i++){ //dont count first bond since the first payment occurs after the first tenor
    price=price+bonds[i].value/tenorHat-bonds[i+1].value*khat;
  }
  return price;
}
double Swap::getRate(double timeToMaturity){ //gets swap rate
  int numToMaturity=(int)(timeToMaturity/tenor);
  double den=0;
  double num=0;
  for(int i=1; i<numToMaturity; i++){
    den+=bonds[i+1].value;
    num+=bonds[i].value;
  }
  return (num/den-1.0)/tenor;
}
double Swap::getBondPrice(double timeToMaturity){
  int numToMaturity=(int)(timeToMaturity/tenor);
  return bonds[numToMaturity+1].value;
}
double Swap::getBondYield(double timeToMaturity){
  int numToMaturity=(int)(timeToMaturity/tenor);
  return (1.0/bonds[numToMaturity].value-1.0)/timeToMaturity; //simple interest, no compounding
}
double getSwapPrice(std::vector<double> &BondPrices, std::vector<double> &BondPricesDelta, double k, double tenor){
  double sPrice=0;
  k=1.0/tenor+k;
  for(int i : BondPrices){
    sPrice+=BondPrices[i]/tenor-k*BondPricesDelta[i];
  }
  return sPrice;
}
double getSwapRate(std::vector<double> &BondPrices, std::vector<double> &BondPricesDelta, double tenor){ //gets swap rate
  double den=0;
  double num=0;
  for(int i : BondPrices){
    num+=BondPrices[i];
    den+=BondPricesDelta[i];
  }
  return (num/den-1.0)/tenor;
}
