#include "Swap.h"

Swap::Swap(const std::vector<double> &libor_, double tenor_){
  //maxT=horizon; //eg 30 (years)
  tenor=tenor_; //eg .5 (half a year)
  libor=libor_; //array of libor rates given as (b_i-b_{i+1})/(b_{i+1}*tenor)
  n=libor.size()+1;
  bonds=std::vector<double>(n);//b(t, t+tenor*i)
  bonds[0]=1;//b(t, t)
  for(int i=1; i<n; i++){
    bonds[i]=bonds[i-1]/(tenor*libor[i-1]+1.0);
    std::cout<<bonds[i]<<" "<<(std::pow(bonds[i], -1.0/(double)i)-1.0)/(tenor*i)<<std::endl;
  }
}
double Swap::getPrice(double timeToMaturity, double k){
  return getPrice(timeToMaturity, 1, k);
}
double Swap::getPrice(double timeToMaturity, double discountToNextPayment, double k){ //k is the swap rate (should be a constant for each swap).  discountToNextPayment is the (unannualized) discount from now to the first payment
  double khat=(1.0/tenor+k)*discountToNextPayment; //in case less than exactly "tenor" to next payment
  double tenorHat=tenor*discountToNextPayment;
  std::vector<double> bonds(n);
  double price=0;
  int numToMaturity=(int)(timeToMaturity/tenor);
  std::cout<<numToMaturity<<std::endl;
  for(int i=1; i<numToMaturity;i++){ //dont count first bond since the first payment occurs after the first tenor
    price=price+bonds[i]/tenorHat-bonds[i+1]*khat;
  }
  return price;
}
