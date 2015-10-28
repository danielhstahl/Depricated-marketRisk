#include "Vasicek.h"

Vasicek::Vasicek(YieldCurve& yield_, VolatilityCurve& vCurve_, Rate r0_){
  yield=yield_;
  vCurve=vCurve_;
  r0=r0_;
  estimateSpeedVolatility();
  createContinuousYield();
//estimateTheta();
}
Vasicek::Vasicek(YieldCurve &yield_, Speed a_, ShortRateSigma sigma_, Rate r0_){
  a=a_;
  yield=yield_;
  sigma=sigma_;
  r0=r0_;
  createContinuousYield();
  //estimateTheta();
}
void Vasicek::createContinuousYield(){ //solves for polynomial that perfectly fits yield curve.  Note that this method may have computational difficulties: consider cubic splines.
  int n=yield.size();
  Eigen::MatrixXd HoldParameters(n, n);
  Eigen::VectorXd ThetaEigen(n);
  Eigen::VectorXd yieldValues(n);
  Date currDate;

  for(int i=0; i<n; ++i){ //n needs to be greater than 2..
    yieldValues(i)=yield[i].value;
    HoldParameters(i, 0)=1;
    HoldParameters(i, 1)=yield[i].date-currDate;
    if(HoldParameters(i, 1)>maxMaturity){
      maxMaturity=HoldParameters(i, 1);
    }
    for(int j=2; j<n; ++j){
      HoldParameters(i, j)=pow(HoldParameters(i, 1), j);
    }
  }
  ThetaEigen=HoldParameters.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(yieldValues);
  theta=Theta(n);
  for(int i=0; i<n; ++i){
    theta[i]=ThetaEigen(i);
    std::cout<<"Theta: "<<theta[i]<<std::endl;
  }
}
double Vasicek::thetaCalc(double t){ //returns theta for given (arbitrary) t.  IS this even called from anywhere?
  int n=theta.size();
  double thetVal=0;
  double thetValD=0;
  thetVal+=theta[0]+theta[1]*t;
  thetValD+=theta[1];
  double getPow=0;
  for(int i=2; i<n; ++i){
    getPow=theta[i]*pow(t, i);
    thetVal+=getPow;
    thetValD+=i*getPow;
  }
  return thetValD/a+thetVal+((sigma*sigma)/(2*a*a))*(1-exp(-2*a*t));
}
double Vasicek::get_yield_function(double t){ //F(0, t)
  int n=theta.size();
  double thetVal=0;
  thetVal+=theta[0]+theta[1]*t;
  double getPow=0;
  for(int i=2; i<n; ++i){
    thetVal+=theta[i]*pow(t, i);
  }
  return thetVal;
}
void Vasicek::estimateSpeedVolatility(){
  estimateSpeedVolatility(.3, .1);
}
void Vasicek::estimateSpeedVolatility(double guessMu, double guessSigma){
  Newton nt;
  std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
  int n=vCurve.size();
  std::vector<double> dataToMinimizeOver(n); //the volatility from above
  std::vector<std::vector<double> > additionalParameters(n, std::vector<double>(2)); //a and sigma
  for(int i=0; i<n; i++){ //populate our functions to minimize over.
    vCurve[i].beginDate.setScale("year");
    vCurve[i].endDate.setScale("year");
    dataToMinimizeOver[i]=vCurve[i].value;
    Date currDate;
    additionalParameters[i][0]=vCurve[i].endDate-currDate;
    additionalParameters[i][1]=vCurve[i].beginDate-currDate;
    meanSquare.push_back([](std::vector<double> &guess, std::vector<double> &additionalParameters){
      return Vasicek_Volatility(guess[0], guess[1], additionalParameters[0], additionalParameters[1]);
    });
  }
  std::vector<double> guess(2);
  guess[0]=guessMu; //these seem decent guesses
  guess[1]=guessSigma;//these seem decent guesses
  nt.optimize(meanSquare, additionalParameters, dataToMinimizeOver, guess);
  a=guess[0];
  sigma=guess[1];
}
/*void Vasicek::estimateTheta(){ //can only be run after Speed, ShortRateSigma are found
  n=yield.size();
  theta=Theta(n);
  times=Times(n+1);
  //std::vector<int> tm(n+1);
  double at=0;
  double coef1=0;
  double coef2=0;
  double cf=sigma*sigma/(2*a*a);
  double plac=0;
  //int maxInt=0;
  minDiffT=10000;
  times[0]=0;
  //tm[0]=0;
  for(int i=0; i<n; i++){
    yield[i].date.setScale("year");
    times[i+1]=yield[i].date-initial_t;
    plac=times[i+1]-times[i];
    if(plac<minDiffT){
      minDiffT=plac;
    }
  }
  //maxInt=(int)(times[n]/minDiffT);

  for(int i=0; i<n; i++){
    //tm[i+1]=times[i+1]*maxInt;
    plac=0;
    at=(1-exp(-a*times[i+1]))/a;
    coef1=(exp(-a*times[i+1])/a)*(exp(a*times[i+1])-exp(a*times[i]))-(times[i+1]-times[i]);
    coef2=cf*(a*.5*at*at-times[i+1]+at);
    for(int j=0; j<i; j++){
      plac+=theta[j]*((exp(a*times[j+1])-exp(a*times[j]))*(exp(-a*times[i+1])/a)-(times[j+1]-times[j]));
    }
    theta[i]=(at*r0-yield[i].value*times[i+1]+coef2-plac)/coef1;
    std::cout<<"theta: "<<theta[i]<<"  times: "<<times[i+1]<<std::endl;
  }

}*/

/*double Vasicek::CtT(double t1, double t2){ //computes C(t, T) based off given t and T...there HAS to be a reasonable method for making this a hashmap (lookup table)
  int j=0;
  int i=0;
  double at=(1-exp(-a*t2))/a;
  double bprice=-(sigma*sigma/(2*a*a))*(a*.5*at*at-t2+t1+at);
  if(t1==0){ //
    i++;
  }
  while(t1>times[i]){ //linear search is still efficient for small n;
      i++;
  }
  if(t2>times[i]){
    bprice+=theta[i-1]*((exp(a*times[i])-exp(a*t1))*exp(-a*(t2-t1))/(a)-(times[i]-t1));
    while(t2>times[i]){
      bprice+=theta[i]*((exp(a*times[i+1])-exp(a*times[i]))*exp(-a*(t2-t1))/(a)-(times[i+1]-times[i]));
      i++;
    }
    bprice+=theta[i]*((exp(a*t2)-exp(a*times[i]))*exp(-a*(t2-t1))/(a)-(t2-times[i]));
  }
  else{
    bprice+=theta[i]*((exp(a*t2)-exp(a*t1))*exp(-a*(t2-t1))/(a)-(t2-t1));
  }

  return bprice;
}*/


Discount Vasicek::Bond_Price(Rate r, FutureTime t1, BondMaturity t2) {
  double ctT=(1-exp(-a*(t2-t1)))/a;
  double atT=-ctT*r;
  double yieldFunc1=get_yield_function(t1);
  double yieldFunc2=get_yield_function(t2);
  ctT=ctT*yieldFunc1-(sigma*sigma/(4*a))*ctT*ctT*(1-exp(-2*a*t1));
  return exp(atT+ctT)*exp(-yieldFunc2*t2)/exp(-yieldFunc1*t1); //see hull and white's seminal paper
}
Discount Vasicek::Bond_Price(BondMaturity t2){
  return exp(-get_yield_function(t2)*t2);
  //return exp(-r0*(1-exp(-a*t2))/a+CtT(0, t2));
}
Price Vasicek::Bond_Put(Rate r, FutureTime t0, Strike k, BondMaturity tM, OptionMaturity t){
  return Vasicek_Put(Bond_Price(r, t0, tM), Bond_Price(r, t0, t), r, a, sigma, k, tM-t0, t-t0);
}
Price Vasicek::Bond_Put(Strike k , BondMaturity tM, OptionMaturity t){
  return Vasicek_Put(Bond_Price(tM), Bond_Price(t), r0, a, sigma, k, tM, t);
}
Price Vasicek::Bond_Call(Rate r, FutureTime t0, Strike k, BondMaturity tM, OptionMaturity t){
  return Vasicek_Call(Bond_Price(r, t0, tM), Bond_Price(r, t0, t), r, a, sigma, k, tM-t0, t-t0);
}

Price Vasicek::Bond_Call(Strike k , BondMaturity tM, OptionMaturity t){
  return Vasicek_Call(Bond_Price(tM), Bond_Price(t), r0, a, sigma, k, tM, t);
}

Price Vasicek::Caplet(Rate r, FutureTime t0, Strike k , Tenor delta, OptionMaturity t){
  return Vasicek_Caplet(Bond_Price(r, t0, t+delta), Bond_Price(r, t0, t), r, a, sigma, k, delta, t-t0);
}
Price Vasicek::Caplet(Strike k , Tenor delta, OptionMaturity t){
  return Vasicek_Caplet(Bond_Price(t+delta), Bond_Price(t), r0, a, sigma, k, delta, t);
}
void Vasicek::setFutureTimes(std::map<double, double> &endingTimes){
  if(!storeParameters.empty()){
    storeParameters.erase(storeParameters.begin()); //time->{drift->dr, vol->v}
  }
  double t=0;
  //int i=0;
  double r=r0;
  auto explicitPath=[&](double r, double t1, double t2){ //
    //double nextR=0;
    double expT=exp(-a*(t2-t1));

    //std::vector<double> driftVol(3);
    r=get_yield_function(t2)-get_yield_function(t1)*expT+((sigma*sigma)/(2*a*a))*(1+exp(-2*a*t2)-expT-exp(a*(t2+t1)));//make sure this is right...


    /*double nextR=0;
    if(t1==0){ //
      i++;
    }
    while(t1>times[i]){ //linear search is still efficient for small n;
        i++;
    }
    if(t2>times[i]){
      nextR=thetaCalc(times[i])*(exp(a*times[i])-exp(a*t1));///(times[i]-t1);
      while(t2>times[i]){
        nextR+=theta[i]*(exp(a*times[i+1])-exp(a*times[i]));///(times[i+1]-times[i]);
        i++;
      }
      nextR+=theta[i]*(exp(a*t2)-exp(a*times[i]));///(t2-times[i]);
    }
    else{
      nextR+=theta[i]*(exp(a*t2)-exp(a*t1));///(t2-t1);
    }

    nextR=nextR*exp(-a*t2)*a;*/
    std::vector<double> driftVol(3);
    driftVol[0]=r; //drift
    driftVol[1]=expT; //damp
    driftVol[2]=sigma*sqrt((1-exp(-2*a*(t2-t1)))/(2*a)); //volatility
    storeParameters[t2]=driftVol;
  };
  for(auto& iter : endingTimes) { //
    explicitPath(r, t,  iter.second);
    t=iter.second;
    r=storeParameters[t][0];
    std::cout<<"r at t: "<<r<<", "<<t<<std::endl;
  }

}
std::unordered_map<double, double> Vasicek::simulate(SimulNorm &nextNorm){ //pasthrough of some random generating thing
  std::unordered_map<double, double> simul;
  double r=r0;
  /*auto explicitPath=[&](double r, double t){
    //std::cout<<storeParameters[t][2]<<std::endl;
    return storeParameters[t][0]+storeParameters[t][1]*r+storeParameters[t][2]*nextNorm.getNorm();//I hope SimulNorm is thread safe...
  };*/
  for(auto& iter : storeParameters) { //
    r=storeParameters[iter.first][0]+storeParameters[iter.first][1]*r+storeParameters[iter.first][2]*nextNorm.getNorm();//I hope SimulNorm is thread safe...
    simul[iter.first]=r;
  }
  return simul;
}

Price Vasicek::Swap_Price(Rate r, FutureTime t, SwapRate k, Tenor delta, Tenor freq, SwapMaturity mat){ //freq is the frequency of the bond payments...typically this is the same as delta, but not always.
  int numDates=(int)((mat-t)/freq);
  double firstDate=mat-numDates*freq;
  std::vector<double> bonds1(numDates);
  std::vector<double> bonds2(numDates);
  for(int i=0; i<numDates; i++){
    bonds1[i]=Bond_Price(r, t, firstDate+i*freq);
    bonds2[i]=Bond_Price(r, t, firstDate+i*freq+delta);
  }
  return getSwapPrice(bonds1, bonds2, k, delta);
}
Price Vasicek::Swap_Price(Rate r, FutureTime t, SwapRate k, Tenor delta, SwapMaturity mat){
  return Swap_Price(r, t, k, delta, delta, mat);
}
Price Vasicek::Swap_Price(SwapRate k, Tenor delta, SwapMaturity mat){
  return Swap_Price(r0, 0, k, delta, delta, mat);
}
Price Vasicek::Swap_Price(SwapRate k, Tenor delta, Tenor freq, SwapMaturity mat){
  return Swap_Price(r0, 0, k, delta, freq, mat);
}
SwapRate Vasicek::Swap_Rate(Tenor delta, SwapMaturity mat){
  return Swap_Rate(r0, 0, delta, delta, mat);
}
SwapRate Vasicek::Swap_Rate(Tenor delta, Tenor freq, SwapMaturity mat){
  return Swap_Rate(r0, 0, delta, freq, mat);
}
SwapRate Vasicek::Swap_Rate(Rate r, FutureTime t, Tenor delta, SwapMaturity mat){
  return Swap_Rate(r, t, delta, delta, mat);
}
SwapRate Vasicek::Swap_Rate(Rate r, FutureTime t, Tenor delta, Tenor freq, SwapMaturity mat){
  int numDates=(int)((mat-t)/freq);
  double firstDate=mat-numDates*freq; //this should always be t+delta
  std::vector<double> bonds1(numDates);
  std::vector<double> bonds2(numDates);
  for(int i=0; i<numDates; i++){
    bonds1[i]=Bond_Price(r, t, firstDate+i*freq); //this algorithm may be getting expensive...
    bonds2[i]=Bond_Price(r, t, firstDate+i*freq+delta);
  }
  return getSwapRate(bonds1, bonds2, delta);
}
Price Vasicek::Swaption(Rate r, FutureTime tFut, Strike k, Tenor delta, SwapMaturity mat, OptionMaturity optMat){
  return Swaption(r, tFut, k, delta, delta, mat, optMat);
}
Price Vasicek::Swaption(Rate r, FutureTime tFut, Strike k, Tenor delta, Tenor freq, SwapMaturity mat, OptionMaturity optMat){
  int i=0;
  auto getTreeDrift=[&](double t, double x){
    while(times[i]<(t+tFut)){ //increment by future time since "t" always starts at zero in tree.
      i++;
    }
    return a*(theta[i]-x)/sigma;
  };
  auto inv=[&](double x){
    return sigma*x;
  };
  auto sig=[&](double t, double x){
    return 0;
  };
  auto payoff=[&](double x){
    double swp=Swap_Rate(x, optMat, delta, freq, mat);//careful...this won't work with american options since "optMat" doesn't change (needs to be a function of t)
    if(swp>k){
      return swp-k;
    }
    else{
      return 0.0;
    }
  };
  return getSwaptionPrice(200, r/sigma, optMat-tFut, getTreeDrift, sig, inv, payoff);
}

Discount Vasicek_Price(Rate r, Speed a, Mu b, ShortRateSigma sigma, BondMaturity t){
  double at=(1-exp(-a*t))/a;
  double ct=sigma*sigma;
  ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
  return exp(-at*r+ct);
}
/*Discount Vasicek_Price(Rate r, Speed a, YieldCurve b, ShortRateSigma sigma, BondMaturity t){
  double at=(1-exp(-a*t))/a;
  double ct=sigma*sigma;
  ct=(b-ct/(2*a*a))*(at-t)-ct*at*at/(4*a);
  return exp(-at*r+ct);
}*/
BSSigma Vasicek_Volatility (Speed a, ShortRateSigma sigma, BondMaturity tM, OptionMaturity t){ //can plug this into Black Scholes for option on bond
  return (1-exp(-a*(tM-t)))*(sigma/a)*sqrt((1-exp(-2*a*t))/(2*a*t));
}

double Vasicek_Put(Rate r, Speed a, Mu b, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  Underlying S0=Vasicek_Price(r, a, b, sigma, tM);
  Discount discount=Vasicek_Price(r, a, b, sigma, t);
  return BSPut(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}

double Vasicek_Put(Underlying S0, Discount discount, Rate r, Speed a, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  return BSPut(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}
double Vasicek_Call(Underlying S0, Discount discount, Rate r, Speed a, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  return BSCall(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}
double Vasicek_Call(Rate r, Speed a, Mu b, ShortRateSigma sigma, Strike k, BondMaturity tM, OptionMaturity t){
  Underlying S0=Vasicek_Price(r, a, b, sigma, tM);
  Discount discount=Vasicek_Price(r, a, b, sigma, t);
  return BSCall(S0, k, discount, Vasicek_Volatility(a, sigma, tM, t), t);
}
double Vasicek_Caplet(Rate r, Speed a, Mu b, ShortRateSigma sigma, Strike k, Tenor delta, OptionMaturity t){
  return (k*delta+1)*Vasicek_Put(r, a, b, sigma, 1/(delta*k+1), t+delta, t);
}

double Vasicek_Caplet(Underlying S0, Discount discount, Rate r, Speed a, ShortRateSigma sigma, Strike k, Tenor delta, OptionMaturity t){
  return (k*delta+1)*Vasicek_Put(S0, discount, r, a, sigma, 1/(delta*k+1), t+delta, t);
}
