#define _USE_MATH_DEFINES
#include <iostream>
//#include "Swap.h"
#include <ctime>
#include <vector>
#include <chrono> //for accurate multithreading time using std::chrono
//#include "EulerSimulation.h"
//#include "Matrix.h"
//#include "Double.h"
#include "Date.h"
#include "BlackScholes.h"
#include "Vasicek.h"
#include "MarketData.h"
#include "Newton.h"


/*struct SpotValue {
  Date date;
  double value;
}*/
/*struct ForwardValue{
  Date beginDate;
  Date endDate;
  Double value;
};
Matrix alpha(double t, Matrix &x){
  return .1*x;
}
Matrix sigma(double t, Matrix &x){

  return .3*x;
}
double callback(Matrix val, double valth){
  return val.l1norm();
}*/

//double meanSquare(std::vector<SpotValue> &volatility){

//}

int main(){
  double a=.4;
  double b=.05;
  double sig=.05;
  double t=1;
  double r0=.03;
  double delta=.25;

  /*Newton Test */
  Newton nt;
  int m=10;//number of volatility dates
  int k=10; //number of volatility strikes per date;
  std::vector<SpotValue> val;
  Date volDate;
  volDate.setScale("year");
  double volDateD=0;
  for(int i=0;i<m;i++){
    volDateD=(i+1)*.25;
    //for(int j=0; j<k;j++){
      val.push_back(SpotValue(volDate+volDateD, Vasicek_Volatility(a, sig,volDateD+delta, volDateD)));//for now, "constant volatility" and constant "delta"
    //}
  }

  std::vector<std::function<double(std::vector<double>&, std::vector<double>&)> > meanSquare;
  int n=val.size();
  std::vector<double> dataToMinimizeOver(n);
  std::vector<std::vector<double> > additionalParameters(n, std::vector<double>(2));
  for(int i=0; i<n; i++){
  //  auto mS=
    Date dt=val[i].date;
    //double valu=val[i].value;
    dataToMinimizeOver[i]=val[i].value;
    dt.setScale("year");
    Date currDate;
    double dtdiff=dt-currDate;
    additionalParameters[i][0]=dtdiff+delta;
    additionalParameters[i][1]=dtdiff;
    meanSquare.push_back([](std::vector<double> &guess, std::vector<double> &additionalParameters){
      return Vasicek_Volatility(guess[0], guess[1], additionalParameters[0], additionalParameters[1]);
    });
  }
  std::vector<double> guess(2);
  guess[0]=.3;
  guess[1]=.1;
  nt.optimize(meanSquare, additionalParameters, dataToMinimizeOver, guess);
  std::cout<<guess[0]<<", "<<guess[1]<<std::endl;
  /*End newton test*/

  std::cout<<Vasicek_Caplet(r0, a, b, sig, .04, delta, t)<<std::endl; //option on simple interest rate
  std::vector<SpotValue> yield;
  Date dt;
  dt.setScale("year");
  double delt=0;
  for(int i=0; i<20; i++){
    delt=(i+1)*.25;
    std::cout<<dt+delt<<": ";
    std::cout<<-log(Vasicek_Price(r0, a, b, sig, delt))/delt<<std::endl;
    yield.push_back(SpotValue(dt+delt, -log(Vasicek_Price(r0, a, b, sig, delt))/delt));
  }
  Vasicek vs(yield, a, sig, r0);
  std::cout<<vs.Bond_Price(t)<<std::endl;
  std::cout<<Vasicek_Price(r0, a, b, sig, t)<<std::endl;





  //BlackScholes bs(s0, bnd.getPrice(r0), k, t);
  //std::cout<<bs.getPut(bnd.getVasicekVolatility(delta))<<std::endl;
  //BlackScholes tst(40, exp(-.03), 40, 1);
  //std::cout<<tst.getCall(.3)<<std::endl;
  /*std::vector< std::vector<double> > sigs(3, std::vector<double>(3));
  sigs[0][0]=.3;
  sigs[1][0]=.2;
  sigs[2][0]=.4;
  sigs[0][1]=0;
  sigs[1][1]=.15;
  sigs[2][1]=.3;
  sigs[0][2]=0;
  sigs[1][2]=0;
  sigs[2][2]=.3;
  Matrix sig(sigs);
  int n=5;
  std::vector<ForwardValue> forwardLIBOR(n);
  //for(int i=0; i<n;i++){
  forwardLIBOR[0]={"1/1/2015", "3/31/2015", .03};
  forwardLIBOR[1]={"2/1/2015", "4/30/2015", .03};
  forwardLIBOR[2]={"3/1/2015", "5/31/2015", .03};
  forwardLIBOR[3]={"4/1/2015", "6/30/2015", .03};*/
  //}
  /*int ncol=5;
  int nrow=4;
  std::vector<std::vector<double> > matP(nrow, std::vector<double>(ncol, 0));
  std::vector<double> vP(ncol, 0);
  for(int i=0; i<ncol; i++){
    for(int j=0; j<nrow; j++){
      matP[j][i]=(i+1)*(j+1)*5;
    }
    vP[i]=i+3;
  }
  Matrix mat(matP);
  //std::cout<<vP[2]<<std::endl;
  Matrix mt=mat*vP;
  //std::cout<<mt.getNumRow()<<std::endl;
  for(int i=0; i<nrow; i++){
    std::cout<<mt[i][0]<<std::endl;
  }*/

/*  Double x=5;
  std::cout<<x*5<<std::endl;
  Date expiration("1/1/2016");
  Date current;
  std::cout<<expiration-current<<std::endl;
  EulerSimulation<Matrix> ms(100, 100000, 1);

  auto start = std::chrono::system_clock::now();
  std::vector<std::vector<double> > init(3, std::vector<double>(1));
  init[0][0]=.8;
  init[1][0]=1.1;
  init[2][0]=1.3;
  ms.simulate(Matrix(init), alpha, [&](double t, Matrix x){return sig;});
  auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);
	std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
  std::cout<<ms.getResult()<<std::endl;
  std::cout<<ms.getMCVol()<<std::endl;*/




  /*int numL=20;
	std::vector<double> libor(numL);
  for(int i=0; i<numL; i++){
    libor[i]=.01*i;
  }
  Swap swp(libor, .5);
  double rate=swp.getRate(5);
  std::cout<<rate<<std::endl;
  std::cout<<swp.getPrice(5, rate)<<std::endl; //this is zero (as it should be)
  std::cout<<swp.getBondYield(5)<<std::endl; //should be nearish the swap rate*/
}
