#include <iostream>
#include <string>
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson
#include "BondUtilities.h"
#include "Vasicek.h"
#include "YieldSpline.h"
#include "MarketData.h"
#include "Date.h"
#include "MC.h"
#include "SimulNorm.h"
#include "ComputePortfolio.h"
//#include <ios>

typedef Vasicek<YieldSpline> VasicekEngine;
int main(){
  std::string yieldData;
  std::string historicalData;
  std::string portfolioData;

  //strategy for yield curve: use future rates (goes out ~10 years) and adjust for forward rate.  then fit spline to this.  Then find implied zero bond prices.  Then get implied zero yield curve.


  for (yieldData; std::getline(std::cin, yieldData);) {
    break;
  }
  for (historicalData; std::getline(std::cin, historicalData);) {
    break;
  }

  rapidjson::Document dYield;
  rapidjson::Document dHistorical;

  dYield.Parse(yieldData.c_str());//yield data
  dHistorical.Parse(historicalData.c_str()); //historical data

  yieldData.clear();
  historicalData.clear();

  int n=dYield.Size();
  int m=dHistorical["observations"].Size();

  //Deal with yield
  std::vector<SpotValue> SwapRates;
  std::vector<SpotValue> LiborRates;
  Date dt;//current date
  dt.setScale("day");
  for(int i=0; i<n; ++i){
    if(strcmp(dYield[i]["type"].GetString(), "Swap")==0){
      SwapRates.push_back(SpotValue(dt+dYield[i]["daysPlus"].GetInt(), std::stod(dYield[i]["value"].GetString())*.01));
    }
    else{
      LiborRates.push_back(SpotValue(dt+dYield[i]["daysPlus"].GetInt(), std::stod(dYield[i]["value"].GetString())*.01));
    }

  }

  //end deal with yield
  //Deal with historical
  std::vector<SpotValue> historical;
  std::string val;
  for(int i=0; i<m; ++i){
    val=dHistorical["observations"][i]["value"].GetString();
    if(val!="."){//nulls show up as "."
      historical.push_back(SpotValue(dHistorical["observations"][i]["date"].GetString(), std::stod(dHistorical["observations"][i]["value"].GetString())*.01));
      //std::cout<<"val: "<<std::stod(dHistorical["observations"][i]["value"].GetString())*.01<<std::endl;
    }

  }
  //end deal with historical

  double a=.4; //made up for demonstrate purposes
  double sig=.02; //made up for demonstration purposes

  //r0=-log(1-r0*shortRatePeriod)/shortRatePeriod;//convert simple libor rate to continuously compounded rate
  //vs->setR(r0);
  double r0=convertLiborToContinuous(historical.back().value, 7.0/360.0);//convert simple libor rate to continuously compounded rate...make the 7 day convention known from the incoming object
  //-log(1/(1+historical.back().value*7.0/360.0))*(360.0/7.0);!
  //std::cout<<"n: "<<n<<std::endl;
  VasicekEngine vsk(a, sig, r0);
  YieldSpline spl;
  //spl.computeSimpleFutureSpline(r0, FuturePrices, vsk);
  spl.computeSimpleSwapSpline(LiborRates, SwapRates);
  vsk.setYield(&spl);


//Deal with portfolio
  for (portfolioData; std::getline(std::cin, portfolioData);) {
    break;
  }

  rapidjson::Document dPortfolio;
  std::vector<AssetFeatures> portfolio;
  dPortfolio.Parse(portfolioData.c_str());//portfolio
  portfolioData.clear();
  int sizePortfolio=dPortfolio.Size();

  for(int i=0; i<sizePortfolio; ++i){
    /*Date Maturity;
    Date UnderlyingMaturity;
    double Strike;
    double Tenor;
    std::string type;*/
    if(strcmp(dPortfolio[i]["type"].GetString(), "bond")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["maturity"].GetString()),
        0,
        0,
        "bond"
      });
    }
    else if(strcmp(dPortfolio[i]["type"].GetString(), "caplet")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["maturity"].GetString()),
        dPortfolio[i]["strike"].GetDouble(),
        dPortfolio[i]["tenor"].GetDouble(),
        "caplet"
      });
    }
    else if(strcmp(dPortfolio[i]["type"].GetString(), "swaption")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["underlyingMaturity"].GetString()),
        dPortfolio[i]["strike"].GetDouble(),
        dPortfolio[i]["tenor"].GetDouble(),
        "swaption"
      });
    }
    else if(strcmp(dPortfolio[i]["type"].GetString(), "call")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["underlyingMaturity"].GetString()),
        dPortfolio[i]["strike"].GetDouble(),
        0,
        "call"
      });
    }
  }
  //end deal with portfolio



  vsk.findHistoricalB(historical);

  SimulNorm simul;
  ComputePortfolio<VasicekEngine> myPortfolio(&vsk, portfolio);
  myPortfolio.setDates(10); //ten day VaR

  MC mc(1000);//1000 simulations

  auto start = std::chrono::system_clock::now();
  mc.simulateDistribution([&](){
    return myPortfolio.execute(vsk.simulate(simul));//simulates path and computes portfolio value along the path
  });
  auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

  std::cout<<"{\"time\": "<<end.count()/1000.0<<", \"VaR\": "<<mc.getVaR(.99)<<", \"MCError\": "<<mc.getError()<<"}"<<std::endl;
}
