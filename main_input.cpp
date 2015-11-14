#include <iostream>
#include <string>
#include "document.h" //rapidjson
#include "writer.h" //rapidjson
#include "stringbuffer.h" //rapidjson
#include "Vasicek.h"
#include "YieldSpline.h"
#include "MarketData.h"
#include "Date.h"
#include "MC.h"
#include "SimulNorm.h"
#include "ComputePortfolio.h"

typedef Vasicek<YieldSpline> VasicekEngine;
//std::string jsonYield(*v8::String::Utf8Value(args[0]->ToString()));
//std::string jsonHistorical(*v8::String::Utf8Value(args[1]->ToString()));
//std::string jsonPortfolio(*v8::String::Utf8Value(args[2]->ToString()));
int main(){
  std::string yieldData;
  std::string historicalData;
  std::string portfolioData;
  for (yieldData; std::getline(std::cin, yieldData);) {
    break;
  }
  for (historicalData; std::getline(std::cin, historicalData);) {
    break;
  }
  for (portfolioData; std::getline(std::cin, portfolioData);) {
    break;
  }
  rapidjson::Document dYield;
  rapidjson::Document dHistorical;
  rapidjson::Document dPortfolio;
  dYield.Parse(yieldData.c_str());//yield data
  dHistorical.Parse(historicalData.c_str()); //historical data
  dPortfolio.Parse(portfolioData.c_str());//portfolio
  yieldData.clear();
  historicalData.clear();
  portfolioData.clear();
  int n=dYield.Size();
  int m=dHistorical["observations"].Size();
  int sizePortfolio=dPortfolio.Size();
  //Deal with yield
  std::vector<SpotValue> yield;
  Date dt;//current date
  dt.setScale("day");
  for(int i=0; i<n; ++i){
    yield.push_back(SpotValue(dt+dYield[i]["daysPlus"].GetInt(), std::stod(dYield[i]["value"].GetString())*.01));
  }
  YieldSpline spl(yield);
  //std::cout<<"test"<<spl.Yield(20)<<std::endl;
  //end deal with yield
  //Deal with historical
  std::vector<SpotValue> historical;
  std::string val;
  for(int i=0; i<m; ++i){
    val=dHistorical["observations"][i]["value"].GetString();
    if(val!="."){//nulls show up as "."
      historical.push_back(SpotValue(dHistorical["observations"][i]["date"].GetString(), std::stod(dHistorical["observations"][i]["value"].GetString())*.01));
    }
  }
  //end deal with historical
  //Deal with portfolio
  //std::cout<<"This is first string date: "<<dPortfolio[0]["maturity"].GetString()<<std::endl;
  //Date date1=Date("2017-11-05");
  //Date date2=Date("2014-10-02");
  //std::cout<<date1<<std::endl;
  //std::cout<<date2<<std::endl;
  std::vector<AssetFeatures> portfolio;
  for(int i=0; i<sizePortfolio; ++i){
    /*Date Maturity;
    Date UnderlyingMaturity;
    double Strike;
    double Tenor;
    std::string type;*/
    if(strcmp(dPortfolio[i]["type"].GetString(), "bond")==0){
      //newAsset.Maturity=Date(dPortfolio[i]["maturity"].GetString());
      //newAsset.type="bond";
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["maturity"].GetString()),
        0,
        0,
        "bond"
      });
      //portfolio.back().Maturity=Date(dPortfolio[i]["maturity"].GetString());
      //portfolio.back().type="bond";
    }
    else if(strcmp(dPortfolio[i]["type"].GetString(), "caplet")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["maturity"].GetString()),
        dPortfolio[i]["strike"].GetDouble(),
        dPortfolio[i]["tenor"].GetDouble(),
        "caplet"
      });
      /*portfolio.push_back(AssetFeatures());
      portfolio.back().Maturity=Date(dPortfolio[i]["maturity"].GetString());
      portfolio.back().type="caplet";
      portfolio.back().Strike=dPortfolio[i]["strike"].GetDouble();
      portfolio.back().Tenor=dPortfolio[i]["tenor"].GetDouble();*/
      //portfolio.push_back(AssetFeatures(dPortfolio[i]["maturity"].GetString(), dPortfolio[i]["strike"].GetDouble(), dPortfolio[i]["tenor"].GetDouble(), "caplet"));
    }
    else if(strcmp(dPortfolio[i]["type"].GetString(), "swaption")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["underlyingMaturity"].GetString()),
        dPortfolio[i]["strike"].GetDouble(),
        dPortfolio[i]["tenor"].GetDouble(),
        "swaption"
      });
      //portfolio.push_back(AssetFeatures(dPortfolio[i]["maturity"].GetString(), dPortfolio[i]["strike"].GetDouble(), dPortfolio[i]["tenor"].GetDouble(),dPortfolio[i]["underlyingMaturity"].GetString(), "swaption"));
    }
    else if(strcmp(dPortfolio[i]["type"].GetString(), "call")==0){
      portfolio.push_back({
        Date(dPortfolio[i]["maturity"].GetString()),
        Date(dPortfolio[i]["underlyingMaturity"].GetString()),
        dPortfolio[i]["strike"].GetDouble(),
        0,
        "call"
      });
    //  portfolio.push_back(AssetFeatures(dPortfolio[i]["maturity"].GetString(), dPortfolio[i]["strike"].GetDouble(), dPortfolio[i]["underlyingMaturity"].GetString(), "call"));
    }

    //std::cout<<dPortfolio[i]["maturity"].GetString()<<std::endl;
    //std::cout<<portfolio[0].Maturity<<std::endl;
  }
  //end deal with portfolio


  double a=.4;
  double sig=.02;
  double r0=.03;
  VasicekEngine vsk(spl, a, sig, r0);
  vsk.findHistoricalB(historical);

  //std::cout<<"swaption Test: "<<vsk.Swaption(.01, 0, .03, .25, 2, .02777)<<std::endl;

  SimulNorm simul;
  ComputePortfolio<VasicekEngine> myPortfolio(&vsk, portfolio);
  myPortfolio.setDates(10); //ten day VaR

  //std::cout<<"bond call: "<<vsk.Bond_Call(.0238855, .02777778, .99, .35, .1)<<std::endl;

  MC mc(1000);//1000 simulations 

  auto start = std::chrono::system_clock::now();
  mc.simulateDistribution([&](){
    return myPortfolio.execute(vsk.simulate(simul));//simulates path and computes portfolio value along the path
  });
  auto end=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start);

  std::cout<<"Time it took: "<<end.count()/1000.0<<std::endl;
  std::cout<<"99VaR: "<<mc.getVaR(.99)<<std::endl;
  std::cout<<"Error: "<<mc.getError()<<std::endl;

  //Local<Object> obj = Object::New(isolate);
  //obj->Set(String::NewFromUtf8(isolate, "VaR"), Number::New(isolate, mc.getVaR(.99)));
  //obj->Set(String::NewFromUtf8(isolate, "Error"), Number::New(isolate, mc.getError()));
}
