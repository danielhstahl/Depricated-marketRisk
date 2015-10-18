#ifndef __SIMULNORM_H_INCLUDED__
#define __SIMULNORM_H_INCLUDED__
#include <random>
#include <chrono>
#include <type_traits>
template<typename numberType>
class SimulNorm{
private:
  std::mt19937_64 generator;
  std::normal_distribution<double> norm;
  int seed;
  int n=0;
  numberType result;
public:
  SimulNorm(int n_){
    seed=std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    n=n_;
    norm=std::normal_distribution<double>(0.0, 1.0);
  }
  SimulNorm(int n_, int seed){
    generator.seed(seed);
    n=n_;
    norm=std::normal_distribution<double>(0.0, 1.0);
  }
  SimulNorm(int n_, double mu, double sigma, int seed){
    generator.seed(seed);
    norm=std::normal_distribution<double>(mu, sigma);
    n=n_;
  }
  SimulNorm(int n_, double mu, double sigma){
    seed=std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    n=n_;
    norm=std::normal_distribution<double>(mu, sigma);
  }
  /*template<typename num=numberType>
  typename std::enable_if<!std::is_fundamental<num>::value, void >::type setN(num init){
    n=init.size();
  }
  template<typename num=numberType>
  typename std::enable_if<std::is_fundamental<num>::value, void >::type setN(num init){
    n=1;
  }
  template<typename num=numberType>
  typename std::enable_if<std::is_fundamental<num>::value, double >::type getNorm(){
    return norm(generator);
  }
  template<typename num=numberType>
  typename std::enable_if<!std::is_fundamental<num>::value, std::vector<double> >::type*/
  void setNorm(){
    if(n==1){
      result=numberType(norm(generator)); //will this work???
    }
    else{
      std::vector<double> v(n);
      for(int i=0; i<n; i++){
        v[i]=norm(generator);
      }
      result=numberType(v);
    }

    //return v; //hopefully optimizer takes care of the copy
  }
  numberType getNorm(){
    setNorm();
    return result;
  }
};
#endif
