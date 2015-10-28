#ifndef __SIMULNORM_H_INCLUDED__
#define __SIMULNORM_H_INCLUDED__
#include <random>
#include <chrono>
class SimulNorm{
private:
  std::mt19937_64 generator;
  std::normal_distribution<double> norm;
  int seed;
public:
  SimulNorm(){
    seed=std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    //n=n_;
    norm=std::normal_distribution<double>(0.0, 1.0);
  }
  SimulNorm(int seed){
    generator.seed(seed);
    //n=n_;
    norm=std::normal_distribution<double>(0.0, 1.0);
  }
  SimulNorm(double mu, double sigma, int seed){
    generator.seed(seed);
    norm=std::normal_distribution<double>(mu, sigma);
  }
  SimulNorm(double mu, double sigma){
    seed=std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
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
  double getNorm(){
    return norm(generator); //will this work???
    /*if(n==1){
      result=numberType(norm(generator)); //will this work???
    }
    else{
      std::vector<double> v(n);
      for(int i=0; i<n; i++){
        v[i]=norm(generator);
      }
      result=numberType(v);
    }*/

    //return v; //hopefully optimizer takes care of the copy
  }
  /*numberType getNorm(){
    setNorm();
    return result;
  }*/
};
#endif
