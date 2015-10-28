#ifndef __SWAP_H_INCLUDED__
#define __SWAP_H_INCLUDED__
#include <vector>
#include <cmath>
#include <iostream> //for debugging
#include "Date.h"
#include "MarketData.h"
#include "Tree.h"
typedef double Tenor;
typedef std::vector<double> Bonds1;
typedef std::vector<double> Bonds2;
typedef double Maturity;
typedef double SwapRate;
class Swap {
	private:
	//	std::vector<double> libor;
		std::vector<ForwardValue> libor;
    double tenor;
    int n;
    std::vector<ForwardValue> bonds;

	public:
		//Swap(const std::vector<double>& ,  double);
    //Swap(const std::vector<SpotValue>&, double);
    Swap(const std::vector<ForwardValue>&);
    Swap(double);
		double getPrice(double, double); //returns price for a given time to maturity and swap rate
		double getPrice(double, double,  double); //returns price at times when interval isn't exact
    double getRate(double);
    double getRate(std::vector<double>&, std::vector<double>&);
    double getPrice(std::vector<double>&, std::vector<double>&, double);
    double getBondPrice(double);
    double getBondYield(double);
};

double getSwapRate(Bonds1&, Bonds2&, Tenor);
double getSwapPrice(Bonds1&, Bonds2&, Tenor, Maturity);
template<typename ALPHA, typename SIGMA, typename FINV, typename PAYOFF>
double getSwaptionPrice(int m, double init, Maturity T, ALPHA&& alph, SIGMA&& sig, FINV&& inv, PAYOFF&& payoff){
	auto discount=[](double t, double x, double dt){
		return exp(-x*dt);
	};
	/*auto payoff=[](double x, k){
		if(x>k){
			return x-k;
		}
		else{
			return 0;
		}
	};*/
	Tree tree(m, T, init);
	tree.execute(alph, sig, inv, payoff, discount);
	return tree.getPrice();
}

#endif
