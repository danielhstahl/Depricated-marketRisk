#include "MC.h"
MC::MC(int m_) {
  m=m_;
}
double MC::getEstimate() {
	return estimate;
}
double MC::getError(){
	return error;
}
std::vector<double> MC::getDistribution(){
	return distribution;
}
double MC::getVaR(double q){ //eg, .99
  std::sort(distribution.begin(), distribution.end());
  return distribution[(int)((1.0-q)*m)];
  //return distribution[(int)((1.0-q)*m)];
}
