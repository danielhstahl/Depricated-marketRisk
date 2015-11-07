#ifndef __SPLINE_H_INCLUDED__
#define __SPLINE_H_INCLUDED__
#include <cmath>
#include <vector>

std::vector<double> spline(std::vector<double>&, std::vector<double>&);
double splint(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);
double splintD(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);


#endif
