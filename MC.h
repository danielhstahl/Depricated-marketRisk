#ifndef __MC_H_INCLUDED__
#define __MC_H_INCLUDED__
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
class MC {
	private:
    int m;
		double estimate;
		double error;
    std::vector<double> distribution;
	public:
		MC(int);
		double getEstimate();
		double getError();
		double getVaR(double);
    std::vector<double> getDistribution();
		template<typename FN>
		void simulate(FN&& fn) {
      estimate=0;
      error=0;
			#pragma omp parallel//multithread using openmp
			{
				#pragma omp for //multithread using openmp
				for(int j=0; j<m; j++){
					//int tid = omp_get_thread_num();
					//std::cout<<tid<<std::endl;
          estimate+=fn();
          error+=estimate*estimate;
				}
			}
      estimate=estimate/(double)m;
      error=(error/(double)m-estimate*estimate)/(double)m;
      error=sqrt(error);
		}
    template<typename FN>
		void simulateDistribution(FN&& fn) {
		  estimate=0;
		  error=0;
		  distribution=std::vector<double>(m);
		  #pragma omp parallel//multithread using openmp
		  {
		    #pragma omp for //multithread using openmp
		    for(int j=0; j<m; j++){
		      distribution[j]=fn();
					//std::cout<<distribution[j]<<std::endl;
		      estimate+=distribution[j];
		      error+=distribution[j]*distribution[j];
		    }
		  }
		  estimate=estimate/(double)m;
		  error=(error/(double)m-estimate*estimate)/(double)m;
		  error=sqrt(error);
		}
};
#endif
