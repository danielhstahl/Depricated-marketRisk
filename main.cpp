#define _USE_MATH_DEFINES
#include <iostream>
#include "Swap.h"
#include <ctime>
#include <vector>
#include <chrono> //for accurate multithreading time using std::chrono


int main(){
  int numL=20;
	std::vector<double> libor(numL);
  for(int i=0; i<numL; i++){
    libor[i]=.01*i;
  }
  Swap swp(libor, .5);
  std::cout<<swp.getPrice(5, .03)<<std::endl;


}
