#include "BondUtilities.h"
double convertLiborToContinuous(double yield, double t){//rate, time
  return -log(1/(1+yield*t))/t;
}
double convertLiborToBond(double libor, double t){//rate, time
  return 1.0/(libor*t+1.0);
}
double convertBondToLibor(double bond, double t){//bond, time
    return (1.0/bond-1.0)/t;
}
