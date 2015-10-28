LDFLAGS=-L ../NewtonOptimization -lNewton -L ../MiscellaniousUtilities -lDate  -L../eigen -L ../BinomialTree -lTree
INCLUDES=-I ../NewtonOptimization -I ../MiscellaniousUtilities -I ../eigen -I ../BinomialTree

marketRisk: main.o Vasicek.o Swap.o BlackScholes.o MC.o 
	g++ -std=c++11 -O3  main.o BlackScholes.o Swap.o Vasicek.o MC.o  $(LDFLAGS) $(INCLUDES) -o marketRisk -fopenmp

main.o: main.cpp Vasicek.h MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Vasicek.o: Vasicek.cpp MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -c Vasicek.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Swap.o: Swap.cpp MarketData.h
	g++ -std=c++11 -O3  -c Swap.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

BlackScholes.o: BlackScholes.cpp
	g++ -std=c++11 -O3  -c BlackScholes.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

MC.o: MC.cpp
	g++ -std=c++11 -O3  -c MC.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

clean:
	     -rm *.o marketRisk
