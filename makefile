LDFLAGS=-L ../NewtonOptimization -lNewton -L ../MiscellaniousUtilities -lDate  -L../eigen
INCLUDES=-I ../NewtonOptimization -I ../MiscellaniousUtilities -I ../eigen
marketRisk: main.o Vasicek.o BlackScholes.o
	g++ -std=c++11 -O3  main.o BlackScholes.o Vasicek.o  $(LDFLAGS) $(INCLUDES) -o marketRisk -fopenmp

main.o: main.cpp Vasicek.h MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
	
Vasicek.o: Vasicek.cpp MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -c Vasicek.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

BlackScholes.o: BlackScholes.cpp 
	g++ -std=c++11 -O3  -c BlackScholes.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o marketRisk
