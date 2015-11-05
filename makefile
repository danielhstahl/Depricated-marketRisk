LDFLAGS=-L../NewtonOptimization -lNewton -L../MiscellaniousUtilities -lDate  -L../eigen -L../BinomialTree -lTree
INCLUDES=-I../NewtonOptimization -I../MiscellaniousUtilities -I../eigen -I../BinomialTree

marketRisk: main.o Vasicek.o Swap.o BlackScholes.o MC.o Spline.o
	g++ -std=c++11 -O3  -w main.o BlackScholes.o Swap.o Vasicek.o MC.o Spline.o $(LDFLAGS) $(INCLUDES) -o marketRisk -fopenmp

main.o: main.cpp Vasicek.h MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -w -c main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Vasicek.o: Vasicek.cpp MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -w -c Vasicek.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Swap.o: Swap.cpp MarketData.h
	g++ -std=c++11 -O3  -w -c Swap.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

BlackScholes.o: BlackScholes.cpp
	g++ -std=c++11 -O3  -w -c BlackScholes.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

MC.o: MC.cpp
	g++ -std=c++11 -O3  -w -c MC.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Spline.o: Spline.cpp
		g++ -std=c++11 -O3  -w -c Spline.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o marketRisk
