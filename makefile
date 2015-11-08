LDFLAGS=-L../NewtonOptimization -lNewton -L../MiscellaniousUtilities -lDate  -L../eigen -L../BinomialTree -lTree
INCLUDES=-I../NewtonOptimization -I../MiscellaniousUtilities -I../eigen -I../BinomialTree

marketRisk: main.o Vasicek.o Swap.o BlackScholes.o MC.o Spline.o YieldSpline.o YieldPolynomial.o YieldNelsonSiegal.o
	g++ -std=c++11 -O3  -w -fPIC main.o BlackScholes.o Swap.o Vasicek.o MC.o Spline.o YieldSpline.o YieldPolynomial.o YieldNelsonSiegal.o $(LDFLAGS) $(INCLUDES) -o marketRisk -fopenmp

main.o: main.cpp Vasicek.h MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -w -c -fPIC main.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Vasicek.o: Vasicek.cpp MarketData.h BlackScholes.h
	g++ -std=c++11 -O3  -w -c -fPIC Vasicek.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Swap.o: Swap.cpp MarketData.h
	g++ -std=c++11 -O3  -w -c -fPIC Swap.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

BlackScholes.o: BlackScholes.cpp
	g++ -std=c++11 -O3  -w -c -fPIC BlackScholes.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

MC.o: MC.cpp
	g++ -std=c++11 -O3  -w -c -fPIC MC.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

Spline.o: Spline.cpp
		g++ -std=c++11 -O3  -w -c -fPIC Spline.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

YieldSpline.o: YieldSpline.cpp
		g++ -std=c++11 -O3  -w -c -fPIC YieldSpline.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

YieldPolynomial.o: YieldPolynomial.cpp
		g++ -std=c++11 -O3  -w -c -fPIC YieldPolynomial.cpp $(LDFLAGS) $(INCLUDES) -fopenmp

YieldNelsonSiegal.o: YieldNelsonSiegal.cpp
		g++ -std=c++11 -O3  -w -c -fPIC YieldNelsonSiegal.cpp $(LDFLAGS) $(INCLUDES) -fopenmp
clean:
	     -rm *.o marketRisk
