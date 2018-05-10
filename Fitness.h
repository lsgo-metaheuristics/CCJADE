#ifndef _FITNESS_H
#define _FITNESS_H

#include <sstream>
#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdint>
using namespace std;


#define _PI (3.141592653589793238462643383279)
#define _E  (2.718281828459045235360287471352)
#define _L(i) ((int64_t)i)
#define _D(i) ((double)i)

#define tFitness long double

class Fitness
{
protected:	
	double minX;
	double maxX;
	unsigned dimension;
	unsigned ID;
public:
	vector< vector<unsigned> > idealGroups;

	Fitness() {};
	virtual ~Fitness() {  };
	virtual tFitness compute(double* x) = 0;
	virtual tFitness compute(vector<double> x) = 0;
	virtual void createIdealGroups() = 0;
	virtual double getMinX() = 0;
	virtual double getMaxX() = 0;
	virtual unsigned getID() = 0;
	virtual unsigned  getDimension() = 0;
	void setDimension(unsigned d) {};
};

#endif