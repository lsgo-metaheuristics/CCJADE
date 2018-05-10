#pragma once
#include "Benchmarks2010.h"
namespace CEC2010 {
	class F20 :public Benchmarks{
	protected:
	public:
		F20();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F20();
	};
}


