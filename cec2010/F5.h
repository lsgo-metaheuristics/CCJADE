#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F5 :public Benchmarks{
	protected:
	public:
		F5();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F5();
	};
}

