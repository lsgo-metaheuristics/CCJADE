#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F11 :public Benchmarks{
	protected:
	public:
		F11();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F11();
	};
}

