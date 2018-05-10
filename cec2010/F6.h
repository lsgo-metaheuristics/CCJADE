#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F6 :public Benchmarks{
	protected:
	public:
		F6();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F6();
	};
}

