#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F12 :public Benchmarks{
	protected:
	public:
		F12();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F12();
	};
}

