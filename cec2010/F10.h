#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F10 :public Benchmarks{
	protected:
	public:
		F10();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F10();
	};
}

