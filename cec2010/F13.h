#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F13 :public Benchmarks{
	protected:
	public:
		F13();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F13();
	};
}


