#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F8 :public Benchmarks{
	protected:
	public:
		F8();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F8();
	};
}

