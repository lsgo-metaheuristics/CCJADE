#pragma once

#include "Benchmarks2010.h"

namespace CEC2010 {
	class F3 :public Benchmarks{
	protected:

	public:
		F3();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F3();
	};
}

