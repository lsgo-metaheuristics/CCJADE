#pragma once

#include "Benchmarks2010.h"

namespace CEC2010 {
	class F2 :public Benchmarks{
	protected:

	public:
		F2();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F2();
	};
}

