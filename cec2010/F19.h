#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F19 :public Benchmarks{
	protected:
	public:
		F19();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F19();
	};
}



