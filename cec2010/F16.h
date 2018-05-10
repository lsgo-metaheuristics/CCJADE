#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F16 :public Benchmarks{
	protected:
		void generateInterArray();
	public:
		F16();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F16();
	};
}


