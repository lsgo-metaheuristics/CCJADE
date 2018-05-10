#pragma once

#include "Benchmarks2010.h"
namespace CEC2010 {
	class F18 :public Benchmarks{
	protected:
		void generateInterArray();
	public:
		F18();
		tFitness compute(double* x);
		tFitness compute(vector<double> x);
		~F18();
	};
}


