//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : CCDE.cpp
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
//=============================================================================================

#include "CCDE.h"
#include "Decomposer.h"
#include "JADE.h"
#include <ctime>

using namespace std;


//******************************************************************************************/
//
//
//
//******************************************************************************************/
CCDE::CCDE()
{
    JADE_c = 0.05;
    JADE_p = 0.05;

    JADE_mutationStrategy = 1;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer * CCDE::createDecomposer(unsigned sizeOfSubcomponents, unsigned individualsPerSubcomponent, DG2 &dg2, bool RG)
{
    vector<unsigned> allCoordinates;
    for(unsigned i=0; i<problemDimension; ++i)
        allCoordinates.push_back(i);

    //if ( RG )
      //  shuffle(allCoordinates.begin(), allCoordinates.end(), eng);

    unsigned seed = unifRandom(eng)*100000;

    return new Decomposer(*this, seed, allCoordinates, sizeOfSubcomponents, individualsPerSubcomponent, population, contextVector, optimum, dg2, RG);
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::initPopulation(unsigned numOfIndividuals)
{
    population.clear();
    fitnessValues.clear();

    for (unsigned i = 0; i<numOfIndividuals; ++i)
    {
        vector< double > position;
        for (unsigned int d = 0; d<problemDimension; ++d)
            position.push_back(lowerLimit + unifRandom(eng)*(upperLimit - lowerLimit));
        population.push_back(position);
        fitnessValues.push_back(0);
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness CCDE::computeFitnessValue(vector<double> &x)
{
    tFitness f = fitness->compute(&x[0]);
    return f;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
CCDE::~CCDE()
{

}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::initContextVector()
{
    contextVector.resize(problemDimension);

    ///Initialize the context vector (just pick the first individual, which is randomly generated)
    for (unsigned i = 0; i < problemDimension; ++i)
        contextVector[i] = population[0][i];

    globalBestFitness = computeFitnessValue(contextVector);
    numberOfEvaluations++;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness CCDE::getFinalFitnessValue()
{
    return globalBestFitness;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::optimize(Fitness* _function, unsigned dim, double domain, tFitness _optimum,
	unsigned int _maxNumberOfEvaluations,
	unsigned _sizeOfSubcomponents,
	unsigned individualsPerSubcomponent,
	vector<ConvPlotPoint> &convergence,
	int seed,
	unsigned numItePerCycle, DG2 &dg2)

{
	eng.seed(seed);

	optimum = _optimum;

	numberOfEvaluations = dg2.getNumEvaluations();

	maxNumberOfEvaluations = _maxNumberOfEvaluations;

	fitness = _function;

	clock_t startTime = clock();

	problemDimension = dim;

	lowerLimit = -domain;

	upperLimit = domain;

	initPopulation(individualsPerSubcomponent);

	initContextVector();

	Decomposer *dec = createDecomposer(_sizeOfSubcomponents, individualsPerSubcomponent, dg2, true);	

	convergence.push_back(ConvPlotPoint(numberOfEvaluations, fabs(globalBestFitness - optimum), 0.0));

	cout << "Cycle=" << 0 << "  NOE=" << numberOfEvaluations << "  err=" << std::scientific << fabs(globalBestFitness - optimum) << endl;

	for (ite = 0; numberOfEvaluations<maxNumberOfEvaluations; ++ite)
	{
		clock_t begin = clock();
		for (unsigned j = 0; j < dec->optimizers.size(); j++)
		{
			JADE *optimizer = dec->optimizers[j];	    
			optimizer->loadIndividuals(dec->population);
		    optimizer->optimize(numItePerCycle);
			optimizer->storeIndividuals(dec->population);
			numberOfEvaluations += optimizer->nfe;
			optimizer->nfe = 0;
		}

		dec->buildContextVector();

		tFitness error = fabs(dec->bestAchievedFitness - optimum);

		clock_t end = clock();
		double timeSec = (end - begin) / static_cast<double>(CLOCKS_PER_SEC);

		cout << "Cycle=" << ite + 1 << "  NOE=" << numberOfEvaluations << "  err=" << std::scientific << error << endl;

		convergence.push_back(ConvPlotPoint(numberOfEvaluations, error, 0.0));
	}
	clock_t stopTime = clock();
	elapsedTime = ((double)(stopTime - startTime)) / CLOCKS_PER_SEC;
	cout << "elapsed time = " << elapsedTime << " s" << endl;
	delete dec;
}

//******************************************************************************************/
//
//
//
//******************************************************************************************/
void CCDE::optimizeSubcomponents(Decomposer *dec, unsigned nGenPerIteration)
{
    for (unsigned j = 0; j < dec->optimizers.size(); ++j)
    {
        JADE *optimizer = dec->optimizers[j];
        optimizer->loadIndividuals(dec->population);
        optimizer->optimize(nGenPerIteration);
        optimizer->storeIndividuals(dec->population);
    }
    dec->buildContextVector();
}

