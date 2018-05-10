//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Decomposer.cpp
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
//=============================================================================================


#include "Decomposer.h"
#include "CCDE.h"
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>

//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer::Decomposer(CCDE &_CCOptimizer,
                       unsigned seed,
                       vector<unsigned> &_coordinates,
                       unsigned _sizeOfSubcomponents,
                       unsigned _individualsPerSubcomponent,
                       vector< vector<double> > &_population,
                       vector<double>  &_contextVector,
                       tFitness _expectedOptimum,
	                   DG2 &_dg2,
                       bool RG,
                       bool allocateOptimizers) : CCOptimizer(_CCOptimizer), sizeOfSubcomponents(_sizeOfSubcomponents),
    individualsPerSubcomponent(_individualsPerSubcomponent), dg2(_dg2), applyRandomGrouping(RG)
{
    eng.seed(seed);
    expectedOptimum = _expectedOptimum;
    bestAchievedFitness = std::numeric_limits<tFitness>::infinity();
    coordinates = _coordinates;

    for (unsigned i = 0; i < individualsPerSubcomponent; ++i)
        population.push_back(_population[i]);

	//create the subcomponents including all separable variables
	int	q = 0;
	int remainingVariables = dg2.seps.size();
	while ( remainingVariables>0 )
	{
		int size = min(sizeOfSubcomponents, remainingVariables);
		remainingVariables -= size;
		if (remainingVariables < sizeOfSubcomponents/2)
		{
			size += remainingVariables;
			remainingVariables = 0;
		}

		cout << "creating a separable subcomponent of size " << size << endl;
		
		unsigned *optCoord = &(dg2.seps[q]);
		
		sizes.push_back(size);
		for (int j = 0; j < individualsPerSubcomponent; j++)
		{
			for (int k = 0; k < size; ++k)
			{
				population[j][optCoord[k]] = CCOptimizer.lowerLimit + unifRandom(eng) * (CCOptimizer.upperLimit - CCOptimizer.lowerLimit);
			}
		}
		JADE *optimizer = new JADE(size, individualsPerSubcomponent, *this);
		optimizer->setCoordinates(optCoord, size);
		optimizers.push_back(optimizer);
		q += size;
		
	}
	cout << "remainingVariables =" << remainingVariables << endl;
	//create the nonseparable subcomponents
	if( dg2.nonseps.size()>0 )
	  cout << "creating " << dg2.nonseps.size() << " nonseparable subcomponenta" << endl;
	for (int i = 0; i < dg2.nonseps.size(); ++i)
	{
		cout << "   creating a nonseparable subcomponent of size " << dg2.nonseps[i].size() << endl;
		sizes.push_back(dg2.nonseps[i].size());
		unsigned *optCoord = &(dg2.nonseps[i][0]);
		for (int j = 0; j < individualsPerSubcomponent; j++)
		{
			for (int k = 0; k < dg2.nonseps[i].size(); ++k)
			{
				population[j][optCoord[k]] = CCOptimizer.lowerLimit + unifRandom(eng) * (CCOptimizer.upperLimit - CCOptimizer.lowerLimit);
			}
		}
		JADE *optimizer = new JADE(dg2.nonseps[i].size(), individualsPerSubcomponent, *this);
		optimizer->setCoordinates(optCoord, dg2.nonseps[i].size());
		optimizers.push_back(optimizer);
	}

    numberOfSubcomponents = sizes.size();

    optimizers.resize(numberOfSubcomponents);

    contextVector.resize(CCOptimizer.problemDimension);
    _contextVector.resize(CCOptimizer.problemDimension);
    for (unsigned i = 0; i < CCOptimizer.problemDimension; ++i)
        contextVector[i] = _contextVector[i] = population[unifRandom(eng)*population.size()][i];


	for (unsigned i = 0; i < optimizers.size(); ++i)
	{
		optimizers[i]->loadIndividuals(population);
		optimizers[i]->evaluateParents();
		CCOptimizer.numberOfEvaluations += optimizers[i]->nfe;
		optimizers[i]->nfe = 0;
	}

}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
vector< JADE* >  Decomposer::allocateOptimizers(vector<unsigned> &indexes)
{
	vector< JADE* > bunchOfOptimizers;

	for (int i = 0; i < indexes.size(); ++i)
	{
		unsigned j = indexes[i];
		JADE *optimizer = new JADE(sizes[j], individualsPerSubcomponent, *this);
		optimizer->setCoordinates(&coordinates[baseCoordIndex[j]], sizes[j]);
		optimizers[j] = optimizer;
		bunchOfOptimizers.push_back(optimizer);
		optimizer->loadIndividuals(population);
	}
	return bunchOfOptimizers;
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
vector< JADE* >  Decomposer::allocateOptimizers(list< vector<unsigned> > &decomposition)
{
	vector< JADE* > bunchOfOptimizers;

	for (int i = 0; i < optimizers.size(); ++i)
		delete optimizers[i];

	numberOfSubcomponents = decomposition.size();
	optimizers.resize(numberOfSubcomponents);

	int i = 0;
	for (auto ite = decomposition.begin(); ite!=decomposition.end(); ++ite)
	{
		JADE *optimizer = new JADE((*ite).size(), individualsPerSubcomponent, *this);
		optimizer->setCoordinates(&(*ite)[0], (*ite).size());
		optimizers[i++] = optimizer;
		bunchOfOptimizers.push_back(optimizer);
	}
	return bunchOfOptimizers;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
vector< JADE* >  Decomposer::allocateOptimizers(vector< vector<unsigned> > &decomposition)
{
	vector< JADE* > bunchOfOptimizers;

	for (int i = 0; i < optimizers.size(); ++i)
		delete optimizers[i];

	numberOfSubcomponents = decomposition.size();
	optimizers.resize(numberOfSubcomponents);

	int i = 0;
	for (auto ite = decomposition.begin(); ite != decomposition.end(); ++ite)
	{
		JADE *optimizer = new JADE((*ite).size(), individualsPerSubcomponent, *this);
		optimizer->setCoordinates(&(*ite)[0], (*ite).size());
		optimizers[i++] = optimizer;
		bunchOfOptimizers.push_back(optimizer);
	}
	return bunchOfOptimizers;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
JADE* Decomposer::allocateOptimizer()
{
    JADE *optimizer = new JADE(sizes[0], individualsPerSubcomponent, *this);
    optimizer->setCoordinates(&coordinates[baseCoordIndex[0]], sizes[0]);
    optimizer->loadIndividuals(population);
    return optimizer;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void  Decomposer::setOptimizersCoordinatesAndEvaluatePopulation(vector<unsigned> &indexes)
{
    for (int i = 0; i < indexes.size(); ++i)
    {
        unsigned j = indexes[i];
        JADE *optimizer = optimizers[j];
        optimizer->setCoordinates(&coordinates[baseCoordIndex[j]], sizes[j]);
        optimizer->updateIndividuals(population);
        optimizer->evaluateParents();
        optimizer->updateIndexOfBest();
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void  Decomposer::setOptimizersCoordinates(vector<unsigned> &indexes)
{
    for (int i = 0; i < indexes.size(); ++i)
    {
        unsigned j = indexes[i];
        JADE *optimizer = optimizers[j];
        optimizer->setCoordinates(&coordinates[baseCoordIndex[j]], sizes[j]);
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setPopulation(vector< vector<double> > &_population)
{
    population.clear();
    for (unsigned i = 0; i < individualsPerSubcomponent; ++i)
        population.push_back(_population[i]);
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setOptimizersCoordinatesAndEvaluatePopulation()
{
    unsigned d = 0, k = 0, size = sizeOfSubcomponents;
    while ( d<coordinates.size() )
    {
        if (d + size > coordinates.size())
            size = coordinates.size() - d;

        JADE *optimizer = optimizers[k];

        optimizer->setCoordinates(&(coordinates[d]), size);

        optimizer->loadIndividuals(population);
        optimizer->evaluateParents();
        optimizer->updateIndexOfBest();

        d += size;
        k++;
    }

}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
Decomposer::~Decomposer()
{
    for (unsigned i = 0; i < optimizers.size(); ++i)
        delete optimizers[i];
};



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setSeed(unsigned seed)
{
    eng.seed(seed);
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setSubcomponentsOfEqualSize(unsigned newSizeOfSubcomponents, int individualsPerSubcomponent)
{
    for (unsigned i = 0; i < optimizers.size(); ++i)
        delete optimizers[i];
    optimizers.clear();

    sizeOfSubcomponents = newSizeOfSubcomponents;
    unsigned numberOfSubcomponents = coordinates.size() / newSizeOfSubcomponents;


    for (unsigned i = 0; i<numberOfSubcomponents; ++i)
    {
        optimizers.push_back(new JADE(sizeOfSubcomponents, individualsPerSubcomponent, *this));

        optimizers[i]->setCoordinates(&(coordinates[i*sizeOfSubcomponents]), sizeOfSubcomponents);

        optimizers[i]->loadIndividuals(population);
		optimizers[i]->evaluateParents();
        optimizers[i]->updateIndexOfBest();
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::setCoordinates(vector<unsigned> &_coordinates)
{
    coordinates = coordinates;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::updateContextVector(JADE *optimizer)
{
    vector<double> v = optimizer->getCollaborator();
    tFitness newBestCandidate = optimizer->calculateFitnessValue(v);
    CCOptimizer.numberOfEvaluations++;
   // if ( newBestCandidate < bestAchievedFitness )
    {
        for (unsigned ld = 0; ld<v.size(); ld++)
            contextVector[optimizer->coordinates[ld]] = v[ld];

        bestAchievedFitness = newBestCandidate;
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::buildContextVector()
{
	
	for (unsigned j = 0; j<optimizers.size(); ++j)
    {	
		optimizers[j]->nfe = 0;
        vector<double> v = optimizers[j]->getCollaborator();
        //tFitness newBestCandidate = optimizers[j]->calculateFitnessValue(v);
        //if ( newBestCandidate<bestAchievedFitness )
        {
            for (unsigned ld = 0; ld<v.size(); ld++)
                contextVector[optimizers[j]->coordinates[ld]] = v[ld];
            //bestAchievedFitness = newBestCandidate;
        }
    }
	bestAchievedFitness = this->CCOptimizer.computeFitnessValue(contextVector);
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void Decomposer::randomGrouping()
{
    if ( optimizers.size() && this->applyRandomGrouping )
    {
        shuffle(coordinates.begin(), coordinates.end(), eng);
        unsigned numOfCoordinatesPerSubgroup = coordinates.size() / optimizers.size();
        for (unsigned i = 0; i < optimizers.size(); ++i)
        {
            optimizers[i]->setCoordinates(&(coordinates[i*numOfCoordinatesPerSubgroup]), numOfCoordinatesPerSubgroup);
            optimizers[i]->loadIndividuals(population);
            optimizers[i]->evaluateParents();
            CCOptimizer.numberOfEvaluations += optimizers[i]->nfe;
            optimizers[i]->nfe = 0;
        }
    }
}


