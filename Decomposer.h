//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Decomposer.h
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
//=============================================================================================

#pragma once
#include "JADE.h"
#include <vector>
#include <queue>
#include <map>
#include <random>
#include <functional>
#include <numeric>
#include "DG2.hpp"

class CCDE;

using namespace std;

typedef mt19937 RandomEngine;


typedef struct
{
	int variable;
	double intensity;
} Interaction;

typedef struct
{
	int variable;
	list<Interaction> interactions;
} Interactions;



class Decomposer
{
public:
	struct doCompareOptimizers
	{
		doCompareOptimizers(const vector<JADE *> &opt) : optimizers(opt) { }
		const vector<JADE*> &optimizers;

		bool operator()(const int & i1, const int & i2)
		{
			return optimizers[i1]->bestFitness > optimizers[i2]->bestFitness;
		}
	};

    vector<unsigned> coordinates;
    CCDE &CCOptimizer;
    vector< JADE* > optimizers;
    int sizeOfSubcomponents;
    bool applyRandomGrouping;
    tFitness bestAchievedFitness;
    RandomEngine eng;
    uniform_real_distribution<double> unifRandom;
    unsigned individualsPerSubcomponent;
    unsigned numberOfSubcomponents;
	DG2 &dg2;
	int baseVariable;
	list<Interactions> allInteractionsB;

    //Current population
    vector< vector<double> > population;

    //Final global best position and context vector
    vector<double> contextVector;

    tFitness expectedOptimum;

    //Fitnesses of population
    vector< tFitness > fitnessValues;    

    vector<unsigned> sizes;
    vector<unsigned> baseCoordIndex;
	vector< list< vector<unsigned>>> allDecompositions;
	vector<double> interationLevels;
	vector<int> minSizes, maxSizes;

    Decomposer(CCDE &_CCOptimizer, unsigned seed, vector<unsigned> &_coordinates,
               unsigned _sizeOfSubcomponents,
               unsigned _individualsPerSubcomponent,
               vector< vector<double> > &_population,
               vector<double>  &_contextVector,
               tFitness _expectedOptimum,
		       DG2 &_dg2,
			   bool RG,
			   bool allocateOptimizers=true);
    ~Decomposer();
    vector< JADE* >  allocateOptimizers(vector<unsigned> &indexes);
	vector< JADE* >  allocateOptimizers(list< vector<unsigned> > &decomposition);
	vector< JADE* >  allocateOptimizers(vector< vector<unsigned> > &decomposition);
    JADE*  allocateOptimizer();
    void setPopulation(vector< vector<double> > &_population);
    void setSubcomponentsOfEqualSize(unsigned newSizeOfSubcomponents, int individualsPerSubcomponent);
    void setCoordinates(vector<unsigned> &_coordinates);
    void updateContextVector(JADE *optimizer);
    void buildContextVector();
    void randomGrouping();
    void setSeed(unsigned seed);
    void setOptimizersCoordinatesAndEvaluatePopulation();
    void setOptimizersCoordinatesAndEvaluatePopulation(vector<unsigned> &indexes);
    void setOptimizersCoordinates(vector<unsigned> &indexes);
};
