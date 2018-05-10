//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Jade.h
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
//=============================================================================================

#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <vector>
#include <deque>
#include <map>
#include <cstddef>
#include "Fitness.h"

class Decomposer;

//using namespace dlib;
using namespace std;

class JADE
{
	Decomposer &decomposer;

	struct doCompareIndividuals
	{
		doCompareIndividuals(const tFitness *_f) : f(_f) { }
		const tFitness *f;

		bool operator()(const int & i1, const int & i2)
		{
			return f[i1] < f[i2];
		}
	};



public:
	JADE(unsigned _dimension, unsigned _numberOfIndividuals, Decomposer &_group);
	~JADE();
	void setCoordinates(vector<unsigned> &_coordinates);
	void setCoordinates(unsigned *coordinates, unsigned numOfCoordinates);
	void update();
	void updateContextVector();
	void updateContextVector(vector<double> &cv, vector<unsigned> &coords, unsigned &vi);
	void sortPopulation(vector<tFitness> &fitness, vector<int> &sortIndex);
	void evaluateOffsprings(vector<bool> &toEvaluate, vector<bool> &hasTrueFitness);	
	int evaluateParents();	
	tFitness calculateFitnessValue(vector<double> &p, bool updateArchive=true);
	void optimize(int iterations);
	tFitness updateIndexOfBest();
	void loadIndividuals(vector< vector<double> > &population);
	void updateIndividuals(vector< vector<double> > &population);
	void storeIndividuals(vector< vector<double> > &population);
	void setParentFitness(vector<tFitness> &fitnessValues);		
	vector<double> &getCollaborator();		
	unsigned nfe;
	vector<unsigned> coordinates;	
	unsigned int dimension;	
	unsigned int ite;
	double	JADE_mu_cr;
	double	JADE_mu_ff;
	tFitness lastFitnessGain;

	///array containing the positions of all individuals
	vector< vector<double> > parents;
	vector< vector<double> > offsprings;
	vector<int> sortIndex;
	vector<double> FF;
	vector<double> CR;
	vector<double> SSFF;  // successful F values
	vector<double> SSCR;  // successful CR values 
	vector<int> binaryVector;
	vector< tFitness > parentsFitness;///array containing the current fitness of all particles 	
	vector< tFitness > offspringsFitness;	
	tFitness bestFitness;
	unsigned indexOfBest;
	unsigned numberOfIndividuals;
	vector< double > xp;
	uniform_real_distribution<double> unifRandom;	
	JADE *mainCollaborator;
};

