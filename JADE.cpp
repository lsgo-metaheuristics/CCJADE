//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : Jade.cpp
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
//=============================================================================================

#include "JADE.h"
#include "CCDE.h"
#include <chrono>

using namespace std;

//******************************************************************************************/
//
//
//
//******************************************************************************************/
JADE::JADE(unsigned _dimension, unsigned _numberOfIndividuals, Decomposer &_group) :
    decomposer(_group), dimension(_dimension), numberOfIndividuals(_numberOfIndividuals)
{
    JADE_mu_cr = 0.5;
    JADE_mu_ff = 0.5;
    nfe = 0;
    parentsFitness.resize(numberOfIndividuals, 0);    
    offspringsFitness.resize(numberOfIndividuals, 0);
    FF.resize(numberOfIndividuals, 0);
    CR.resize(numberOfIndividuals, 0);
    SSFF.resize(numberOfIndividuals, 0);
    SSCR.resize(numberOfIndividuals, 0);
    indexOfBest = (unsigned)(numberOfIndividuals * unifRandom(decomposer.eng));
    binaryVector.resize(dimension, 0);
    xp.resize(decomposer.CCOptimizer.problemDimension);    
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
JADE::~JADE()
{    
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
vector<double> &JADE::getCollaborator()
{
    return parents[indexOfBest];
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateContextVector()
{
    for (unsigned ld = 0; ld < coordinates.size(); ld++)
        decomposer.contextVector[coordinates[ld]] = parents[indexOfBest][ld];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateContextVector(vector<double> &cv, vector<unsigned> &coords, unsigned &vi)
{
    for (unsigned ld = 0; ld < coordinates.size(); ld++)
    {
        cv[vi + ld] = parents[indexOfBest][ld];
        coords[vi + ld] = coordinates[ld];
    }
    vi += coordinates.size();
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setCoordinates(vector<unsigned> &_coordinates)
{
    coordinates = _coordinates;
};


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setCoordinates(unsigned *_coordinates, unsigned numOfCoordinates)
{
    if (coordinates.size() == numOfCoordinates)
    {
        for (unsigned i = 0; i < numOfCoordinates; ++i)
            coordinates[i] = _coordinates[i];
        return;
    }

    coordinates.clear();
    for (unsigned i = 0; i < numOfCoordinates; ++i)
        coordinates.push_back(_coordinates[i]);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::loadIndividuals(vector< vector<double> > &population)
{
    parents.clear();

    for (unsigned i = 0; i < numberOfIndividuals && i < population.size(); ++i)
    {
        vector< double > position;

        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            position.push_back(population[i][coordinates[ld]]);

        parents.push_back(position);
    }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::updateIndividuals(vector< vector<double> > &population)
{
    for (unsigned i = 0; i < numberOfIndividuals && i < population.size(); ++i)
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            parents[i][ld] = population[i][coordinates[ld]];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::storeIndividuals(vector< vector<double> > &population)
{
    for (unsigned i = 0; i < numberOfIndividuals && i < population.size(); ++i)
        for (unsigned ld = 0; ld < coordinates.size(); ld++)
            population[i][coordinates[ld]] = parents[i][ld];
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::sortPopulation(vector<tFitness> &fitness, vector<int> &sortIndex)
{
    sortIndex.resize(fitness.size());

    for (unsigned int j = 0; j < fitness.size(); j++)
        sortIndex[j] = j;

    sort(&sortIndex[0], &sortIndex[0] + sortIndex.size(), doCompareIndividuals(&fitness[0]));
}




//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::update()
{
    //Sort the population from best to worst
    sortPopulation(parentsFitness, sortIndex);


    //Generate the CR and F values based on Gaussian and Cauchy distribution, respectively
    cauchy_distribution<double> cauchy(JADE_mu_ff, 0.1);

    normal_distribution<double> gaussian(JADE_mu_cr, 0.1);

    for (unsigned int i = 0; i < parents.size(); i++)
    {
        do
        {
            FF[i] = cauchy(decomposer.eng);
        } while (FF[i] <= 0.0);

        if (FF[i] > 1.0) FF[i] = 1.0;

        CR[i] = gaussian(decomposer.eng);

        if (CR[i] < 0.0) CR[i] = 0.0;

        if (CR[i] > 1.0) CR[i] = 1.0;
    }

    offsprings.clear();
    vector<bool> toEvaluate;
    for (unsigned int i = 0; i < parents.size(); i++)
    {
        unsigned r1, r2, r3;
        //Generate the mutant vector
        //Randomly choose the p_best individual
        unsigned p_index = unifRandom(decomposer.eng) * parents.size() * decomposer.CCOptimizer.JADE_p;
        p_index = sortIndex[p_index];

        //Select three parents randomly
        do
        {
            r1 = unifRandom(decomposer.eng) * parents.size();
        } while (r1 == i);

        do
        {
            r2 = unifRandom(decomposer.eng) * parents.size();
        } while (r2 == i || r2 == r1);

        do
        {
            r3 = unifRandom(decomposer.eng) * parents.size();
        } while (r3 == i || r3 == r2 || r3 == r1);

        vector<double> child(dimension, 0);

        for (unsigned int j = 0; j < dimension; j++)
        {
            if (decomposer.CCOptimizer.JADE_mutationStrategy == 1)
            {
                child[j] = parents[i][j] +
                           FF[i] * (parents[p_index][j] - parents[i][j]) +
                           FF[i] * (parents[r1][j] - parents[r2][j]);
            }
            else if (decomposer.CCOptimizer.JADE_mutationStrategy == 2)
            {
                child[j] = parents[r1][j] +
                           FF[i] * (parents[p_index][j] - parents[r1][j]) +
                           FF[i] * (parents[r2][j] - parents[r3][j]);
            }

            if (child[j] < decomposer.CCOptimizer.lowerLimit || child[j] > decomposer.CCOptimizer.upperLimit)
                child[j] = decomposer.CCOptimizer.lowerLimit + unifRandom(decomposer.eng) * (decomposer.CCOptimizer.upperLimit - decomposer.CCOptimizer.lowerLimit);
        }

        offsprings.push_back(child);
        toEvaluate.push_back(false);

        //Generate the binary vector based on the binomial crossover
        unsigned j_rnd = dimension * unifRandom(decomposer.eng);

        for (unsigned j = 0; j < dimension; j++)
        {
            if (unifRandom(decomposer.eng) < CR[i] || j == j_rnd)
                binaryVector[j] = 1;
            else
                binaryVector[j] = 0;
        }

        //Repair the crossover rate with its binary vector generated above
        unsigned int tt = 0;

        for (unsigned int j = 0; j < dimension; j++)
            tt += binaryVector[j];

        CR[i] = (double)tt / ((double)dimension);

        //Generate the trial vector based on the binary vector and the mutant vector
        for (unsigned int j = 0; j < dimension; j++)
            if (binaryVector[j] == 0)
                offsprings[i][j] = parents[i][j];
            else
                toEvaluate[i] = true;
    }

    //Evaluate the child population
    vector<bool> offspringHasTrueFitness(offsprings.size(), true);
    evaluateOffsprings(toEvaluate, offspringHasTrueFitness);


    //Selection and save the successful parameters
    SSFF.clear();
    SSCR.clear();

    for (unsigned int i = 0; i < parents.size(); i++)
    {
        if (offspringsFitness[i] <= parentsFitness[i])
        {
            for (unsigned int j = 0; j < dimension; j++)
                parents[i][j] = offsprings[i][j];

            parentsFitness[i] = offspringsFitness[i];

            //Save the successful CR and F values
            SSFF.push_back(FF[i]);
            SSCR.push_back(CR[i]);
        }
    }

    //Update mu_CR and mu_F based on the successful CRs and Fs
    if (SSCR.size() )
    {
        //Update the mu_CR values
        double mean_cr = 0.0;

        for (unsigned int i = 0; i < SSCR.size(); i++)
            mean_cr += SSCR[i];

        mean_cr = mean_cr / ((double)SSCR.size());
        JADE_mu_cr = (1 - decomposer.CCOptimizer.JADE_c) * JADE_mu_cr + decomposer.CCOptimizer.JADE_c * mean_cr;

        //Update the mu_F value
        double mean_ff = 0.0;
        double t1 = 0.0;
        double t2 = 0.0;

        for (unsigned int i = 0; i < SSFF.size(); i++)
        {
            t1 += SSFF[i] * SSFF[i];
            t2 += SSFF[i];
        }

        mean_ff = t1 / t2;

        //Lehmer mean
        JADE_mu_ff = (1 - decomposer.CCOptimizer.JADE_c) * JADE_mu_ff + decomposer.CCOptimizer.JADE_c * mean_ff;
    }
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::setParentFitness(vector<tFitness> &fitnessValues)
{
    parentsFitness = fitnessValues;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::evaluateOffsprings(vector<bool> &toEvaluate, vector<bool> &hasTrueFitness)
{
    //use real fitness function
    for (unsigned d = 0; d < decomposer.CCOptimizer.problemDimension; ++d)
        xp[d] = decomposer.contextVector[d];

    offspringsFitness.resize(offsprings.size());

    for (unsigned i = 0; i < offsprings.size(); i++)
    {
        if ( !toEvaluate[i] )
        {
            offspringsFitness[i] = parentsFitness[i];            
            continue;
        }

        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = offsprings[i][ld];

        offspringsFitness[i] = decomposer.CCOptimizer.computeFitnessValue(xp);
        nfe++;        

        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = decomposer.contextVector[coordinates[ld]];
    }
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
int JADE::evaluateParents()
{
    for (unsigned d = 0; d < decomposer.CCOptimizer.problemDimension; ++d)
        xp[d] = decomposer.contextVector[d];

    parentsFitness.resize(parents.size());

    for (unsigned i = 0; i < parents.size(); i++)
    {
        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = parents[i][ld];

        parentsFitness[i] = decomposer.CCOptimizer.computeFitnessValue(xp);
        nfe++;            

        for (unsigned ld = 0; ld < dimension; ld++)
            xp[coordinates[ld]] = decomposer.contextVector[coordinates[ld]];
    }

    updateIndexOfBest();

    return parents.size();
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness JADE::calculateFitnessValue(vector<double> &p, bool updateArchive)
{
    for (unsigned d = 0; d < decomposer.CCOptimizer.problemDimension; ++d)
        xp[d] = decomposer.contextVector[d];

    for (unsigned ld = 0; ld < coordinates.size(); ld++)
        xp[coordinates[ld]] = p[ld];

    tFitness f = decomposer.CCOptimizer.computeFitnessValue(xp);

    nfe++;
    return f;
}



//******************************************************************************************/
//
//
//
//******************************************************************************************/
tFitness JADE::updateIndexOfBest()
{
    bestFitness = std::numeric_limits<tFitness>::infinity();
    for (unsigned i = 0; i < parents.size(); ++i)
        if (parentsFitness[i] < bestFitness)
        {
            indexOfBest = i;
            bestFitness = parentsFitness[i];
        }
	return bestFitness;
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
void JADE::optimize(int iterations)
{
	tFitness prevFitness = updateIndexOfBest();
	for (ite = 0; ite < iterations; ++ite)
    {
        if ( decomposer.CCOptimizer.numberOfEvaluations + this->nfe >= decomposer.CCOptimizer.maxNumberOfEvaluations )
            break;
        update();
    }
	lastFitnessGain = fabs(prevFitness - updateIndexOfBest());
}
