//============================================================================================
// SACCJADE - Surrogate-Assisted Cooperative Coevolutionary version of JADE
//============================================================================================
// Name        : main.cpp
// Authors     : Giuseppe A. Trunfio - trunfio@uniss.it
//               Ivanoe De Falco
//               Antonio Della Cioppa
// Version     : v1.0
// Created on  : Gen 16, 2017
//
//=============================================================================================

#define _USE_MATH_DEFINES
#include "CCDE.h"
#include "DG2.hpp"
#include "tclap/CmdLine.h"
#include "Benchmarks2010.h"
#include "header.h"

using namespace TCLAP;

using namespace std;


// create new object of class with default setting
Fitness* generateFuncObj(int funcID)
{
	Fitness *fp = NULL;

	using namespace CEC2010;
	if (funcID == 1)       fp = new F1();
	else if (funcID == 2)  fp = new F2();
	else if (funcID == 3)  fp = new F3();
	else if (funcID == 4)  fp = new F4();
	else if (funcID == 5)  fp = new F5();
	else if (funcID == 6)  fp = new F6();
	else if (funcID == 7)  fp = new F7();
	else if (funcID == 8)  fp = new F8();
	else if (funcID == 9)  fp = new F9();
	else if (funcID == 10) fp = new F10();
	else if (funcID == 11) fp = new F11();
	else if (funcID == 12) fp = new F12();
	else if (funcID == 13) fp = new F13();
	else if (funcID == 14) fp = new F14();
	else if (funcID == 15) fp = new F15();
	else if (funcID == 16) fp = new F16();
	else if (funcID == 17) fp = new F17();
	else if (funcID == 18) fp = new F18();
	else if (funcID == 19) fp = new F19();
	else if (funcID == 20) fp = new F20();
	else
	{
		cerr << "Fail to locate Specified Function Index" << endl;
		exit(-1);
	}

	return fp;
}

//******************************************************************************************/
//
//
//
//******************************************************************************************/
void optimization(int argc, char* argv[])
{   
	unsigned int functionIndex;
	unsigned int numRep;
	unsigned int numItePerCycle;
	unsigned int problemDimension;
	unsigned int sizeOfSubcomponents;
	unsigned int numOfIndividuals;
	unsigned int numberOfEvaluations;
    
	vector<int> seeds;
    unsigned maxNumRep = 100;
    for (unsigned i = 0; i < maxNumRep; ++i)
        seeds.push_back(i);

	//parsing command line
	try
	{
		CmdLine cmd("CCJADE - Cooperative Coevolution for Large-Scale Optimization of Computationally Expensive Objective Functions", ' ', "1.0");

		ValueArg<unsigned int> functionArg("f", "function", "function to optimize [1-20]", false, 12, "int");
		cmd.add(functionArg);

		ValueArg<unsigned int> repArg("r", "repetitions", "number of independent repetitions [1-100]", false, 1, "int");
		cmd.add(repArg);

		ValueArg<unsigned int> iteArg("i", "iterations", "number of JADE iterations per cycle", false, 50, "int");
		cmd.add(iteArg);		

		ValueArg<unsigned int> dimArg("d", "dimension", "problem dimension", false, 1000, "int");
		cmd.add(dimArg);		

		ValueArg<unsigned int> sdimArg("s", "subdim", "size of subcomponents", false, 250, "int");
		cmd.add(sdimArg);		

		ValueArg<unsigned int> npopArg("p", "popsize", "number of individuals in each subcomponent", false, 50, "int");
		cmd.add(npopArg);		

		ValueArg<unsigned int> feArg("e", "fevals", "allowed number of exact fitness evaluations", false, 3000000, "int");
		cmd.add(feArg);		

		cmd.parse(argc, argv);
		
		functionIndex = functionArg.getValue();
		numRep = repArg.getValue();
		numItePerCycle = iteArg.getValue();
		problemDimension = dimArg.getValue();
		sizeOfSubcomponents = sdimArg.getValue();
		numOfIndividuals = npopArg.getValue();
		numberOfEvaluations = feArg.getValue();

	}
    catch (ArgException& e)
    {
	   cout << "ERROR: " << e.error() << " " << e.argId() << endl;
    }

		
	if ( functionIndex < 1 || functionIndex>20 )
	{
		cerr << "function index out of allowed bounds [1..20]" << endl;
		exit(1);
	}


	if (problemDimension<1 || problemDimension>1000)
	{
		cerr << "problem dimension must be in [1..1000]" << endl;
		exit(1);
	}

	if (sizeOfSubcomponents<1 || sizeOfSubcomponents>problemDimension)
	{
		cerr << "problem dimension must be in [1..problem dimension]" << endl;
		exit(1);
	}



	//from: Marwa Sharawi and Mohammed El-Abd, A Cooperative Co-evolutionary LSHADE Algorithm for Large - Scale Global Optimization
	// ©2017 IEEE
	switch (functionIndex) 
	{
		case 1: sizeOfSubcomponents = 10; break;
		case 2: sizeOfSubcomponents = 10; break;
		case 3: sizeOfSubcomponents = 50; break;
		case 4: sizeOfSubcomponents = 100; break;
		case 5: sizeOfSubcomponents = 100; break;
		case 6: sizeOfSubcomponents = 50; break;
		case 7: sizeOfSubcomponents = 250; break;
		case 8: sizeOfSubcomponents = 250; break;
		case 9: sizeOfSubcomponents = 200; break;
		case 10: sizeOfSubcomponents = 10; break;
		case 11: sizeOfSubcomponents = 25; break;
		case 12: sizeOfSubcomponents = 100; break;
		case 13: sizeOfSubcomponents = 100; break;
		case 14: sizeOfSubcomponents = 50; break;
		case 15: sizeOfSubcomponents = 50; break;
		case 16: sizeOfSubcomponents = 50; break;
		case 17: sizeOfSubcomponents = 50; break;
		case 18: sizeOfSubcomponents = 50; break;
		case 19: numOfIndividuals = 50; break;
		case 20: numOfIndividuals = 50;  break;
	}



	//create the fitness function object
	Fitness* f = generateFuncObj(functionIndex);

    cout << "Optimizing f" << functionIndex  << endl;
	cout << "Problem dimension = " << problemDimension << endl;
    cout << "Number of iterations per cycle = " << numItePerCycle << endl;    	
	//cout << "The problem is decomposed in " << problemDimension / sizeOfSubcomponents << " equal subcomponents of size " << sizeOfSubcomponents << endl;
	cout << "The problem is decomposed using Differential Grouping" << endl;
	
	cout << "Number of repetitions = " << numRep << endl;
	cout << "Allowed number of exact function evaluations = " << numberOfEvaluations << endl;

	

	//create the decomposition using Differential Grouping
	DG2 dg2(f, problemDimension, -f->getMaxX(), f->getMaxX());
	dg2.ism();
	dg2.dsm();	
	cout << "FEs used in DG: " << dg2.getNumEvaluations() << endl;
	
	//find the maximum size of nonseparable subcomponent
	unsigned ms = 0;
	for (int i = 0; i < dg2.nonseps.size(); i++)
		if (ms < dg2.nonseps[i].size())
			ms = dg2.nonseps[i].size();
	ms = max(ms, sizeOfSubcomponents);
	numOfIndividuals = max(min((unsigned)100, ms), (unsigned)25);


	cout << "Number of individuals per subcomponent = " << numOfIndividuals << endl;

    double time = 0;
    vector< vector<ConvPlotPoint> > convergences;
    for (unsigned k = 0; k < numRep; ++k)
    {
        vector<ConvPlotPoint> convergence;
        CCDE ccde;
        int seed = seeds[k];
        ccde.optimize(f, problemDimension, f->getMaxX(), 0.0, numberOfEvaluations, sizeOfSubcomponents, numOfIndividuals,
                      convergence, seed, numItePerCycle, dg2);
        time = ccde.elapsedTime;
        convergences.push_back(convergence);
    }

    char fName[256];
    FILE *file;

    if (convergences.size() == 0) return;
    if (convergences[0].size() == 0) return;

	double avg = 0.0;
	double stdDev = 0.0;
	for (unsigned q = 0; q < convergences.size(); ++q)
		avg += convergences[q][convergences[q].size() - 1].f;
	avg /= convergences.size();
	for (unsigned q = 0; q < convergences.size(); ++q)
		stdDev += (avg - convergences[q][convergences[q].size() - 1].f)*(avg - convergences[q][convergences[q].size() - 1].f);
	stdDev /= convergences.size();
	stdDev = sqrt(stdDev);

	
	fopen_s(&file, "results.csv", "at");
	fprintf(file, "f%d; %.8lE; %.8lE\n", functionIndex, avg, stdDev);
	fclose(file);


    sprintf_s(fName, "convplot_f%d_dec%d_popsize%d_CCJADE.csv", functionIndex, sizeOfSubcomponents, numOfIndividuals);
    
    fopen_s(&file, fName, "wt");
    int maxSize = 0;
    int idOfMaxSize = 0;
    for (unsigned q = 0; q < convergences.size(); ++q)
        if (convergences[q].size() > maxSize)
        {
            maxSize = convergences[q].size();
            idOfMaxSize = q;
        }

    for (unsigned q = 0; q < maxSize; ++q)
    {
        for (unsigned k = 0; k < numRep; ++k)
        {
            if (q < convergences[k].size())
            {
                fprintf(file, "%d; %.8Le; ", convergences[k][q].nfe, convergences[k][q].f);
            }
            else
            {
                fprintf(file, "; ;");
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


//******************************************************************************************/
//
//
//
//******************************************************************************************/
int main(int argc, char* argv[])
{
    optimization(argc, argv);
}

