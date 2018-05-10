//
//  main.cpp
//  DG2
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Function.hpp"
#include "DG2.hpp"
#include "CEC2013.hpp"
#include <vector>
#include <random>
#include <numeric>
using namespace std;
///Pseudorandom generator
typedef mt19937 RandomEngine;
RandomEngine eng;

uniform_real<double> unifRandom;

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


void removeVariable(int v, list<Interactions> &allInteractions)
{
	allInteractions.erase(remove_if(allInteractions.begin(), allInteractions.end(), [v](Interactions& i) { return i.variable == v; }), allInteractions.end());
	
	list<Interactions>::iterator iterator = allInteractions.begin();
	for (; iterator != allInteractions.end(); ++iterator)
	{
		iterator->interactions.erase(remove_if(iterator->interactions.begin(), iterator->interactions.end(), [v](Interaction& i) { return i.variable == v; }), iterator->interactions.end());
	}
	
}


void normalizeSoftmax(double tau, list<Interaction> &I)
{
	double den = 0.0;
	
	for (list<Interaction>::iterator iterator = I.begin(), end = I.end(); iterator != end; ++iterator)
	{
		iterator->intensity = exp(iterator->intensity / tau);
		den += iterator->intensity;
	}
	transform(I.begin(), I.end(), I.begin(), [den](Interaction d) -> Interaction { d.intensity /= den; return d; });
	
}


int main(int argc, char **argv) {
    
    char filename[256];
    //for (int i = 1 ; i < argc ; i++){
      //  int funID = atoi(argv[i]);
	    int funID = 15;
        cout<<"Function: "<<funID<<endl;
        CEC2013 f(funID);
        DG2 dg2(&f);
        dg2.ism();
        dg2.dsm();
        sprintf(filename, "./results/theta%02d.txt", funID);
        dg2.save(filename);
        cout<<"FEs: "<<dg2.getNumEvaluations()<<endl;

    //}


    //parameters
	int sizeOfSubcomponents = 100;
	double tau = 0.1;
	
	list<Interactions> allInteractions;

	vector<vector<double>> lambda = dg2.getLambda();
	
	for (int i = 0; i < f.getDim(); ++i)
	{
		Interactions Ii;
		for (int j = 0; j < f.getDim(); ++j)
		  if( i!=j )
		  {
			  Ii.interactions.push_back(Interaction{ j, lambda[i][j] });
		  }
		Ii.variable = i;
		allInteractions.push_back(Ii);
	}
	

	
	for(list<Interactions>::iterator iterator=allInteractions.begin(); iterator!= allInteractions.end(); ++iterator)
	{
		list<Interaction>::iterator maxInteraction = max_element(iterator->interactions.begin(), iterator->interactions.end(),
			[](Interaction const& a, Interaction const& b) {return a.intensity < b.intensity; });
		double norm = maxInteraction->intensity;
		transform(iterator->interactions.begin(), iterator->interactions.end(), iterator->interactions.begin(), [norm](Interaction d) -> Interaction { d.intensity /= norm; return d; });
	}
	

	//creating subcomponents
	vector< vector<int> > decomposition;

	//random_shuffle(allInteractions.begin(), allInteractions.end());

	for (int i = 0; i < f.getDim()/ sizeOfSubcomponents; i++)
	{
		vector<int> variables;

		//pick the first variable
		int v = allInteractions.begin()->variable;
		list<Interaction> interactions = allInteractions.begin()->interactions;
		
		variables.push_back(v);

		removeVariable(v, allInteractions);
		
		//randomly select the next variables
		if (interactions.size() < sizeOfSubcomponents)
		{
			for (auto iterator = interactions.begin(); iterator != interactions.end(); ++iterator)
				variables.push_back(iterator->variable);
		}
		else
		for (int j = 0; j < sizeOfSubcomponents - 1; ++j)
		{
			normalizeSoftmax(tau, interactions);
			
			double p = unifRandom(eng);
			double cumProb = 0;
			
			list<Interaction>::iterator iterator = interactions.begin();
			for(; iterator != interactions.end(); ++iterator)
			{
				cumProb += iterator->intensity;
				if (cumProb >= p)
					break;
			}
			
			variables.push_back(iterator->variable);
			
			interactions.erase(iterator);
			
			removeVariable(iterator->variable, allInteractions);
		}

		decomposition.push_back(variables);

		cout << "Subcomponent" << i + 1 << endl;
		for (int k = 0; k < variables.size(); ++k)
			cout << variables[k] << " ";
		cout << endl;
	}

	
	
	
	//check decomposition
	/*
	vector<int> v1;
	for (int i = 0; i < f.getDim(); ++i)
		v1.push_back(i);

	vector<int> vv;
	for (int j = 0; j < decomposition.size(); j++)
		for (int k = 0; k < decomposition[j].size(); k++)
			vv.push_back(decomposition[j][k]);
	sort(vv.begin(), vv.end());
	cout << "size=" << vv.size() << endl;
	for (int i = 0; i < vv.size(); ++i)
		cout << vv[i] << " ";

	for (int i = 0; i < f.getDim(); ++i)
		if (vv[i] != v1[i])
			cout << "cazzo!" << endl;
	*/

	system("pause");
    return 0;
}
