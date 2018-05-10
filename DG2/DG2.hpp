//
//  DG.hpp
//  DG
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#ifndef DG2_hpp
#define DG2_hpp

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include "Fitness.h"


typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, float > Graph;

class DG2 {
    
    public:
    
        DG2(Fitness* f, unsigned _dim, tFitness _lb, tFitness _ub);
        ~DG2();
        void ism (void);
        void dsm (void);
        void dg (double);
        bool save (std::string);
        bool load (std::string);
        int getNumEvaluations();
		std::vector< std::vector<double> > getLambda();
		std::vector< std::vector<double> > getEInf();
		std::vector< std::vector<double> > getESup();
    
 
    
        Graph g;
		Fitness *func;
        double **delta1;
        double **delta2;
        double **lambda;
        double **theta;
        double **epsilon;
		double **eInfV;
		double **eSupV;
        int eval_count;
        double **archive_f;
        double *archive_fhat;
        double archive_base;
        std::vector<unsigned int> seps;
        std::vector< std::vector<unsigned int> > nonseps;
		unsigned dim;
		tFitness lb, ub;

        void arrayCopy (double *dest, double *src, int size);
        double gammaFunc (double d);
    
};
#endif /* DG_hpp */
