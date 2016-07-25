//
//  generate_list.h
//  RedLibsMPI
//
//  Updated on 25/07/2016
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2015. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016
//

#ifndef __conslibMPI__generate_list__
#define __conslibMPI__generate_list__

#include <iostream>
#include "rawdata.h"
#include <mpi.h>

struct cons_props {
    double distance;                //Kolmogorov-Smirnov distance to the target distribution
    std::string seq;                //degenerate sequence
    int libsize;                    //total degeneracy of the degenerate sequence
    
    cons_props () ;
    cons_props ( double d, std::string s, int l );
    cons_props ( const cons_props& );
    
    void set ( double d, std::string s, int l );
    void set ( const cons_props& c );
};

void generate_list ( RAWDAT rawdat, int numnodes, int mynode, MPI_Status Stat );
unsigned int combrange ( int mynode, int numnodes, unsigned int combs );
void get_cons_props ( cons_props * temp_cons, RAWDAT rawdat );
void sort_best ( std::vector<cons_props> * best_cons, cons_props * temp_cons );

#endif /* defined(__conslibMPI__generate_list__) */
