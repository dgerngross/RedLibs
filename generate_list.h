//
//  generate_list.h
//  RedLibsMPI
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
//

#ifndef __conslibMPI__generate_list__
#define __conslibMPI__generate_list__

#include <iostream>
#include "rawdata.h"
#include "mpi.h"

struct dSequenceProperties {
    double distance;                //Kolmogorov-Smirnov distance to the target distribution
    std::string dSequence;          //degenerate sequence
    int degeneracy;                 //total degeneracy of the degenerate sequence
    
    dSequenceProperties () ;
    dSequenceProperties ( double d, std::string s, int l );
    dSequenceProperties ( const dSequenceProperties& );
    
    void set ( double d, std::string s, int l );
    void set ( const dSequenceProperties& c );
};

void generate_list ( RAWDAT rawdat, int numnodes, int mynode, MPI_Status Stat );
long int combrange ( int mynode, int numnodes, long int combs );
void getDSequenceProperties ( dSequenceProperties * tempDSequence, RAWDAT rawdat );
void sortBest ( std::vector<dSequenceProperties> * bestDSequence, dSequenceProperties * tempDSequence );

#endif /* defined(__conslibMPI__generate_list__) */
