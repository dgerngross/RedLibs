//
//  main.cpp
//  RedLibsMPI
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
//

#include <iostream>
#include <string>
#include <sstream>
#include "rawdata.h"
#include "logfile.h"
#include "generate_list.h"
#include "mpi.h"                                        //MPI needs to be available on the cluster

using namespace std;

int main(int argc, char* argv[]){
    
    MPI_Init( &argc, &argv);
    int numnodes;                                       //number of nodes on the cluster
    int mynode;                                         //current node
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    
    RAWDAT rawdat ( &argc, argv, numnodes, mynode );    //load and prepare raw data
    
    generate_list ( rawdat, numnodes, mynode );   //generate reduced libraries
    
    
    if ( mynode == 0 ) {
        std::cout << "Calculations completed.\n\n";
        writelog ( rawdat.output, " Calculations completed.\n" );
    }
    MPI_Finalize();
    return 0;
}

