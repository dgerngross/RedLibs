//
//  rawdata.h
//  RedLibsMPI
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
//

#ifndef __conslibMPI__rawdata__
#define __conslibMPI__rawdata__

#include <iostream>
#include <vector>

struct sequenceData {
    sequenceData ( double level, std::string sequence ) : level(level), sequence(sequence) {}
    
    double level;                           //numeric value to which a specific DNA sequence is assigned
    std::string sequence;                   //corresponding DNA sequence
};

struct degenerateData {
    degenerateData ( int totalDegeneracy, unsigned int combinations ) : totalDegeneracy(totalDegeneracy), combinations(combinations) {}
    
    int totalDegeneracy;                    //total degeneracy of a degenerate sequence
    unsigned int combinations;              //number of degenerate sequences with a identical total degeneracy
};

class RAWDAT {
public:
    std::string dataPath;                   //path of the input file
    std::string dataBase;                   //path of the degenerate sequences database
    std::string output;                     //name of output folder
    std::string projectName;                //name of project
    std::string distribution;               //type of target distribution: 'uniform' and 'normal' available
    int startPosition;                      //start position of the degenerate sequence in the input sequence
    int endPosition;                        //end position of the degenerate sequence in the input sequence
    int length;                             //length of the degenerate sequence
    int degeneracy;                         //target total degeneracy
    long int combinations;                  //variable containing the number of combinations of degenerate sequences for the given total degeneracy
    double minLevel;                        //lower bound for uniform distribution
    double maxLevel;                        //upper bound for uniform distribution
    double mu;                              //mu of normal distribution
    double sigma;                           //sigma of normal distribution
    std::vector<sequenceData> data;         //variable containing input data

    RAWDAT ( int*, char*[], int numnodes, int mynode );
};

void help();                                //help function


#endif /* defined(__conslibMPI__rawdata__) */
