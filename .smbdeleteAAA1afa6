//
//  rawdata.h
//  RedLibsMPI
//
//  Updated on 07/10/2015
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2015. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, [submission in progress]
//

#ifndef __conslibMPI__rawdata__
#define __conslibMPI__rawdata__

#include <iostream>
#include <vector>

struct seq_data {
    seq_data ( double level, std::string seq ) : level(level), seq(seq) {}
    
    double level;           //numeric value to which a specific DNA sequence is assigned
    std::string seq;        //corresponding DNA sequence
};

struct cons_data {
    cons_data ( int libsize, unsigned int combs ) : libsize(libsize), combs(combs) {}
    
    int libsize;            //total degeneracy of a degenerate sequence
    unsigned int combs;     //number of degenerate sequences with a identical total degeneracy
};

class RAWDAT {
public:
    std::string datapath;               //path of the input file
    std::string output;                 //name of output folder
    std::string name;                   //name of project
    std::string distr;                  //type of target distribution: 'uniform' and 'normal' available
    int startpos;                       //start position of the degenerate sequence in the input sequence
    int endpos;                         //end position of the degenerate sequence in the input sequence
    int length;                         //length of the degenerate sequence
    int minsize;                        //minimal target total degeneracy
    int maxsize;                        //maximal taget total degeneracy
    double minlev;                      //lower bound for uniform distribution
    double maxlev;                      //upper bound for uniform distribution
    double mu;                          //mu of normal distribution
    double sigma;                       //sigma of normal distribution
    std::vector<seq_data> data;         //variable containing input data
    std::vector<cons_data> libsizes;    //variable containing the number of combinations of degenerate sequences for each total degeneracy

    RAWDAT ( int*, char*[], int numnodes, int mynode );
};

void help();                            //help function


#endif /* defined(__conslibMPI__rawdata__) */
