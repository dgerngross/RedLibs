//
//  logfile.h
//  RedLibsMPI
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
//

#ifndef __conslibMPI__logfile__
#define __conslibMPI__logfile__

#include <iostream>

void writelog ( std::string output, std::string message , int t = 1 );

#endif /* defined(__conslibMPI__logfile__) */
