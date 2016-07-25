//
//  rawdata.cpp
//  RedLibsMPI
//
//  Updated on 25/07/2016
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2015. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016
//

#include "rawdata.h"
#include "logfile.h"
#include <sstream>
#include <string>
#include <time.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/stat.h>

RAWDAT::RAWDAT ( int* argc, char* argv[], int numnodes, int mynode ) {
    int count = *argc - 1;
    std::stringstream err;                              //error output string
    std::stringstream log;                              //log string to be written into the lof file
    
    int check[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  //array for checking which input variables were given
                                                        // check: name, output, distr, startpos, endpos, minsize, maxsize, minlev, maxlev, mu, sigma
    char firstop = '\0';                                //indicates short form of input handle
    int secondop = 1;                                   //indicates long from of input handle
    
    if ( *argc < 2 ) {
        std::cout << "Error: No data file provided.\n\n";
        help();
        exit ( EXIT_FAILURE );
    }
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Read variable from input variables and keep track which were not given */
    while ( count > 0 ) {
        std::stringstream ss;
        std::stringstream op;
        secondop = 1;
        ss.str("");
        op.str("");
        if ( argv[count][0]=='-' ) {            //initially assume short form of input handle
            ss << argv[count];
            if ( count < *argc - 1 ) op << argv[count+1];
            else op << "empty";
            if ( op.str()[0] == '-' ) {
                err << "Error: illegal or missing option argument \"" << op.str() << "\" for option \"" << ss.str() << "\".\n\n";
                std::cout << err.str();
                if ( mynode == 0 ) writelog ( output, err.str() );
                err.str("");
                exit ( EXIT_FAILURE );
            }
            firstop = ss.str()[1];
            while ( secondop == 1 ) {
            switch ( firstop ) {
                case 'n':
                    secondop = 0;
                    op >> name;
                    check[0] = 1;
                    break;
                case 'o': {
                    secondop = 0;
                    op >> output;
                    struct stat sb;
                    if ( output.at(output.length()-1) != '/' ) output = output + "/";
                    if ( stat ( output.c_str(), &sb ) != 0 ) {
                        system ( ("mkdir " + output).c_str() );
                    }                                              //check whether output folder exists and create it if it does not exist
                    check[1] = 1;
                    break;
                }
                case 'd':
                    secondop = 0;
                    if ( op.str() == "uniform" ) {
                        op >> distr;
                        check[2] = 1;
                    }
                    else if ( op.str() == "normal" ){
                        op >> distr;
                        check[2] = 1;
                    }
                    else{
                        err << "Error: illegal option argument (" << op.str() << ").\n\n";
                        std::cout << err.str();
                        if ( mynode == 0 ) writelog ( output, err.str() );
                        err.str("");
                        exit ( EXIT_FAILURE );
                    }
                    break;
                case 's':
                    secondop = 0;
                    op >> startpos;
                    check[3] = 1;
                    break;
                case 'e':
                    secondop = 0;
                    op >> endpos;
                    check[4] = 1;
                    break;
                case 'm':
                    secondop = 0;
                    op >> minsize;
                    check[5] = 1;
                    if ( minsize < 2 ){
                        err << "Warning: illegal option argument (" << op.str() << "). Degeneracy must be at least 2.\nMinimal degeneracy will be set to 2.\n";
                        std::cout << err.str();
                        minsize = 2;
                        check[5] = 1;
                        if ( mynode == 0 ) writelog ( output, err.str() );
                        err.str("");
                    }
                    break;
                case 'M':
                    secondop = 0;
                    op >> maxsize;
                    check[6] = 1;
                    break;
                case 'l':
                    secondop = 0;
                    op >> minlev;
                    check[7] = 1;
                    break;
                case 'L':
                    secondop = 0;
                    op >> maxlev;
                    check[8] = 1;
                    break;
                case 'u':
                    secondop = 0;
                    op >> mu;
                    check[9] = 1;
                    break;
                case 'g':
                    secondop = 0;
                    op >> sigma;
                    check[10] = 1;
                    break;
                case 'h':
                    secondop = 0;
                    help();
                    break;
                case '-':
                         if ( ss.str() == "--output" )     firstop = 'o';
                    else if ( ss.str() == "--name" )       firstop = 'n';
                    else if ( ss.str() == "--distr" )      firstop = 'd';
                    else if ( ss.str() == "--startpos" )   firstop = 's';
                    else if ( ss.str() == "--endpos" )     firstop = 'e';
                    else if ( ss.str() == "--minsize" )    firstop = 'm';
                    else if ( ss.str() == "--maxsize" )    firstop = 'M';
                    else if ( ss.str() == "--minlev" )     firstop = 'l';
                    else if ( ss.str() == "--maxlev" )     firstop = 'L';
                    else if ( ss.str() == "--mu" )         firstop = 'u';
                    else if ( ss.str() == "--sigma" )      firstop = 's';
                    else if ( ss.str() == "--help" )       firstop = 'h';
                    else {
                        err << "Error: illegal option (" << ss.str() << ").\n\n";
                        std::cout << err.str();
                        if ( mynode == 0 ) writelog ( output, err.str() );
                        err.str("");
                        exit ( EXIT_FAILURE );
                    }
                default:
                    secondop = 0;
                    err << "Error: illegal option (" << ss.str() << ").\n\n";
                    std::cout << err.str();
                    if ( mynode == 0 ) writelog ( output, err.str() );
                    err.str("");
                    exit ( EXIT_FAILURE );
                    break;
            }
            }
        }
        --count;
    }
    datapath = argv[*argc - 1];
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Set all variables that were not set with input options to default values */
    if ( check[0] == 0 ) {
        int sub_beg = 0;
        int sub_end = int ( datapath.length() - 1 );
        for ( unsigned int i = 0; i < datapath.length(); ++i ) {
            if ( datapath.at(i) == '/' ) sub_beg = i + 1;
            if ( datapath.at(i) == '.' ) sub_end = i;
        }
        name = datapath.substr(sub_beg, sub_end);
    }
    if ( check[1] == 0 ) {
        char timestring[20];
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        strftime ( timestring, 20, "%Y%m%d_%H%M_", timeinfo );
        
        output = timestring + name + "/";
        struct stat sb;
        int foldersuffix = 0;
        if ( output.at(output.length()-1) != '/' ) output = output + "/";
        while ( stat ( output.c_str(), &sb ) == 0 ) {
            std::stringstream append;
            append << foldersuffix;
            output = output.substr( 0, output.size()-1 );
            output = output + append.str();
            ++foldersuffix;
            append.str("");
        }
        if ( output.at(output.length()-1) != '/' ) output = output + "/";
        system ( ("mkdir " + output).c_str() );
    }
    if ( check[2] == 0 ) {
        distr = "uniform";
    }
    if ( check[3] == 0 ) {
        startpos = 1;
    }
    if ( check[4] == 0 ) {
        err << "Error: end position of degenerate sequence in data file has to be specified.\n\n";
        std::cout << err.str();
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
    if ( check[5] == 0 ) {
        minsize = 2;
    }
    if ( check[6] == 0 ) {
        maxsize = 4;
        for ( int i = 0; i < endpos - startpos + 1; ++i ) maxsize = maxsize * 4;
        err << "Warning: Missing maximal library size.\nMaximal library size will be set to maximal possible library size.\n";
        std::cout << err.str();
    }
    if ( maxsize < minsize ){
        err << "Warning: illegal option argument (" << maxsize << "). Maximal library size must be larger or equal to minimal library size.\nMaximal library size will be set equal to minimal library size.\n";
        std::cout << err.str();
        maxsize = minsize;
        check[6] = 1;
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
    }
    if ( check[7] == 0 ) {
        std::ifstream indata ( datapath.c_str() );
        std::string line;
        std::string dumpseq;
        double lastlev;
        
        if ( indata.is_open() ) {
            while ( indata.good() ) {
                std::stringstream readout;
                std::getline ( indata, line );
                readout.str(line);
                readout >> dumpseq >> lastlev;
                if ( dumpseq.length() > 0 ) minlev = lastlev;
            }
            indata.close();
        }
        else {
            err << "Error: unable to open \"" << datapath << "\".\n\n";
            std::cout << err.str();
            if ( mynode == 0 ) writelog ( output, err.str() );
            err.str("");
            exit ( EXIT_FAILURE );
        }
    }
    if ( check[8] == 0 ) {
        std::ifstream indata ( datapath.c_str() );
        std::string line;
        std::stringstream readout;
        std::string dumpseq;
        
        if ( indata.is_open() ) {
            std::getline ( indata, line );
            readout.str(line);
            readout >> dumpseq >> maxlev;
            indata.close();
        }
        else {
            err << "Error: unable to open \"" << datapath << "\".\n\n";
            std::cout << err.str();
            if ( mynode == 0 ) writelog ( output, err.str() );
            err.str("");
            exit ( EXIT_FAILURE );
        }
    }
    if ( check[9] == 0 ) {
        mu = 0;
    }
    if ( check[10] == 0 ) {
        sigma = 1;
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Set length of degenerate sequence according to start and end in
       given data set                                                    */
    length = endpos - startpos + 1;
    if ( length <= 0 ) {
        err << "Error: End of degenerated sequence smaller than its start (" << endpos << " < " << startpos << ").\n\n";
        std::cout << err.str();
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Read in sequences and their corresponding data level
       while cropping the sequence to the degenerated part                   */
    std::ifstream indata ( datapath.c_str() );
    std::string line;
    std::string temp_seq;
    double temp_level;
    if ( indata.is_open() ){
        while ( indata.good() ) {
            std::stringstream readout;
            std::getline ( indata, line );
            readout.str(line);
            readout >> temp_seq >> temp_level;
            if( int(temp_seq.length()) >= length ) temp_seq = temp_seq.substr( startpos-1, length );
            data.push_back( seq_data ( temp_level, temp_seq ) );
            readout.str("");
        }
        indata.close();
    } else {
        err << "Error: unable to open \"" << datapath << "\".\n\n";
        std::cout << err.str();
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Get list of possible library sizes within given range */
    std::stringstream sizespath;
    
/*!!INPUT NEEDED!!*/
    sizespath << "/links/grid/shared/gerdanie/conslib0/_conslib-" << length << "-combinations.txt";     //the data base has to be placed into a folder available for programs running on the cluster. The path to this folder needs to be written here. In this example "/links/grid/shared/gerdanie/conslib0/"
    std::ifstream sizesdata ( sizespath.str().c_str() );
    unsigned int sum_combs = 0;             //sum of number of combinations of degenerate sequences with the chosen degeneracies
    int temp_size;
    unsigned int temp_combs;
    line = "";
    if ( sizesdata.is_open() ) {
        while ( sizesdata.good() ){
            std::stringstream readout;
            std::getline ( sizesdata, line );
            readout << line;
            readout >> temp_size >> temp_combs;
            if ( temp_size >= minsize && temp_size <= maxsize ){
                libsizes.push_back( cons_data ( temp_size, temp_combs ) );
                sum_combs = sum_combs + temp_combs;
            }
            readout.str("");
        }
        sizesdata.close();
        if ( libsizes.size() == 0 ){
            err << "Error: no library size in selected range " << "(" << minsize << " - " << maxsize << ").\n\n";
            std::cout << err.str();
            if ( mynode == 0 ) writelog ( output, err.str() );
            err.str("");
            exit ( EXIT_FAILURE );
        }
        libsizes.push_back( cons_data ( 0, sum_combs ) );
    } else {
        err << "Error: unable to open \"" << sizespath.str() << "\".\n\n";
        std::cout << err.str();
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Write summary of values used to calculate the reduced libraries into log */
    if ( mynode == 0 ) {
    log << "Project preferences:";
    log << "\n\tProject name:\t\t\t\t" << name;
    log << "\n\tData:\t\t\t\t\t" << datapath;
    log << "\n\tOutput path:\t\t\t\t" << output;
    log << "\n\tDegenerated sequence range:\t" << startpos << " - " << endpos;
    log << "\n\tLibrary size range:\t\t\t" << minsize << " - " << maxsize;
    log << "\n\t\t->\tNumber of consensus sequences to be checked:\n\t\t\t" << libsizes.at(libsizes.size()-1).combs;
    unsigned int total = 15;
    unsigned int total_base = 4;
    for ( int i = 1; i < length; ++i ) {
        total = total * 15;
        total_base = total_base * 4;
    }
    log << "\n\t\t\t(" << int(double(sum_combs) / double(total - total_base) * 100 + 0.5) << " % of all possible combinations)";
    if ( distr == "uniform" ) {
        log << "\n\tUniform target distribution between:\t" << minlev << " - " << maxlev;
    }
    if ( distr == "normal" ) {
        log << "\n\tNormal target distribution with:\nmean = \t" << mu << "\nstandard deviation =\t" << sigma;
    }
    log << "\n\n";
    
    writelog ( output, log.str(), 0 );
    log.str("");
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* Function printing the help text */
void help() {
    std::string help = "[conslib]  version 0.2.0\nCopyright (C) 2015 by Daniel Gerngross\n[conslib] is a program to...\n\nUsage: conslib [OPTION]... RAWDATAPATH\n\nOptions\n-o, --output\tpath for results\n-n, --name\tname of project\n-s, --startpos\tstart position of degenerated sequence in raw data\n\t\tdefault value: 1\n-e, --endpos\tend position of degenerated sequence in raw data\n\t\thas to be specified\n-m, --minsize\tminmal library size\n-M, --maxsize\tmaximal library size\n-l, --minlev\tminimal sequence data level\n-L, --maxlev\tmaximal sequence data level\n-h, --help\tshow this help\n\nSee http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.\n";
    std::cout << help;
    exit ( EXIT_SUCCESS );
}
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */