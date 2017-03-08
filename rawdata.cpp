//
//  rawdata.cpp
//  RedLibsMPI
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
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
#include <algorithm>

RAWDAT::RAWDAT ( int* argc, char* argv[], int numnodes, int mynode ) {
    int count = *argc - 1;
    std::stringstream err;                              //error output string
    std::stringstream log;                              //log string to be written into the lof file
    
    int check[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};     //array for checking which input variables were given
                                                        // check: projectName, output, distribution, startPosition, endPosition, degeneracy, minLevel, maxLevel, mu, sigma
    char firstop = '\0';                                //indicates short form of input handle
    int secondop = 1;                                   //indicates long from of input handle
    
    if ( *argc < 3 ) {
        std::cout << "Error: No path to degenerate library or 'conslib' database provided.\n\n";
        help();
        exit ( EXIT_FAILURE );
    }
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Read variable from input variables and keep track which were not given */
    dataPath = argv[*argc - 1];
    dataBase = argv[*argc - 2];

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
                    op >> projectName;
                    check[0] = 1;
                    break;
                case 'o':
                    secondop = 0;
                    op >> output;
                    struct stat sb;
                    if ( output.at(output.length()-1) != '/' ) output = output + "/";
                    if ( stat ( output.c_str(), &sb ) != 0 ) 
					{
						if (mynode == 0)
						{
							//check whether output folder exists and create it if it does not exist
#ifdef WIN32
							system(("mkdir " + output.substr(0, output.length() - 1)).c_str());
#else
							system(("mkdir " + output).c_str());
#endif
						}
                    }
                    check[1] = 1;
                    break;
                case 'd':
                    secondop = 0;
                    if ( op.str() == "uniform" ) {
                        op >> distribution;
                        check[2] = 1;
                    } else if ( op.str() == "normal" ){
                        op >> distribution;
                        check[2] = 1;
                    } else {
                        err << "Error: illegal option argument (" << op.str() << ").\n\n";
                        std::cout << err.str();
                        if ( mynode == 0 ) writelog ( output, err.str() );
                        err.str("");
                        exit ( EXIT_FAILURE );
                    }
                    break;
                case 's':
                    secondop = 0;
                    op >> startPosition;
                    if( check[4] == 1 ){
                        length = endPosition - startPosition + 1;
                        if ( length <= 0 ) {
                            err << "Error: End of degenerated sequence smaller than its start (" << endPosition << " < " << startPosition << ").\n\n";
                            std::cout << err.str();
                            if ( mynode == 0 ) writelog ( output, err.str() );
                            err.str("");
                            exit ( EXIT_FAILURE );
                        }
                    }
                    check[3] = 1;
                    break;
                case 'e':
                    secondop = 0;
                    op >> endPosition;
                    if( check[3] == 1 ){
                        length = endPosition - startPosition + 1;
                        if ( length <= 0 ) {
                            err << "Error: End of degenerated sequence smaller than its start (" << endPosition << " < " << startPosition << ").\n\n";
                            std::cout << err.str();
                            if ( mynode == 0 ) writelog ( output, err.str() );
                            err.str("");
                            exit ( EXIT_FAILURE );
                        }
                    }
                    check[4] = 1;
                    break;
                case 'g': {
                    secondop = 0;
                    op >> degeneracy;
                    check[5] = 1;
                    if ( degeneracy < 2 ){
                        err << "Error: illegal option argument (" << op.str() << "). Degeneracy must be at least 2.";
                        std::cout << err.str();
                        degeneracy = 2;
                        check[5] = 1;
                        if ( mynode == 0 ) writelog ( output, err.str() );
                        err.str("");
                        exit ( EXIT_FAILURE );
                    }
                    std::stringstream conslibFilePathTest;
                    conslibFilePathTest << dataBase << "/conslib-" << length << "-" << degeneracy << ".txt";
                    std::ifstream conslibFile(conslibFilePathTest.str());
                    if ( !conslibFile ) {
                        std::stringstream availableDegeneracies;
                        std::stringstream conslibFilePathTemp;
                        int maxDegeneracy = 4;
                        for( int i; i < length-1; i++) maxDegeneracy = maxDegeneracy * 4;
                        for ( int i; i <= maxDegeneracy; i++ ){
                            conslibFilePathTemp << dataBase << "/conslib-" << length << "-" << i << ".txt";
                            std::ifstream conslibFileTemp(conslibFilePathTemp.str());
                            if (conslibFileTemp) {
                                availableDegeneracies << i << "\t";
                            }
                            conslibFilePathTemp.str("");
                            conslibFileTemp.close();
                        }
                        err << "Error: degeneracy (" << op.str() << ") not available or 'conslib' database file '" << dataBase << "/conslib-" << length << "-" << degeneracy << ".txt' was not created.\n\n" << "Available degeneracies are:\n" << availableDegeneracies.str() << "\n";
                        conslibFilePathTest.str("");
                        conslibFile.close();
                        std::cout << err.str();
                        if ( mynode == 0 ) writelog ( output, err.str() );
                        err.str("");
                        exit ( EXIT_FAILURE );
                    }
                    break;
                }
                case 'l':
                    secondop = 0;
                    op >> minLevel;
                    check[6] = 1;
                    break;
                case 'L':
                    secondop = 0;
                    op >> maxLevel;
                    check[7] = 1;
                    break;
                case 'u':
                    secondop = 0;
                    op >> mu;
                    check[8] = 1;
                    break;
                case 'a':
                    secondop = 0;
                    op >> sigma;
                    check[9] = 1;
                    break;
                case 'h':
                    secondop = 0;
                    help();
                    break;
                case '-':
                         if ( ss.str() == "--output" )     firstop = 'o';
                    else if ( ss.str() == "--projectName" )       firstop = 'n';
                    else if ( ss.str() == "--distribution" )      firstop = 'd';
                    else if ( ss.str() == "--startPosition" )   firstop = 's';
                    else if ( ss.str() == "--endPosition" )     firstop = 'e';
                    else if ( ss.str() == "--degeneracy" )    firstop = 'g';
                    else if ( ss.str() == "--minLevel" )     firstop = 'l';
                    else if ( ss.str() == "--maxLevel" )     firstop = 'L';
                    else if ( ss.str() == "--mu" )         firstop = 'u';
                    else if ( ss.str() == "--sigma" )      firstop = 'a';
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
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Set all variables that were not set with input options to default values */
    if ( check[0] == 0 ) {
        int substrinBeginning = 0;
        int substringEnd = 0;
        for ( int i = int(dataPath.length()-1); i > 0; i-- ) {
            if ( dataPath.at(i) == '.' ) substringEnd = i;
            if ( dataPath.at(i) == '/' ){
                substrinBeginning = i + 1;
                i = 0;
            }
        }
        projectName = dataPath.substr(substrinBeginning, substringEnd - substrinBeginning);
    }
    if ( check[1] == 0 ) {
        char timestring[20];
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        strftime ( timestring, 20, "%Y%m%d_%H%M_", timeinfo );
        
        output = timestring + projectName + "/";
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
        if( output.at(output.length()-1) != '/' ) output = output + "/";
		if (mynode == 0)
		{
#ifdef WIN32
			system(("mkdir " + output.substr(0, output.length() - 1)).c_str());
#else
			system(("mkdir " + output).c_str());
#endif
		}
    }
    if ( check[2] == 0 ) {
        distribution = "uniform";
    }
    if ( check[3] == 0 ) {
        startPosition = 1;
    }
    if ( check[4] == 0 ) {
        err << "Error: end position of degenerate sequence in data file has to be specified.\n\n";
        std::cout << err.str();
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
    if ( check[5] == 0 ) {
        degeneracy = 18;
    }
    if ( check[6] == 0 ) {
        std::ifstream indata ( dataPath.c_str() );
        std::string line;
        std::string dumpseq;
        
        if ( indata.is_open() ) {
            while ( indata.good() ) {
                std::stringstream readout;
                std::getline ( indata, line, ',' );
                readout.str(line);
                readout >> minLevel;
            }
            indata.close();
        }
        else {
            err << "Error: unable to open \"" << dataPath << "\".\n\n";
            std::cout << err.str();
            if ( mynode == 0 ) writelog ( output, err.str() );
            err.str("");
            exit ( EXIT_FAILURE );
        }
    }
    if ( check[7] == 0 ) {
        std::ifstream indata ( dataPath.c_str() );
        std::string line;
        std::string element;
        std::stringstream readout;
        std::stringstream readoutElement;
        
        if ( indata.is_open() ) {
            std::getline( indata, line );
            if( line.back() == '\r' ){
                indata.close();
                indata.open(dataPath.c_str());
                std::getline( indata, line, '\r' );
            }
            readout.str(line);
            while( readout.good() ) {
                std::getline( readout, element, ',' );
            }
            readoutElement.str(element);
            readoutElement >> maxLevel;
            indata.close();
        }
        else {
            err << "Error: unable to open \"" << dataPath << "\".\n\n";
            std::cout << err.str();
            if ( mynode == 0 ) writelog ( output, err.str() );
            err.str("");
            exit ( EXIT_FAILURE );
        }
    }
    if ( check[8] == 0 ) {
        mu = 0.0;
    }
    if ( check[9] == 0 ) {
        sigma = 1.0;
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Read in sequences and their corresponding data level
       while cropping the sequence to the degenerated part                   */
    std::ifstream indata ( dataPath.c_str() );
    std::string line;
    std::string tempSequence;
    std::string tempLevelString;
    std::stringstream tempLevelStringStream;
    double tempLevel;
    
    if ( indata.is_open() ){
        while ( indata.good() ) {
            std::stringstream readout;
            std::getline ( indata, line );
            if( line.back() == '\r' ) line.pop_back();
            readout.str(line);
            while( readout.good() ) {
                std::getline( readout, tempSequence, ',' );
                std::getline( readout, tempLevelString, ',' );
            }
            tempLevelStringStream.str(tempLevelString);
            tempLevelStringStream >> tempLevel;
            tempLevelStringStream.clear();
            if( int(tempSequence.length()) >= length ) tempSequence = tempSequence.substr( startPosition-1, length );
            data.push_back( sequenceData ( tempLevel, tempSequence ) );
            readout.str("");
        }
        indata.close();
    } else {
        err << "Error: unable to open \"" << dataPath << "\".\n\n";
        std::cout << err.str();
        if ( mynode == 0 ) writelog ( output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Get number of combinations of degenerate sequences with the chosen degeneracy */
    std::stringstream conslibFilePath;
    conslibFilePath << dataBase << "/conslib-" << length << "-" << degeneracy << ".txt";
    
    std::ifstream conslibFile(conslibFilePath.str());
    combinations = std::count(std::istreambuf_iterator<char>(conslibFile), std::istreambuf_iterator<char>(), '\n');
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Write summary of values used to calculate the reduced libraries into log */
    if ( mynode == 0 ) {
    log << "Project preferences:";
    log << "\n\tProject name:\t\t\t\t" << projectName;
    log << "\n\tData:\t\t\t\t\t" << dataPath;
    log << "\n\tOutput path:\t\t\t\t" << output;
    log << "\n\tDegenerated sequence range:\t" << startPosition << " - " << endPosition;
    log << "\n\tTotal degeneracy:\t\t\t" << degeneracy;
    log << "\n\t\t->\tNumber of consensus sequences to be checked:\n\t\t\t" << combinations;
    unsigned int total = 15;
    unsigned int total_base = 4;
    for ( int i = 1; i < length; ++i ) {
        total = total * 15;
        total_base = total_base * 4;
    }
    if ( distribution == "uniform" ) {
        log << "\n\tUniform target distribution between:\t" << minLevel << " - " << maxLevel;
    }
    if ( distribution == "normal" ) {
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
    std::string help = "RedLibs  version 1.1.0\nCopyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.\nRedLibs is an algorithm that allows for the rational design of smart libraries for pathway optimization thereby minimizing the use of experimental resources.\n\nUsage: RedLibs.o [OPTIONS]... [path of conslib database] [path of degenerate input library file]\n\nExample: RedLibs.o -g 18 -s 1 -e 6 /Volumes/user/RedLibs/conslib /Volumes/user/RedLibs/results/myProject\n\nOptions\n-o, --output\t\tpath for results\n\t\t\t\t\t\tdefault: [working directory]/[date]_[project name]\n-n, --projectName\tname of project\n\t\t\t\t\t\tdefault: name of degenerate library file\n-d, --distribution\ttype of target distribution: 'uniform' or 'normal'\n\t\t\t\t\t\tdefault: uniform\n-s, --startPosition\tstart position of degenerate sequence in raw data\n\t\t\t\t\t\tdefault value: 1\n-e, --endPosition\tend position of degenerate sequence in raw data has to be specified\n-g, --degeneracy\ttotal degeneracy\n\t\t\t\t\t\tdefault value: 2\n-l, --minLevel\t\tminimal sequence data level\n\t\t\t\t\t\tdefault value: last value in input library file (sorted ascending)\n-L, --maxLevel\t\tmaximal sequence data level\n\t\t\t\t\t\tdefault value: first value in input library file (sorted ascending)\n-h, --help\t\t\tshow help\n\nSee http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.\n";
    std::cout << help;
    exit ( EXIT_SUCCESS );
}
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */