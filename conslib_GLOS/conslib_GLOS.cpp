//
//  conslib_GLOS.cpp
//  RedLibs Database Generator for GLOS libraries
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite Sabine Oesterle, Daniel Gerngross, Steven Schmitt, Tania Michelle Roberts, and Sven Panke, [submission in progress], 2017
//
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

/////////////////////////////
//	Function Declarations  //
/////////////////////////////

int TotalDegeneracy(int n, int length);
string BaseSymbol(string inputSequence, int position);
string BaseSequence(string inputSequence, int n, int length);

////////////////////
//	Main Program  //
////////////////////

int main ( int argc, char* argv[] ){
    // Generate Files with size and sequence of all possible consensus combinations based on a 6 bases long partially degenerate sequence
    
    // Global variables
    char databaseFileName[100];
    std::stringstream err;
    int length = 6; // Length of degenerate sequence according to GLOS rules
    int maxn = 7;
    for (int i = 1; i < length; i++) maxn = maxn*7;
    
    // Check input
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " BDHVBD" << std::endl;
        exit ( EXIT_FAILURE );
    }
    
    if ( strlen(argv[1]) != length ) {
        err << "Error: the input sequence is not " << length << " bases long (don't use spaces within the sequence).\nExample: BDHVBD\n\n";
        std::cout << err.str();
        err.str("");
        exit ( EXIT_FAILURE );
    }
    
    // Generate output files
    for (int n = 0; n < maxn; n++){
        sprintf(databaseFileName, "conslib-%d-%d.txt", length, TotalDegeneracy(n, length)); // Choose directory for database (directory needs to exist)
        ofstream outputFile;
        outputFile.open (databaseFileName, std::ios::app);
        outputFile << BaseSequence(argv[1], n, length) << "\n";
        outputFile.flush();
        outputFile.close();

    }
    return 0;
}

int TotalDegeneracy(int n, int length) {
    int radix   = 7;
    int digit[] = {1,1,1,2,2,2,3};
    int factor;
    int temp    = n;
    int remainder;
    int product = 1;
    if(n == 0) return product;
    else {
        for(int i = 0; i < length; i++){
            remainder = temp % radix;
            factor    = digit[remainder];
            temp      = ( temp - remainder ) / radix;
            product   = product * factor;
        }
        return product;
    }
}

string BaseSymbol(string inputSequence, int position) {
    std::stringstream err;
    switch ( inputSequence[position] ) {
        case 'B':
            return "CGTSYKB";
            break;
        case 'D':
            return "AGTRWKD";
            break;
        case 'H':
            return "ACTMWYH";
            break;
        case 'V':
            return "ACGMRSV";
            break;
        default:
            err << "Error: the input sequence (" << inputSequence << ") does not consist of triple degenerate bases only, i.e. B, D, H, or V. Use only capital letters.\n\n";
            std::cout << err.str();
            err.str("");
            exit ( EXIT_FAILURE );
            break;
    }
}

string BaseSequence(string inputSequence, int n, int length) {
    string sequence;
    int radix = 7;
    int position = length - 1;
    int temp = n;
    int remainder;
    while (temp >= 0 && position >= 0) {
        remainder = temp % radix;
        sequence  = BaseSymbol(inputSequence, position)[remainder] + sequence;
        temp      = ( temp - remainder ) / radix;
        position--;
    }
    return sequence;
}