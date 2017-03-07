//
//  conslib.cpp
//  RedLibs Database Generator for fully degenerate libraries
//
//  Updated on 02/03/2017
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
//
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstring>

using namespace std;

/////////////////////////////
//	Function Declarations  //
/////////////////////////////

int TotalDegeneracy(long int n, int length);
string BaseSequence(string inputSequence, long int n, int length);

////////////////////
//	Main Program  //
////////////////////

int main ( int argc, char* argv[] ){
    // Generate Files with size and sequence of all possible consensus combinations based on a n bases long fully degenerate ("N") sequence

    // Check input
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " 4" << std::endl;
        exit ( EXIT_FAILURE );
    }
    
    // Global variables
    char databaseFileName[100];
    std::stringstream err;
    int length = atoi(argv[1]);
    long int maxN = 15;
    for (int i = 1; i < length; i++) maxN = maxN*15;
    
    // Generate output files
    for (long int n = 0; n < maxN; n++){
        sprintf(databaseFileName, "conslib-%d-%d.txt", length, TotalDegeneracy(n, length)); // Choose directory for database (directory needs to exist)
        ofstream outputFile;
        outputFile.open (databaseFileName, std::ios::app);
        outputFile << BaseSequence(argv[1], n, length) << "\n";
        outputFile.flush();
        outputFile.close();
        
    }
    return 0;
}

int TotalDegeneracy(long int n, int length) {
    int radix   = 15;
    int digit[] = {1,1,1,1,2,2,2,2,2,2,3,3,3,3,4};
    int factor;
    long int temp    = n;
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

string BaseSequence(string inputSequence, long int n, int length) {
    string sequence;
    string baseSymbol = "ACGTYRWSKMBDHVN";
    int radix = 15;
    int position = length - 1;
    long int temp = n;
    int remainder;
    while (temp >= 0 && position >= 0) {
        remainder = temp % radix;
        sequence  = baseSymbol[remainder] + sequence;
        temp      = ( temp - remainder ) / radix;
        position--;
    }
    return sequence;
}