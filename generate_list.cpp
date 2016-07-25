//
//  generate_list.cpp
//  RedLibsMPI
//
//  Updated on 25/07/2016
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2015. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016
//

#include "generate_list.h"
#include "logfile.h"
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <time.h>
#include <float.h>
#include <cmath>
#include <time.h>
#include <mpi.h>

void generate_list ( RAWDAT rawdat, int numnodes, int mynode, MPI_Status Stat ){
    
    if ( mynode == 0 ) writelog ( rawdat.output, " Starting calculations.\n" );
    
    int rc;                                         //MPI communication
    
    unsigned int my_combmin;                        //lower index bound of degenerate sequences assigned to current node
    unsigned int my_combmax;                        //upper index bound of degenerate sequences assigned to current node
    unsigned int my_combrange;                      //size of corresponding index range
    unsigned int totprogcount = 0;                  //used to calculate progress
    std::stringstream progresspath;                 //path to progress file
    progresspath << rawdat.output << rawdat.name << "_progress.txt";
    
    std::vector<cons_props> best_cons_tot;          //top ten list of degenerate sequences with lowest Kolmogorov-Smirnov distance
    for (int i = 0; i < 10; ++i ) best_cons_tot.push_back ( cons_props() );
    
    time_t * start_time;                            //used to estimate time left from progress
    start_time = new time_t[2];
        
    time_t * current_time;
    current_time = new time_t;
    
    double elapsed_time[2], remaining_time[2], progproc[2];
    int iremaining_t[2];
    std::string t_unit[2];
    
    if ( mynode == 0 ) start_time[1] = time(NULL);
    
    for ( int i_size = 0; i_size < ( int(rawdat.libsizes.size()) - 1 ); ++i_size) {
        
        my_combmin = combrange ( mynode, numnodes, rawdat.libsizes[i_size].combs );
        my_combmax = combrange ( mynode+1, numnodes, rawdat.libsizes[i_size].combs );
        my_combrange = my_combmax - my_combmin;
        
        std::vector<cons_props> best_cons;
        for (int i = 0; i < 10; ++i ) best_cons.push_back ( cons_props() );
        std::stringstream conspath;
/*!!INPUT NEEDED!!*/
        conspath << "/links/grid/shared/gerdanie/conslib0/conslib-" << rawdat.length << "-" << rawdat.libsizes[i_size].libsize << ".txt"; //the data base has to be placed into a folder available for programs running on the cluster. The path to this folder needs to be written here. In this example "/links/grid/shared/gerdanie/conslib0/"
        std::ifstream consdata ( conspath.str().c_str() );
        
        std::string * consline;
        consline = new std::string[my_combrange];
        int i_line = 0;
        int wait_msg = 0;
        unsigned int lineno = 0;
        for ( int i_node = 0; i_node < numnodes; ++i_node ) {
            wait_msg = 0;
            if ( mynode == i_node ) {
                i_line = 0;
                lineno = 0;
                while ( consdata.good() && lineno < my_combmax ) {
                    std::getline ( consdata, consline[i_line] );
                    ++lineno;
                    if ( lineno > my_combmin ) {
                        ++i_line;
                    }
                }
                wait_msg = 1;
                MPI_Bcast ( &wait_msg, 1, MPI_INT, i_node, MPI_COMM_WORLD);

            } else {
                while ( wait_msg == 0 && i_node < numnodes ) MPI_Bcast ( &wait_msg, 1, MPI_INT, i_node, MPI_COMM_WORLD);
            }
        }
        
        if ( mynode == 0 ) start_time[0] = time(NULL);
        
        for ( unsigned int i_comb = 0; i_comb < my_combrange; ++i_comb ) {

            cons_props * temp_cons;
            temp_cons = new cons_props;
            temp_cons->set ( temp_cons->distance, consline[i_comb], rawdat.libsizes[i_size].libsize );

            if ( int(temp_cons->seq.length()) == rawdat.length ) get_cons_props ( temp_cons, rawdat );          //generate distributions and calculate Kolmogorov-Smirnov distace, see below

            if ( best_cons.back().distance > temp_cons->distance ) sort_best ( &best_cons, temp_cons );
            
            
            /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
            /* Report progress (aprrox. from node 0 progress) */
            if ( ( (i_comb - 1) % 1000 == 0 || i_comb == (my_combrange - 1 ) ) && mynode == 0 ) {
                progproc[0] = 100.0 * double (i_comb) / double (my_combrange);
                progproc[1] = 100.0 * double (numnodes) * double (totprogcount) / double (rawdat.libsizes.back().combs);
                
                *current_time = time(NULL);
                for ( int i_t = 0; i_t < 2; ++i_t ) {
                    elapsed_time[i_t] = difftime(*current_time, start_time[i_t]);
                    remaining_time[i_t] = (100 / progproc[i_t] - 1) * (elapsed_time[i_t]);
                    
                    if ( remaining_time[i_t] < 60 ) {
                        iremaining_t[i_t] = int(remaining_time[i_t] + 0.5);
                        t_unit[i_t] = " seconds";
                    }
                    else if ( remaining_time[i_t] < 3600 ) {
                        iremaining_t[i_t] = int(remaining_time[i_t]/60 + 0.5);
                        t_unit[i_t] = " minutes";
                    }
                    else if ( remaining_time[i_t] < 86400 ) {
                        iremaining_t[i_t] = int(remaining_time[i_t]/3600 + 0.5);
                        t_unit[i_t] = " hours";
                    }
                    else {
                        iremaining_t[i_t] = int(remaining_time[i_t]/86400 + 0.5);
                        t_unit[i_t] = " days";
                    }
                }
                
                std::ofstream progress;
                progress.open ( progresspath.str().c_str() );
                progress << "Library size\t" << rawdat.libsizes[i_size].libsize << "\t";
                progress.precision(3);
                progress << progproc[0] << " %\t";
                progress << "Remaining time: " << iremaining_t[0] << t_unit[0] << "\n";
                progress << "All library sizes\t";
                progress << progproc[1] << " %\t";
                progress << "Remaining time: " << iremaining_t[1] << t_unit[1] << "\n";
                progress.close();
            }
            /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
            
            delete temp_cons;
            if ( mynode == 0 ) ++totprogcount;
        }
        /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
        /* Create top 10 list combining all nodes */
        MPI_Barrier (MPI_COMM_WORLD);
        for ( int i_top = 0; i_top < 10; ++i_top ) {
            int * mpi_cons_lib;
            mpi_cons_lib = new int;
            double * mpi_cons_dis;
            mpi_cons_dis = new double;
            char * mpi_cons_seq;
            mpi_cons_seq = new char[rawdat.length];
            std::string mpi_cons_seq_string;
            
            MPI_Barrier (MPI_COMM_WORLD);
            if ( mynode != 0 ) {
                *mpi_cons_lib = best_cons.at(i_top).libsize;
                *mpi_cons_dis = best_cons.at(i_top).distance;
                for ( int i_pos = 0; i_pos < rawdat.length; ++i_pos ) {
                    mpi_cons_seq[i_pos] = best_cons.at(i_top).seq[i_pos];
                }
                rc = MPI_Send ( mpi_cons_lib, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                rc = MPI_Send ( mpi_cons_dis, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
                for ( int i_pos = 0; i_pos < rawdat.length; ++i_pos ) {
                    rc = MPI_Send ( &mpi_cons_seq[i_pos], 1, MPI_CHAR, 0, i_pos + 3, MPI_COMM_WORLD);
                }
            }
            MPI_Barrier (MPI_COMM_WORLD);
            if ( mynode == 0 ) {
                for ( int i_node = 1; i_node < numnodes; ++i_node ) {
                    rc = MPI_Recv ( mpi_cons_lib, 1, MPI_INT, i_node, 1, MPI_COMM_WORLD, &Stat);
                    rc = MPI_Recv ( mpi_cons_dis, 1, MPI_DOUBLE, i_node, 2, MPI_COMM_WORLD, &Stat);
                    for ( int i_pos = 0; i_pos < rawdat.length; ++i_pos ) {
                        rc = MPI_Recv ( &mpi_cons_seq[i_pos], 1, MPI_CHAR, i_node, i_pos + 3, MPI_COMM_WORLD, &Stat);
                        mpi_cons_seq_string.push_back( mpi_cons_seq[i_pos] );
                    }
                    cons_props * temp_cons;
                    temp_cons = new cons_props;
                    temp_cons->set( *mpi_cons_dis, mpi_cons_seq_string, *mpi_cons_lib );
                    if ( best_cons.back().distance > temp_cons->distance ) sort_best ( &best_cons, temp_cons );
                    mpi_cons_seq_string.clear();
                    delete temp_cons;
                }
            }
            MPI_Barrier (MPI_COMM_WORLD);
            delete mpi_cons_lib;
            delete mpi_cons_dis;
            delete mpi_cons_seq;
        }
        MPI_Barrier (MPI_COMM_WORLD);
        /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
        
        
        /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
        /* Write top 10 list of this library size to file */
        if ( mynode == 0 ) {
            std::stringstream resultpath;
            resultpath << rawdat.output << rawdat.name << "_" << rawdat.libsizes[i_size].libsize << ".txt";
            std::ofstream resultfile;
            resultfile.open ( resultpath.str().c_str() );
            for ( int i = 0; i < 10; ++i ) {
                resultfile << best_cons.at(i).libsize << "\t";
                resultfile << best_cons.at(i).seq << "\t";
                resultfile << best_cons.at(i).distance;
                if ( i < 9 ) resultfile << "\n";
            }
            resultfile.close();
            
            std::stringstream log;
            log << " Library size " << rawdat.libsizes[i_size].libsize << " completed.\n";
            writelog ( rawdat.output, log.str() );
            log.str("");
        }
        /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
        
        
        /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
        /* Create top 10 list combining all library sizes */
        if ( mynode == 0 ) {
            for ( int i_best = 0; i_best < 10; ++i_best ) {
                if ( best_cons_tot.back().distance > best_cons.at(i_best).distance ) sort_best ( &best_cons_tot, &best_cons[i_best] );
            }
        }
        /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
        
        conspath.str("");
        delete [] consline;
    }
    
    /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    /* Write top 10 list of all library sizes to file */
    if ( mynode == 0 ) {
        std::stringstream resultpath;
        resultpath << rawdat.output << rawdat.name << "_total.txt";
        std::ofstream resultfile;
        resultfile.open ( resultpath.str().c_str() );
        for ( int i = 0; i < 10; ++i ) {
            resultfile << best_cons_tot.at(i).libsize << "\t";
            resultfile << best_cons_tot.at(i).seq << "\t";
            resultfile << best_cons_tot.at(i).distance;
            if ( i < 9 ) resultfile << "\n";
        }
        resultfile.close();
        
        std::stringstream log;
        log << " All library sizes completed.\n";
        writelog ( rawdat.output, log.str() );
        log.str("");
    }
    /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* Function to split calculations among nodes */
unsigned int combrange ( int mynode, int numnodes, unsigned int combs ) {
    int nodecount = 0;
    unsigned int index = 0;
    unsigned int rest = combs % numnodes;
    unsigned int part = ( combs - rest ) / numnodes;
    
    while ( nodecount < mynode ) {
        index = index + part;
        ++nodecount;
    }
    if ( mynode == numnodes ) index = index + rest;
    
    return index;
}
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* Function definitions for cons_props class */
cons_props::cons_props () {
    distance = DBL_MAX;
    seq = "X";
    for (int i = 1; i < 8; i++) seq = seq + "X";
    libsize = 0;
}
cons_props::cons_props ( double d, std::string s, int l ) {
    distance = d;
    seq = s;
    libsize = l;
}
cons_props::cons_props ( const cons_props& c ){
    distance = c.distance;
    seq = c.seq;
    libsize = c.libsize;
}
void cons_props::set ( double d, std::string s, int l ) {
    distance = d;
    seq = s;
    libsize = l;
}
void cons_props::set ( const cons_props& c ) {
    distance = c.distance;
    seq = c.seq;
    libsize = c.libsize;
}
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* Generate distribution from consensus and raw data, and evaluate Kolmogorov-Smirnov distance */
void get_cons_props ( cons_props * temp_cons, RAWDAT rawdat ) {
    std::stringstream err;
    std::string translation;                //containing translation of degenerate DNA code to string of specific DNA bases
    std::vector<double> distribution;
    std::vector<int> indices;
    std::vector<int> indices_temp;
    int check = 0;
    for ( unsigned int i = 0; i < rawdat.data.size (); ++i ) indices.push_back (i);
    
    // Generating a reduced sequence-value pair list from input data corresponding to a degenerate sequence:
    for ( int j = 0; j < rawdat.length; ++j ) {
        for ( unsigned int i = 0; i < indices.size (); ++i ) {
            // Translation of degenerate DNA code to string of specific DNA bases:
            switch ( temp_cons->seq[j] ) {
                case 'A':
                    translation = "A";
                    break;
                case 'C':
                    translation = "C";
                    break;
                case 'G':
                    translation = "G";
                    break;
                case 'T':
                    translation = "T";
                    break;
                case 'Y':
                    translation = "CT";
                    break;
                case 'R':
                    translation = "AG";
                    break;
                case 'W':
                    translation = "AT";
                    break;
                case 'S':
                    translation = "CG";
                    break;
                case 'K':
                    translation = "GT";
                    break;
                case 'M':
                    translation = "AC";
                    break;
                case 'B':
                    translation = "CGT";
                    break;
                case 'D':
                    translation = "AGT";
                    break;
                case 'H':
                    translation = "ACT";
                    break;
                case 'V':
                    translation = "ACG";
                    break;
                case 'N':
                    translation = "ACGT";
                break;
            }
            for ( unsigned int k = 0; k < translation.length (); ++k ) {
                if ( rawdat.data[indices[i]].seq[j] == translation[k] ) ++check;
            }
            if ( check > 0 ) indices_temp.push_back (indices[i]);
            check = 0;
        }
        indices = indices_temp;
        indices_temp.clear();
    }
    for ( unsigned int i = 0; i < indices.size (); ++i ) distribution.push_back( rawdat.data[indices[i]].level );

// Calculation of Kolmogorov-Smirnov distances:
    // For uniform distributions:
    if ( rawdat.distr == "uniform" ) {
        double step = 1.0 / double ( temp_cons->libsize );
        double skip = 1.0;
        double F = 0.0;     //F: ecdf of the evaluated sequence
        double Dsup = 0.0;  //Dsup: Kolmogorov-Smirnov distance
        double D1 = 0.0;
        double D2 = 0.0;
        double ucdf = 0.0;
        double range = rawdat.maxlev - rawdat.minlev;
        for ( int m = temp_cons->libsize - 1; m >= 0; --m ) {
            F = F + step;
            skip = 1.0;
            while ( distribution[m] == distribution[m - 1] && m > 0 ) {
                F = F + step;
                --m;
                skip = skip + 1.0;
            }
            if ( distribution.at(m) < rawdat.minlev ) {
                ucdf = 0.0;
            } else if ( distribution.at(m) > rawdat.maxlev ) {
                ucdf = 1.0;
            } else {
                ucdf = ( distribution.at(m) - rawdat.minlev ) / range;
            }
            D1 = std::abs ( F - ucdf );
            D2 = std::abs ( F - skip * step - ucdf );
            if ( D1 > Dsup ) Dsup = D1;
            if ( D2 > Dsup ) Dsup = D2;
        }
        temp_cons->set ( Dsup, temp_cons->seq, temp_cons->libsize );        
    // For normal distributions:
    } else if ( rawdat.distr == "normal" ) {
        double step = 1.0 / double ( temp_cons->libsize );
        double skip = 1.0;
        double F = 0.0;     //F: ecdf of the evaluated sequence
        double Dsup = 0.0;  //Dsup: Kolmogorov-Smirnov distance
        double D1 = 0.0;
        double D2 = 0.0;
        double ncdf = 0.0;
        for ( int m = temp_cons->libsize - 1; m >= 0; --m ) {
            F = F + step;
            skip = 1.0;
            while ( distribution[m] == distribution[m - 1] && m > 0 ) {
                F = F + step;
                --m;
                skip = skip + 1.0;
            }
            ncdf = 0.5 * ( 1 + erf ( ( distribution.at(m) - rawdat.mu ) / ( rawdat.sigma * sqrt (2) )  ) );
            D1 = std::abs ( F - ncdf );
            D2 = std::abs ( F - skip * step - ncdf );
            if ( D1 > Dsup ) Dsup = D1;
            if ( D2 > Dsup ) Dsup = D2;
        }
        temp_cons->set ( Dsup, temp_cons->seq, temp_cons->libsize );
    } else {
        err << "Error: only uniform and normal distribution implemented in this version.\n\n";
        std::cout << err.str();
        writelog ( rawdat.output, err.str() );
        err.str("");
        exit ( EXIT_FAILURE );
    }
}
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
/* Sort top 10 consensus sequences ascending according to KS distance */
void sort_best ( std::vector<cons_props> * best_cons, cons_props * temp_cons ){
    int pos = 9;
    int index = pos;
    cons_props mem2 ( *temp_cons );
    cons_props mem1;
    while ( temp_cons->distance < best_cons->at(index).distance && pos >= 0 ) {
        --pos;
        if ( pos >= 0 ) index = pos;
    }
    while ( pos < 9 ) {
        ++pos;
        mem1.set ( best_cons->at(pos) );
        best_cons->at(pos).set ( mem2 );
        mem2.set ( mem1 );
    }
}
/* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */





