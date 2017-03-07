![RedLibs logo](http://www.bpl.ethz.ch/software/redlibs/RedLibs200px.png)
# RedLibs v1.1.0
## General information:
RedLibs is an algorithm that allows for the rational design of smart libraries for pathway optimization thereby minimizing the use of experimental resources.

This software is free software, licensed under the GNU GPL license v3.0. Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved. See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.

*Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 ([DOI:10.1038/ncomms11163](http://www.nature.com/articles/ncomms11163))*

## Generating a degenerate sequence database:
In order to generate the “conslib” database, the *C++* file *conslib.cpp* or *conslib_GLOS.cpp* needs to be compiled.

The file *conslib.cpp8 in folder *conslib* provides the code to generate the sequence database needed to run RedLibs based on fully degenerate (N) input libraries. Example usage for a degenerate sequence length of four (NNNN): `conslib.o 4`

The file *conslib_GLOS.cpp* in folder *conslib_GLOS* provides the code to generate the sequence database needed to run RedLibs based on genome library optimized sequences (GLOS) libraries, i.e. a six bases long triple degenerate sequence (e.g. BDHVBD). Example usage: `conslib_GLOS.o BDHVBD`

(*Please cite Sabine Oesterle, Daniel Gerngross, Steven Schmitt, Tania Michelle Roberts, and Sven Panke, [submission in progress], 2017*)

The database is then generated in the working directory. The corresponding path, where the database is located, will be one of the necessary inputs for the RedLibs algorithm (see below).

## Running the main RedLibs algorithm:
The source code provided here is written in *C++* using *Open MPI v1.4.5* for parallel computing. The files main.cpp, logfile.cpp, rawdata.cpp, and generate_list.cpp constitute together with their header files (.h) the program RedLibs. To generate a executable version, compile them together, e.g. `mpic++ -o ../RedLibsMPI.o -std=c++11 generate_list.cpp logfile.cpp rawdata.cpp main.cpp`.

The user needs to provide a degenerate input library file (csv format). This file contains all explicit sequences based on the degenerate sequence and a corresponding numerical value for each sequence sorted ascending in respect to these values, e.g.:

```
TAGGGA,9166.6
TGGGGA,4227.0
CAGGGA,3898.1
[...],[...]
TGCTCA,2.9
CGCTCA,2.9
```

### Usage:
`RedLibs.o [OPTIONS]... [path of conslib database] [path of degenerate input library file]`

### Example:
`RedLibs.o -g 18 -s 1 -e 6 /Volumes/user/RedLibs/conslib /Volumes/user/RedLibs/results/myProject`

### Options
```
-o, --output		path for results
                    default: [working directory]/[date]_[project name]
-n, --projectName	name of project
                    default: name of degenerate library file
-d, --distribution	type of target distribution: 'uniform' or 'normal'
                    default: uniform
-s, --startPosition	start position of degenerate sequence in raw data
                    default value: 1
-e, --endPosition	end position of degenerate sequence in raw data has to be specified
-g, --degeneracy	total degeneracy
                    default value: 2
-l, --minLevel		minimal sequence data level
                    default value: last value in input library file (sorted ascending)
-L, --maxLevel		maximal sequence data level
                    default value: first value in input library file (sorted ascending)
-h, --help          show help
```

## Generating graphs and spreadsheets from the output:
The file *RedLibs_graphs.R* provides a *R*-script to graphically evaluate the output of RedLibs for uniform target distributions. As output the user receives a pdf containing a graphical representation of the reduced libraries and a xlsx file containing the corresponding data. At the beginning of the file (lines 17-24, see below) the user has to specify paths for the output, file containing the input degenerate library, and RedLibs output file. If the range of the target distribution is not between the maximum and minimum of the input degenerate library, this also has to be specified.
### Example input:
```
#############
## Inputs: ##
#############
#install.packages(xlsx) #Install the package 'xlsx' if not already done
outputPath          <- "/Volumes/user/RedLibs/results/myProject"            #Path to output the graphs and spreadsheet
dataPath            <- "/Volumes/user/RedLibs/results/myProject/data.csv"   #Path to degenerate input library file
RedLibsOutputPath   <- "/Volumes/user/RedLibs/results/myProject/output.txt" #Path to degenerate input library file
name                <- "gene"                           #Name of project
level               <- "Translation Initiation Rate"    #Name data type assigned to sequences
distributionMin     <- TRUE #Lower margin of target uniform distribution, TRUE if absolut minimum of input data
distributionMax     <- TRUE #Upper margin of target uniform distribution, TRUE if absolut maximum of input data
#############
```
### Example output:
![RedLibs example output](http://www.bpl.ethz.ch/software/redlibs/RedLibs_output_example.png)