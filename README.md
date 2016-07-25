<img src="https://www.bsse.ethz.ch/bpl/software/redlibs/_jcr_content/par/textimage/image.imageformat.lightbox.1971066330.png" alt="RedLibs logo" >
# RedLibs
RedLibs is an algorithm that allows for the rational design of smart combinatorial libraries for pathway optimization thereby minimizing the use of experimental resources.

This software is free software, licensed under the GNU GPL license v3.0.
Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2015. All rights reserved.
See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016

The source code provided here is written in C++ using Open MPI v1.4.5 for parallel computing. The files main.cpp, logfile.cpp, rawdata.cpp, and generate_list.cpp constitute together with their header files (.h) the program RedLibs. To generate a executable version, compile them together, e.g. "mpicxx -o ../RedLibsMPI.o -Wall main.cpp logfile.cpp rawdata.cpp generate_list.cpp -lstdc++".

The file "conslib.cpp" in folder "conslib" provides the code to generate the sequence database needed to run RedLibs. After generating the database, the filepath needs to be specified in the files "generate_list.cpp" after line 62 and "rawdata.cpp" after line 345.

The file "RedLibs_graphs.R" provides a R-script to graphically evaluate the output of RedLibs.
