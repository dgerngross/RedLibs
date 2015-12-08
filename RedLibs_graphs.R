//
//  RedLibs_graphs.R
//  RedLibsMPI
//
//  Updated on 07/10/2015
//
//  This software is free software, licensed under the GNU GPL license v3.0.
//  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2015. All rights reserved.
//  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
//  Please cite M. Jeschek, D. Gerngross, S. Panke, [submission in progress]
//


#############
## Inputs: ##
#############
name     <- "gene"          #Name of project, input file, and folder
full.deg <- "NNNNNNN AACTCGAG NTG"   #Fully degenerate sequence from prediction
level    <- "Translation Initiation Rate"            #Name data type assigned to sequences
length   <- 8               #Length of fully degenerate sequence
start    <- 1               #Start position of degenerate sequence
end      <- 8               #End position of degenerate sequence
size     <- 81              #Target degeneracy
dis_min  <- TRUE            #Lower margin of target uniform distribution, TRUE if absolut minimum of input data
dis_max  <- TRUE            #Upper margin of target uniform distribution, TRUE if absolut maximum of input data
#############

########################################
## Ininitialization and Reading Data: ##
########################################
library(xlsx)
library(grImport2)
library(grConvert)
setwd("/Volumes/gerdanie$/RedLibs")
ETH.col <- c( "#1f407a", "#3c5a0f", "#0069b4",    #ETH corporate design colors:
              "#72791c", "#91056a", "#6f6f6e",    #1:dark blue, 2:dark green, 3:light blue,
              "#a8322d", "#007a92", "#956013" )   #4:grass green, 5:violet, 6:grey,
                                                  #7:red, 8:turquoise, 9:khaki
logo <- readPicture("Logo/RedLibs_version-pdf.svg")

data <- read.table ( paste("Results/", name, "/", name, ".txt", collapse="", sep=""), header=FALSE)


data.seq  <- as.character(data[,1])               #Read list of sequences
m         <- length(data.seq)                     #Length of data set
l         <- length                               #Length of degenerated sequence
data.seq <- data.frame(matrix(
  unlist(strsplit(data.seq, split='')),
  nrow=m, byrow=TRUE))                            #Convert character strings to character vectors as matrix rows 
data.deg <- data.seq[, seq ( start, end, 1 ) ]    #Save only degenerated part of sequences
data.lev <- data[, 2]                             #Read "level" values

data.dat <- data.frame(data.deg, level=data.lev)  #Combine needed data into data.frame

data_max   <- data[1,2]                           #Get data value maximum
data_min   <- data[m,2]                           #Get data value minimum
data_range <- data_max - data_min

if ( dis_max == TRUE ) dis_max <- data_max        #If no custum maximum provided, target distribution maximum equals data maximum
if ( dis_min == TRUE ) dis_min <- data_min        #If no custum minimum provided, target distribution minimum equals data minimum
dis_range  <- dis_max - dis_min

acgt      <- c('A','C','G','T')                   #Vector for base to number translation

translate <-function(x){
  switch(x,
         'A'=c(TRUE,FALSE,FALSE,FALSE),
         'C'=c(FALSE,TRUE,FALSE,FALSE),
         'G'=c(FALSE,FALSE,TRUE,FALSE),
         'T'=c(FALSE,FALSE,FALSE,TRUE),
         'Y'=c(FALSE,TRUE,FALSE,TRUE),
         'R'=c(TRUE,FALSE,TRUE,FALSE),
         'W'=c(TRUE,FALSE,FALSE,TRUE),
         'S'=c(FALSE,TRUE,TRUE,FALSE),
         'K'=c(FALSE,FALSE,TRUE,TRUE),
         'M'=c(TRUE,TRUE,FALSE,FALSE),
         'B'=c(FALSE,TRUE,TRUE,TRUE),
         'D'=c(TRUE,FALSE,TRUE,TRUE),
         'H'=c(TRUE,TRUE,FALSE,TRUE),
         'V'=c(TRUE,TRUE,TRUE,FALSE),
         'N'=c(TRUE,TRUE,TRUE,TRUE)
  )
}

uniform_cdf <- data.frame ( c(data_min - 0.05*data_range, dis_min, dis_max, data_max + 0.05*data_range), c(0,0,1,1))

consensus_file <- paste ( "results/", name, "/", name, "_", size, ".txt", collapse = "", sep = "" )
try ( best_list <- read.table ( consensus_file, header=FALSE ) )

if(best_list[1,1]==size){

pdf_file <- paste ( "results/", name, "/", name, "_", size, ".pdf", collapse = "", sep = "")
xlsx_file <- paste( "results/", name, "/", name, "_", size, ".xlsx" , collapse="", sep="")
pdf ( file = pdf_file, family = "Helvetica", pointsize = 10, width=8.3, height=5.8)
for ( i in 1:10 ) {
  if ( best_list[i,1] == size ) {
  consensus <- unlist ( strsplit ( as.character(best_list[i,2]), split='' ) )
  reduced.deg <-  unlist ( strsplit ( full.deg, split='' ) )
  reduced.deg[which(reduced.deg == "N")] <- consensus
  reduced.deg <- paste( reduced.deg, collapse = "", sep = "" )
  distr  <- data.dat
  for ( j in 1:l ) {
    distr  <- distr[distr[, j] %in% acgt[translate(consensus[j])], ]
  }
  xlsx_distr <- data.frame( Sequence=rep(NA,size), Level=rep(NA,size))
  for ( j in 1:size ) {
    temp.sequence <- unlist ( strsplit ( full.deg, split='' ) )
    temp.sequence[which(temp.sequence == "N")] <- as.character( unlist( distr[j,1:l] ) )
    xlsx_distr$Sequence[j] <- paste( temp.sequence, collapse = "", sep = "" )
    xlsx_distr$Level[j] <- distr$level[j]
  }
  xlsx_pars <- data.frame( c( "Sequence", "Degeneracy",
                              "dKS", "Distribution",
                              "Target Minimum", "Target Maximum" ),
                           c( reduced.deg, as.character(size),
                              as.character(best_list[i,3]), "uniform",
                              as.character(dis_min), as.character(dis_max) ) )
  if ( i==1 ) {
    append <- FALSE
  } else append <- TRUE
  write.xlsx( xlsx_distr, xlsx_file, sheetName=reduced.deg, 
              row.names=FALSE, col.names=TRUE, append=append )
  write.xlsx( xlsx_pars, xlsx_file, sheetName=paste( reduced.deg, "_parameters", collapse = "", sep = "" ), 
              row.names=FALSE, col.names=FALSE, append=TRUE )
  
  if ( max(distr$level) > dis_max ) {
    breaks_max <- max(distr$level)
  } else breaks_max <- dis_max
  if ( size <= 20 ) {
    breaks <- seq ( data_min, breaks_max, (breaks_max - data_min)/19 )
  } else if (size > 200 ) {
    breaks <- seq ( data_min, breaks_max, (breaks_max - data_min)/199 )
  } else breaks <- seq ( data_min, breaks_max, (breaks_max - data_min)/(size-1) )
  
  par ( family="" )
  layout( matrix( c( 1, 1, 5,
                     1, 1, 5,
                     2, 2, 5,
                     2, 2, 3,
                     2, 2, 3,
                     2, 2, 3,
                     2, 2, 3,
                     2, 2, 4,
                     2, 2, 4,
                     2, 2, 4,
                     2, 2, 4), nrow = 11, ncol = 3, byrow=TRUE) )
  
  par( mar = c(0,5,5,1) )
  boxplot ( distr[,l+1], seq ( dis_min, dis_max, dis_range/(size-1) ), horizontal = TRUE, boxwex = 0.5, main = name,
            ylim = c(dis_min, breaks_max), axes = FALSE, col = c(ETH.col[7], ETH.col[1]), lwd = 1.5, bty = "1")
  box ( bty = "7" )
  axis ( 2, at = c(0,3), tck = 0, col.axis = "white" )
  abline ( v = dis_min, col = "grey", lty = 2 )
  abline( v = dis_max, col = "grey", lty = 2)
  
  par( mar = c(5,5,0,1) )
  hist ( distr[,l+1], breaks = breaks, main = "", xlab = level, col = ETH.col[7], border = "#6D0712",
         yaxs = "i", xlim = c(data_min, breaks_max) )
  abline ( h = 0, col = "grey" )
  abline ( v = dis_min, col = "grey", lty = 2 )
  abline ( v = dis_max, col = "grey", lty = 2 )
  box ( bty = "u" )
  
  par ( mar = c(5,4.2,1,1) )
  xval1 <- c( data_min - 0.05*data_range, data_min )
  xval2 <- knots ( ecdf(distr[,l+1]) )
  xval3 <- c( data_max, data_max + 0.05*data_range )
  xval <- c( xval1, xval2, xval3 )
  
  ecdf_y <- ecdf(distr[,l+1])(xval)
  unif_y <- NA
  for ( k in 1:length(xval) ) {
    if ( xval[k] <= dis_min ) {
      unif_y[k] <- 0
    } else if ( xval[k] >= dis_max ) {
      unif_y[k] <- 1
    } else unif_y[k] <- (xval[k] - dis_min) / dis_range
  }
  
  diff1 <- abs ( ecdf(distr[,l+1])(xval) - unif_y )
  diff2 <- abs( ecdf(distr[,l+1])(xval)[1:length(xval)-1] - unif_y[2:length(xval)] )
  diff <- c( diff1, diff2 )
  
  xval_diff <- c( xval, xval[2:length(xval)] )
  diff_max <- max(diff)
  dsup_index <- which( diff == diff_max )
  dsup_x <- xval_diff[dsup_index]
  
  if ( dsup_x <= dis_min ) {
    dsup_y1 <- 0
  } else if ( dsup_x >= dis_max ) {
    dsup_y1 <- 1
  } else dsup_y1 <- ( dsup_x - dis_min ) / dis_range

  if( unif_y[which(xval == dsup_x)] + diff_max == ecdf_y[which(xval == dsup_x)] ){
    dsup_y2 <- dsup_y1 + diff_max
  } else dsup_y2 <- dsup_y1 - diff_max
  
  plot (ecdf(distr[,l+1]), main = "", pch = "", verticals = TRUE, xlab = level, ylab = "cdf",
      col = ETH.col[7], lwd = 1.5, xlim = c(dis_min-0.05*(breaks_max - dis_min), breaks_max+0.05*(breaks_max - dis_min)), xaxs = "i" )
  lines ( uniform_cdf, col = ETH.col[1], lwd = 1.5 )
  arrows ( dsup_x, dsup_y1, dsup_x, dsup_y2, length = 0.05, lty = 3, col = rgb(0,0,0,0.8), code = 3 )
  if ( dsup_x > data_range/2 ) pos <- 2
  else pos <- 4
  text ( dsup_x, min(c(dsup_y1, dsup_y2))+0.5*abs(dsup_y1-dsup_y2),
         bquote(d[KS] == .(round(max(diff), digits = 4))), pos = pos, col = rgb(0,0,0,0.8) )
    
  boxplot( c(diff1[3:(length(diff1)-2)], diff2[2:(length(diff2)-2)]), ylab = "d" )
  abline( h=max(diff), lty=2 )
  text( 0.6, diff_max, expression(d[KS]), pos=1 )

  par ( mar = c(1,2.5,2.5,1), family = "mono"  )
  text_result <- paste("Sequence     = ", reduced.deg, "\nDegeneracy   = ", best_list[i,1],
                       "\ndKS          = ", best_list[i,3], "\nDistribution = uniform",
                       "\nMinimum      = ", round(dis_min), "\nMaximum      = ", round(dis_max), collapse="", sep="" )
  plot ( c(0,100), c(0,100), pch=NA, axes=FALSE, xlab="", ylab="")
  points ( 5, 10, pch = 15, col = ETH.col[1], cex = 2 )
  points ( 5, 1, pch = 15, col = ETH.col[7], cex = 2 )
  text ( 7, 10, "uniform", pos = 4 )
  text ( 7, 1, "rationally reduced library", pos = 4 )
  legend ( "topleft", text_result, box.col = "white", bg = rgb(1,1,1,0.5) )
  grid.picture(logo, width=0.15, x=0.915, y=0.05, hjust="right", vjust="bottom")
  par ( family = "" )
  }
}
invisible(dev.off())
}