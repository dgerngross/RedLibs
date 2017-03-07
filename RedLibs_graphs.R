##
##  RedLibs_graphs.R
##  RedLibsMPI
##
##  Updated on 02/03/2017
##
##  This software is free software, licensed under the GNU GPL license v3.0.
##  Copyright (c) ETH Zurich, D-BSSE, BPL, Daniel Gerngross 2017. All rights reserved.
##  See http://www.bsse.ethz.ch/bpl/software/redlibs for updates, documentation, questions, and answers.
##  Please cite M. Jeschek, D. Gerngross, S. Panke, Nature Communications, 2016 (DOI:10.1038/ncomms11163)
##


#############
## Inputs: ##
#############
#install.packages(xlsx) #Install the package 'xlsx' if not already done
outputPath          <- "/Volumes/user/RedLibs/results/myProject"            #Path to output the graphs and spreadsheet
dataPath            <- "/Volumes/user/RedLibs/results/myProject/data.csv"   #Path to degenerate input library file
RedLibsOutputPath   <- "/Volumes/user/RedLibs/results/myProject/output.txt" #Path to RedLibs output file
name                <- "gene"                           #Name of project
level               <- "Translation Initiation Rate"    #Name data type assigned to sequences
distributionMin     <- TRUE #Lower margin of target uniform distribution, TRUE if absolut minimum of input data
distributionMax     <- TRUE #Upper margin of target uniform distribution, TRUE if absolut maximum of input data
#############

########################################
## Ininitialization and Reading Data: ##
########################################
library(xlsx)
setwd(outputPath)

data <- read.csv( dataPath, header=FALSE )

dataSequence  <- as.character(data[,1])               #Read list of sequences
m         <- length(dataSequence)                     #Length of data set
l         <- nchar(dataSequence[1])                   #Length of degenerate sequence
dataSequence <- data.frame(matrix(
  unlist(strsplit(dataSequence, split='')),
  nrow=m, byrow=TRUE))                            #Convert character strings to character vectors as matrix rows 
dataLevel <- data[, 2]                             #Read "level" values

dataCombined <- data.frame(dataSequence, level=dataLevel)  #Combine needed data into data.frame

dataMax   <- data[1,2]                           #Get data value maximum
dataMin   <- data[m,2]                           #Get data value minimum
dataRange <- dataMax - dataMin

if ( distributionMax == TRUE ) distributionMax <- dataMax        #If no custum maximum provided, target distribution maximum equals data maximum
if ( distributionMin == TRUE ) distributionMin <- dataMin        #If no custum minimum provided, target distribution minimum equals data minimum
dis_range  <- distributionMax - distributionMin

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

uniform_cdf <- data.frame ( c(dataMin - 0.05*dataRange, distributionMin, distributionMax, dataMax + 0.05*dataRange), c(0,0,1,1) )

RedLibsOutput <- read.table ( RedLibsOutputPath, header=FALSE )
degeneracy <- RedLibsOutput[1,1]

pdfFile <- paste ( outputPath, "/", name, "_", degeneracy, ".pdf", collapse="", sep="" )
xlsxFile <- paste( outputPath, "/", name, "_", degeneracy, ".xlsx", collapse="", sep="" )
pdf ( file=pdfFile, family="Helvetica", pointsize=10, width=8.3, height=5.8)
for ( i in 1:dim(RedLibsOutput)[1] ) {
  if ( RedLibsOutput[i,1] == degeneracy ) {
  consensus <- unlist ( strsplit ( as.character(RedLibsOutput[i,2]), split='' ) )
  distr  <- dataCombined
  for ( j in 1:l ) {
    distr  <- distr[distr[, j] %in% acgt[translate(consensus[j])], ]
  }
  xlsxDistribution <- data.frame( Sequence=rep(NA,degeneracy), Level=rep(NA,degeneracy))
  for ( j in 1:degeneracy ) {
    xlsxDistribution$Sequence[j] <- paste( as.character( unlist( distr[j,1:l] ) ), collapse="", sep="" )
    xlsxDistribution$Level[j] <- distr$level[j]
  }
  xlsxParameters <- data.frame( c( "Sequence", "Degeneracy",
                              "dKS", "Distribution",
                              "Target Minimum", "Target Maximum" ),
                           c( as.character(RedLibsOutput[i,2]), as.character(degeneracy),
                              as.character(RedLibsOutput[i,3]), "uniform",
                              as.character(distributionMin), as.character(distributionMax) ) )
  if ( i==1 ) {
    append <- FALSE
  } else append <- TRUE
  write.xlsx( xlsxDistribution, xlsxFile, sheetName=as.character(RedLibsOutput[i,2]), row.names=FALSE, col.names=TRUE, append=append )
  write.xlsx( xlsxParameters, xlsxFile, sheetName=paste( as.character(RedLibsOutput[i,2]), "_parameters", collapse="", sep="" ), row.names=FALSE, col.names=FALSE, append=TRUE )
  
  if ( max(distr$level) > distributionMax ) {
    breaksMax <- max(distr$level)
  } else breaksMax <- distributionMax
  if ( degeneracy <= 20 ) {
    breaks <- seq ( dataMin, breaksMax, (breaksMax - dataMin)/19 )
  } else if (degeneracy > 200 ) {
    breaks <- seq ( dataMin, breaksMax, (breaksMax - dataMin)/199 )
  } else breaks <- seq ( dataMin, breaksMax, (breaksMax - dataMin)/(degeneracy-1) )
  
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
                     2, 2, 4), nrow=11, ncol=3, byrow=TRUE) )
  
  par( mar=c(0,5,5,1) )
  boxplot( distr[,l+1], seq( distributionMin, distributionMax, dis_range/(degeneracy-1) ), horizontal=TRUE, boxwex=0.5, main=name,
            ylim=c(distributionMin, breaksMax), axes=FALSE, col=c("#a8322d", "#1f407a"), lwd=1.5, bty="1")
  box( bty="7" )
  axis( 2, at=c(0,3), tck=0, col.axis="white" )
  abline( v=distributionMin, col="grey", lty=2 )
  abline( v=distributionMax, col="grey", lty=2 )
  
  par( mar=c(5,5,0,1) )
  hist( distr[,l+1], breaks=breaks, main="", xlab=level, col="#a8322d", border="#6D0712", yaxs="i", xlim=c(dataMin, breaksMax) )
  abline( h=0, col="grey" )
  abline( v=distributionMin, col="grey", lty=2 )
  abline( v=distributionMax, col="grey", lty=2 )
  box( bty="u" )
  
  par( mar=c(5,4.2,1,1) )
  xval1 <- c( dataMin - 0.05*dataRange, dataMin )
  xval2 <- knots( ecdf(distr[,l+1]) )
  xval3 <- c( dataMax, dataMax + 0.05*dataRange )
  xval <- c( xval1, xval2, xval3 )
  
  ecdf_y <- ecdf(distr[,l+1])(xval)
  unif_y <- NA
  for( k in 1:length(xval) ) {
    if( xval[k] <= distributionMin ) {
      unif_y[k] <- 0
    } else if ( xval[k] >= distributionMax ) {
      unif_y[k] <- 1
    } else unif_y[k] <- (xval[k] - distributionMin) / dis_range
  }
  
  diff1 <- abs( ecdf(distr[,l+1])(xval) - unif_y )
  diff2 <- abs( ecdf(distr[,l+1])(xval)[1:length(xval)-1] - unif_y[2:length(xval)] )
  diff <- c( diff1, diff2 )
  
  xval_diff <- c( xval, xval[2:length(xval)] )
  diff_max <- max(diff)
  dsup_index <- which( diff == diff_max )
  dsup_x <- xval_diff[dsup_index]
  
  if( dsup_x <= distributionMin ) {
    dsup_y1 <- 0
  } else if( dsup_x >= distributionMax ) {
    dsup_y1 <- 1
  } else dsup_y1 <- ( dsup_x - distributionMin ) / dis_range

  if( unif_y[which(xval == dsup_x)] + diff_max == ecdf_y[which(xval == dsup_x)] ){
    dsup_y2 <- dsup_y1 + diff_max
  } else dsup_y2 <- dsup_y1 - diff_max
  
  plot(ecdf(distr[,l+1]), main="", pch="", verticals=TRUE, xlab=level, ylab="cdf",
      col="#a8322d", lwd=1.5, xlim=c(distributionMin-0.05*(breaksMax - distributionMin), breaksMax+0.05*(breaksMax - distributionMin)), xaxs="i" )
  lines( uniform_cdf, col="#1f407a", lwd=1.5 )
  arrows( dsup_x, dsup_y1, dsup_x, dsup_y2, length=0.05, lty=3, col=rgb(0,0,0,0.8), code=3 )
  if( dsup_x > dataRange/2 ) pos <- 2 else pos <- 4
  text( dsup_x, min(c(dsup_y1, dsup_y2))+0.5*abs(dsup_y1-dsup_y2),
         bquote(d[KS] == .(round(max(diff), digits=4))), pos=pos, col=rgb(0,0,0,0.8) )
    
  boxplot( c(diff1[3:(length(diff1)-2)], diff2[2:(length(diff2)-2)]), ylab="d" )
  abline( h=max(diff), lty=2 )
  text( 0.6, diff_max, expression(d[KS]), pos=1 )

  par( mar=c(1,2.5,2.5,1), family="mono"  )
  text_result <- paste("Sequence     = ", as.character(RedLibsOutput[i,2]), "\nDegeneracy   = ", RedLibsOutput[i,1],
                       "\ndKS          = ", RedLibsOutput[i,3], "\nDistribution = uniform",
                       "\nMinimum      = ", round(distributionMin), "\nMaximum      = ", round(distributionMax), collapse="", sep="" )
  plot( c(0,100), c(0,100), pch=NA, axes=FALSE, xlab="", ylab="")
  points( 5, 10, pch=15, col="#1f407a", cex=2 )
  points( 5, 1, pch=15, col="#a8322d", cex=2 )
  text( 7, 10, "uniform", pos=4 )
  text( 7, 1, "rationally reduced library", pos=4 )
  legend( "topleft", text_result, box.col="white", bg=rgb(1,1,1,0.5) )
  par( family="" )
  }
}
invisible( dev.off() )
