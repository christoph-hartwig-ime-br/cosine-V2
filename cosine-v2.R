# Version 2.0
# Written 2019 by Dipl.-Ing.(FH) Christoph Hartwig (Fraunhofer IME-BR)
# Based on first implementation of the workflow by Dr. Florian Zubeil (Fraunhofer IME-BR) in 2017.
#
# Version history:
# V 1.0 2017 Florian Zubeil (Fraunhofer IME-BR) - https://dx.doi.org/10.5281/zenodo.3911715 - https://github.com/fzubeil-IME-BR/metabolomics_cosine
#   FeatureFinding and Bucketing via xcmsSet and xSet (import of bucket table from ProfileAnalysis/MetaboScape also implemented), cosine calculation with own code in C, heatmap as widget with d3heatmap and saveWidget
# V 2.0 2019 Christoph Hartwig (Fraunhofer IME-BR) - moving all calculations to R-packets - including export of results

# Data Processing is handled in DataAnalysis (Bruker) for feature finding [Reprocessed files are stored and can be used for future bucketings, so that reprocessig is necessary only once]
# Bucketing is done via ProfileAnalysis (Bruker) and that output file "XYZ_HPlus.txt" is the input for this script [ProfileAnalysis allows huge sample numbers, only limited by memory.]
# The following steps are performed in this script:
# Calculation of cosine similarities and generation of a heatmap (using cosine dis-similarity). Export of heatmap as .png .
# Export of the cosine similarity matrix [allowing further processing in e.g. Excel]
# Calculation and export of statistics about buckets and sample numbers, uniqueness etc. [feedback information to evaluate feature finding and bucketing parameters]
# Readout and export of clustering sequence and pairwise similarities [to find "jumps" and therefore finding borders of metabolomics clades]

library(readr)
library(coop) 
library(gplots) 
require(data.table) 
library(parallelDist) 


# Define functions for graphics generation
my_palette <- colorRampPalette(c("white","blue"), bias=10)(n = 500) # color scheme for heatmap

calcall<-function(project){
  
  heatpng<-function(mt, graphname){
    png(filename=paste(c(project,graphname), collapse=""),    # create PNG for the heat map        
        width = 50*400,       
        height = 50*400,
        res = 1000,           
        pointsize = 5)        # small font size accomodating for large sample numbers
    test2<-heatmap.2(mt, symm=TRUE, distfun=function(x) as.dist((1-x)/2), notecol="black", # as.dist((1-x)/2) inverts similarity to dissimilarity
                     main = paste(c(project,graphname), collapse=""), # heat map title
                     col=my_palette,       
                     margins =c(25,25),     # widens margins around plot allowing to read sample IDs
                     trace="none")
    dev.off()
    write.table(mt[rev(test2$rowInd),0], file=paste(c(project,"_cosclust_result.txt"), collapse=""), dec=".",sep="\t",row.names=TRUE,col.names=TRUE) # saving the order of similarity matrix in an .txt file
    clusstable2<<-(mt[rev(test2$rowInd),0])
    clusttable<<-rev(test2$rowInd)
  }
  
  
  table.gesamt<-read_tsv(paste(c(project,".txt"), collapse=""),col_names=TRUE, col_types = NULL, na = c("", "NA"), trim_ws = FALSE, skip = 0, n_max = Inf, progress = show_progress(), skip_empty_rows = TRUE)  # load project file
  
  
  table.gesamt<-as.data.frame(table.gesamt) # converts imported project file to data frame
  names<-table.gesamt[,1] # reading the filenames aka sample names
  n<-ncol(table.gesamt) # table length, necessary for matrix generation
  table.buckets<-table.gesamt[,2:n] # cutting out the filenames
  
  # readout of all bucketnames for search and export
  bucketnames<<-colnames(table.gesamt) # readout of bucketnames
  matrix.bucketnames<<-as.matrix(bucketnames) # transforming bucket name table to matrix
  rm(table.gesamt) # frees up memory, helpful for large datasets
  
  table.buckets<<-as.matrix(table.buckets) # transform buckettable to matrix
  n2<-nrow(table.buckets) # calculates the number of samples
  rownames(table.buckets)<-names # adding the filenames as rownames
  bucketnames2<-matrix.bucketnames[-1,] # extracts the bucketnames
  colnames(table.buckets)<-bucketnames2 # adding the bucketnames as colnames
  table.rotated<-t(table.buckets) # transformation of the buckettable for cosine calculation
  n.buckets.filled<<-nrow(counter<-as.matrix(table.buckets[!table.buckets==0])) # write filled buckets in vector and count lines
  ausgabe<-summary(table.buckets[!table.buckets==0]) # saving results in variable for output in summary
  rm(table.buckets) # frees up memory, helpful for large datasets
  
  # Analysis of filled buckets per sample
  table.buckets.logical<-table.rotated
  table.buckets.logical[table.buckets.logical > 1] <- 1
  n.buckets.sample<-as.matrix(colSums(table.buckets.logical)) # how many buckets are filled for every sample
  n.samples.bucket<-as.matrix(rowSums(table.buckets.logical)) # how many samples contain a certain feature, for every feature
  
  # for every sample the number of buckets which are only present in this sample
  n.uniquebuckets.sample<-as.matrix(table.buckets.logical)
  n.uniquebuckets.sample[rowSums(n.uniquebuckets.sample) > 1] <-0
  n.uniquebuckets.sample.result<-as.matrix(colSums(n.uniquebuckets.sample))
  n.uniquebuckets.bucket.result<-as.matrix(rowSums(n.uniquebuckets.sample))
  
  # Calculation of similarities
  table.cosine<<-cosine(table.rotated) # calculation of cosinus similarities based on bucket table
  rm(table.rotated)
  
  # Export of tables
  write.table(table.cosine, file=paste(c(project,"_cossim_result.txt"), collapse=""), dec=".",sep="\t",row.names=TRUE,col.names=NA) # saving the cosine similarity matrix in an .txt file
  write.table(n.buckets.sample, file=paste(c(project,"_numfilledbuckets.txt"), collapse=""), dec=".",sep="\t",row.names=TRUE,col.names=NA) # saving the number of filled buckets per sample in an .txt file
  write.table(matrix.bucketnames, file=paste(c(project,"_bucketnames.txt"), collapse=""), dec=".",sep="\t",row.names=TRUE,col.names=NA) # saving the bucket names in an .txt file
  write.table(n.samples.bucket, file=paste(c(project,"_numsamplesperbucket.txt"), collapse=""), dec=".",sep="\t",row.names=TRUE,col.names=NA) # saving the number of sample containing each bucket in an .txt file
  
  heatpng(table.cosine,"_cossim_result.png") # generates the heatmap output
  
  clusttable3<<-as.matrix(clusttable) # grabs the clustering sequence from the heatmap and uses that to read the pairwise similarities in that sequence. Output allows the generation of metabolic groups (based on clustering)
  groupingoutput <<- matrix(0 , ncol=2 ,nrow=nrow(clusttable3))
  for(i in 1:nrow(clusttable3)){
    groupingoutput[i,1]<<-(colnames(table.cosine, do.NULL)[clusttable3[i]])
  }
  for(i in 2:nrow(clusttable3)){
    groupingoutput[i,2]<<-(table.cosine[clusttable3[i-1],clusttable3[i]])
    
  }
  write.table(groupingoutput, file=paste(c(project,"_grouping_result.txt"), collapse=""), dec=".",sep="\t",row.names=FALSE,col.names=FALSE) # saving the cosine similarity matrix in an .txt file
  
  # generates a summary of the different outputs
  capture.output(
    print(project),
    print("summary(table.buckets[!table.buckets==0])"),
    ausgabe,
    print(""),
    print("Number of Buckets"),
    n,
    print(""),
    print("Number of Samples"),
    n2,
    print(""),
    print("total number of filled buckets"),
    n.buckets.filled,
    print(""),
    print("Number of Samples per bucket"),
    summary(n.samples.bucket),
    print(""),
    print("Number of Buckets per sample"),
    summary(n.buckets.sample),
    print(""),
    print("Number of unique buckets pers sample"),
    summary(n.uniquebuckets.sample.result),
    print(""),
    print("Number of unique buckets pers bucket"),
    summary(n.uniquebuckets.bucket.result),
    
    file=paste(c(project,"_summary.txt"), collapse=""), append=FALSE, type="output", split=FALSE )
  }

setwd("E:/") # Working directory
experiment<-"Dummy-Data"  #Experiment name, has to be Filename_HPlus.txt without _HPlus.txt as input, i.e. "Dummy-Data" for Dummy-Data_HPlus.txt
project<-(paste(c(experiment,"_HPlus"), collapse="")) #adds the _HPlus to adress the correct filename.

calcall(paste(c(experiment,"_HPlus"), collapse="")) #starts the complete script
