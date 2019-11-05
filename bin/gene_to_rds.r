#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
infile	<- args[1]
outfile	<- args[2]

if (length(args)<2) {stop("[USAGE] R --slave --args <infile> <outfile> < gene_to_rds")}

a<-read.table(infile)
b<-as.vector(a$V1)
saveRDS(b, outfile)
