# ![indrop](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) Indrop-Flow

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Nextflow version](https://img.shields.io/badge/nextflow-%E2%89%A50.31.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Build Status](https://img.shields.io/docker/automated/biocorecrg/indrops.svg)](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/indrops/builds)


Indrops analysis pipeline at BioCore@CRG

The pipeline is based on the DropEST tool:
https://github.com/hms-dbmi/dropEst

## We support versions v1, v2 and v3
### Version 1 and 2
**Parameter version: "V1-2"**
* File 1: barcode reads. Structure:
  * Cell barcode, part 1
  * Spacer
  * Cell barcode, part 2
  * UMI
* File 2: gene reads

### Version 3
**Parameter version: "V3_3"**
* File 1: cell barcode
* File 2: 
  * cell barcode 
  * UMI 
* File 3: gene read

**Parameter version: "V3_4"**
* File 1: cell barcode
* File 2: 
  * cell barcode 
  * UMI 
* File 3: gene read
* File 4: library_tag

The parameter **library_tag** is only needed with version **V3_4**

## Install
1. install docker or singularity.
1. git clone https://github.com/biocorecrg/indrop.git; cd indrop
1. sh INSTALL.sh for checking Nextflow and installing bioNextflow

# Parameters
1. Parameters are specified within the **params.config** file

## The pipeline
1. QC: Run FastQC on raw reads. It stores the results within **QC** folder.
1. Indexing: It makes the index of the genome by using STAR.
1. dropTag: It creates a "tagged" fastq file with information about the single cell that originated that read in the header. 
1. Alignment: It aligns tagged reads to the indexed genome by using STAR. Reasults are stored in **Alignments** folder.
1. dropEst: It provides the estimation of read counts per gene per single cell. The results are in **Estimated_counts** folder and consists of an R data object, a file with a list of cells (aka barcode combinations), another with a list of genes and a matrix in Matrix Market format (https://en.wikipedia.org/wiki/Matrix_Market_exchange_formats).
1. dropReport: It reads the R data oject produced by the dropEst step to produce a quality report. It needs a list of mitochondrial genes. 
1. multiQC: It wraps the QC from fastQC and STAR mapping in a single output. 

