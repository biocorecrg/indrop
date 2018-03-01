# ![BiocoreRNAseq](https://github.com/CRG-CNAG/BioCoreMiscOpen/blob/master/logo/biocore-logo_small.png) indrop

Official indrops analysis pipeline at BioCore@CRG

Test data are in **MISSING**

Docker image is in indrops_docker repository

Once the pipeline is finished you can receive a mail with attached the MultiQC report.


Here you have an example of a report:
**MISSING**

-----
To run the pipeline you have to clone this repository and the corresponding docker FILE for creating the docker / singularity image.
image. 
If you are using the CRG cluster you don't need to create the singularity image since it is already available.
The config file contains information about location of the singularity image and whether to use or not singularity and requirements (like memory, CPUs etc) for every step.

**Important!! Check if your nextflow is updated to the latest version!!! (type nextflow self-update)**

## Parameters
To check the required parameters you can type nextflow run 

*nextflow run indrop.nf --help*

|parameter name         | value|
|---------------------------------|------------------------|
|pairs         | /nfs/software/bi/biocore_tools/git/nextflow/indrop/indrop/test/SRR1784313_{1,2}.fastq.gz|
|genome|/db/ensembl/release-88/mus_musculus/genome/Mus_musculus.GRCm38.dna.chrom.fa|
|annotation|/db/ensembl/release-88/mus_musculus/gtf/Mus_musculus.GRCm38.88.chr.gtf|
|config|/nfs/software/bi/biocore_tools/git/nextflow/indrop/indrop/test/indrop_v1_2.xml|
|barcode_list|/nfs/software/bi/biocore_tools/git/nextflow/indrop/indrop/test/indrop_v1_2.txt|
|output (output folder)|output|
|email for notification|yourmail@crg.eu|


### Reads
**Important** when specifying the parameters **reads** you should use **" "** if not the * will be transalted in the first file. Be careful with the way you group together the files since filenames can vary among facilities, machines etc.

### Genome
The genome can be also gzipped (it will be decompressed on the fly and indexed using STAR).

### Annotation
The annotation must be in gtf format. 

### Config
The config file is and xml file required by the program **dropEst** for working. Have a look at the one oin the test folder.

### barcode_list
A file containing in two rows the list of barcodes used. For some unknown reason **dropEst** requires their sequence to be **reverse complement** 

### Output
Is a parameter that specify the output folder. It is useful in case you want to have different run changing some parameters. Default is **output**.

### Email
This parameter is useful to receive a mail once the process is finished / crashed. Default is a fake mail.



## The pipeline
1. QC: Run FastQC on raw reads.
1. dropTag: It uses dropTag tool to tag the sequences looking at barcodes, junction sequence and polyA tail. Joining the output fastq files in a single one.
1. Indexing: It makes the index using the estimated read length, the annotation file and the genome sequence. It stores the results within **Index** folder.
1. QC of filtered reads: it performs a round of fastQC runs on tagged reads.
1. Alignment: It aligns either tagged reads to the reference transcriptome (i.e. to the genome but considering also the splicing sites) and it stores the results within the **Alignments** folder. 
1. dropEst: It uses dropEst tool for estimating the number of reads per gene per cell. It produces a rds R file and a matrix in MatrixMarket format.
1. multiQC: mapped reads are analyzed by qualimap and the outcomes grouped together with the other QC results by multiQC. The results are stored in **QC** folder.
Change information from config.yaml for multiQC to fetch information about user, date, contact email, etc.
1. Intermediate files are stored within the **work** folder.
