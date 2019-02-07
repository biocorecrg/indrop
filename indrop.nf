#!/usr/bin/env nextflow

/* 
 * Copyright (c) 2019, Centre for Genomic Regulation (CRG) 
 * 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 */
 
/* 
 * Define the pipeline parameters
 * 
 * Author Luca Cozzuto
 */

// Pipeline version
version = '1.0'

params.help            = false
params.resume          = false

log.info """

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ┬┌┐┌┌┬┐┬─┐┌─┐┌─┐╔═╗╔═╗╔═╗ 
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ││││ ││├┬┘│ │├─┘╚═╗║╣ ║═╬╗
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ┴┘└┘─┴┘┴└─└─┘┴  ╚═╝╚═╝╚═╝╚
                                                                                       
====================================================
BIOCORE@CRG indropSEQ - N F  ~  version ${version}
====================================================
pairs                         : ${params.pairs}
genome                        : ${params.genome}
annotation                    : ${params.annotation}
config                        : ${params.config}
barcode_list                  : ${params.barcode_list}
email                         : ${params.email}
mtgenes                       : ${params.mitocgenes}
output (output folder)        : ${params.output}
storeIndex
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

genomeFile          = file(params.genome)
annotationFile      = file(params.annotation) 
configFile          = file(params.config) 
barcodeFile         = file(params.barcode_list) 
mitocgenesFile      = file(params.mitocgenes)

outputfolder    = "${params.output}"
outputQC        = "${outputfolder}/QC"
outputMultiQC   = "${outputfolder}/multiQC"
outputMapping   = "${outputfolder}/Alignments"
filt_folder     = "${outputfolder}/Tagged_reads"
est_folder      = "${outputfolder}/Estimated_counts"
rep_folder      = "${outputfolder}/Reports"

if( !barcodeFile.exists() ) exit 1, "Missing barcode file: ${barcodeFile}"
if( !genomeFile.exists() ) exit 1, "Missing genome file: ${genomeFile}"
if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"
if( !mitocgenesFile.exists() ) exit 1, "Missing mitocondrial genes file: ${mitocgenesFile}"

/*
* if (params.strand == "yes") qualiOption = "strand-specific-forward"
* else if (params.strand != "no") qualiOption = "non-strand-specific"
*/

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .into { read_pairs; fastq_files_for_size_est}

Channel
    .fromPath( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { reads_for_fastqc}    


/*
 * Step 0. Run FastQC on raw data
*/
process QConRawReads {
    publishDir outputQC

    tag { read }

    input:
    file(read) from reads_for_fastqc

    output:
    file("*_fastqc.*") into raw_fastqc_files

    script:
    def qc = new QualityChecker(input:read, cpus:task.cpus)
    qc.fastqc()
}

    
/*
 * Step 1. Launch droptag for tagging your files
 */
process dropTag {    
    publishDir filt_folder
    label 'indrop'

    tag { pair_id }

    input:
    set pair_id, file(reads) from read_pairs
    file configFile
    
    output:
    set pair_id, file("*.tagged.1.fastq.gz") into tagged_files_for_alignment
    file("*.tagged.1.fastq.gz") into tagged_files_for_fastqc
    set pair_id, file("*.tagged.params.gz") into params_files_for_estimation
    set pair_id, file("*.tagged.rds") into tagged_rds_for_report
  
    //zcat *.tagged.*.gz >> ${pair_id}_tagged.fastq
    //gzip ${pair_id}_tagged.fastq 
    """
        droptag -S -s -p ${task.cpus} -c ${configFile} ${reads}
    """
}   


/*
 * Step 2. FastQC of your trimmed files
 */

process QCFiltReads {
    publishDir outputQC

    tag { filtered_read }

    input:
    file(filtered_read) from tagged_files_for_fastqc.flatten()

    output:
    file("*_fastqc.*") into trimmed_fastqc_files

    script:
    def qc = new QualityChecker(input:filtered_read, cpus:task.cpus)
    qc.fastqc()
   }

/*
 * Step 3 extract read length of filtered reads?
*/

process getReadLength {   
    tag { pairid }

    input: 
    
    set pairid, file(fastq_files) from fastq_files_for_size_est.first()
 
    output: 
    stdout into read_length

    script:
    def right_pair = fastq_files.last()
    def qc = new QualityChecker(input:right_pair)
    qc.getReadSize()
} 


 /*
 * Step 4. Builds the genome index required by the mapping process
 */

    
process buildIndex {
    tag { genomeFile }
    label 'index_mem_cpus'


    input:
    file genomeFile
    file annotationFile
    val read_size from read_length.map {  it.trim().toInteger() } 

    output:
    file "STARgenome" into STARgenomeIndex
    
    script:
    def aligner = new NGSaligner(reference_file:genomeFile, index:"STARgenome", annotation_file:annotationFile, read_size:read_size-1, cpus:task.cpus)
    aligner.doIndexing("STAR")
}

/*
 * Step 5. Mapping with STAR
 */

process mapping {
    label 'big_mem_cpus'
    publishDir outputMapping
    tag { pair_id }

    input:
    file STARgenome from STARgenomeIndex
    set pair_id, file(reads) from tagged_files_for_alignment    

    output:
    set pair_id, file("STAR_${pair_id}/${pair_id}Aligned.sortedByCoord.out.bam") into STARmappedTags_for_est
    set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_qualimap
    set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_multiQC

    script:
    def aligner = new NGSaligner(id:pair_id, reads:reads, index:STARgenome, cpus:task.cpus, output:"STAR_${pair_id}") 
    aligner.doAlignment("STAR")  

}


process dropEst {
    label 'indrop'
    publishDir est_folder
    tag { pair_id }

    input:
    set pair_id, file(tags), file(params_est) from STARmappedTags_for_est.join(params_files_for_estimation)
    file ("barcode_file.txt") from barcodeFile
    file annotationFile
    file configFile

    output:
    set pair_id, file ("*.rds")  into estimates_rds
    set pair_id, file ("*.tsv")  into estimates_tsv
    set pair_id, file ("*.mtx")  into estimates_mtx_for_plots
    set pair_id, file ("*.mtx")  into estimates_mtx

    script:     
    """
    dropest -r ${params_est} -m -w -g ${annotationFile} -c ${configFile} -o ${pair_id}.est ${tags} 
    """
    

}

/*
*
*/ 
process histCell {
    publishDir plots_folder
    tag { pair_id }

    input:
    set pair_id, file(estimates) from estimates_mtx_for_plots
    file ("barcode_file.txt") from barcodeFile
    file annotationFile
    file configFile

    output:
    set pair_id, file ("*.rds")  into estimates_rds

    script:     
    """
    tail -n +3  ${estimates} |awk '{if (\$3>0) print \$2}'|sort|uniq -c >  ${estimates}_stats.txt
    plot.histogram.r ${estimates}_stats.txt ${pair_id}.est.pdf ${pair_id}
    """


/*
*
*/
process dropReport {
    label 'dropReport'
    publishDir rep_folder, mode: 'copy'
    tag { pair_id }

    input:
    set pair_id, file(estimate), file (droptag) from estimates_rds.join(tagged_rds_for_report)
    file(mitocgenes) from mitocgenesFile
    
    output:
    set pair_id, file ("${pair_id}_report.html")  into outreport

    script:     
    """
    dropReport.Rsc -t ${droptag} -o ${pair_id}_report.html -m ${mitocgenes} ${estimate} 
    """
}



/*
 * Step 7. Multi QC.
*/
    process multiQC_unfiltered {
        publishDir outputMultiQC

        input:
        file '*' from raw_fastqc_files.collect()
        file '*' from STARmappedFolders_for_multiQC.collect()

        output:
        file("multiqc_report.html") into multiQC 
    
        script:
         //
         // multiqc
         // check_tool_version.pl -l fastqc,star,skewer,qualimap,ribopicker,bedtools,samtools > tools_mqc.txt
         //
        """
        multiqc .
        """
}

/*
* send mail
*/
workflow.onComplete {
    def subject = 'indropSeq execution'
    def recipient = "${params.email}"
    def attachment = "${outputMultiQC}/multiqc_report.html"

    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}




