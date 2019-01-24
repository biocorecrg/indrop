#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
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
output (output folder)        : ${params.output}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

genomeFile          = file(params.genome)
annotationFile      = file(params.annotation) 
configFile          = file(params.config) 
barcodeFile         = file(params.barcode_list) 

outputfolder    = "${params.output}"
outputQC        = "${outputfolder}/QC"
outputMultiQC   = "${outputfolder}/multiQC"
outputMapping   = "${outputfolder}/Alignments";
filt_folder     = "${outputfolder}/Tagged_reads";
est_folder      = "${outputfolder}/Estimated_counts";
rep_folder      = "${outputfolder}/Reports";

if( !barcodeFile.exists() ) exit 1, "Missing barcode file: ${barcodeFile}"
if( !genomeFile.exists() ) exit 1, "Missing genome file: ${genomeFile}"
if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"

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
    .set { read_pairs }

Channel
    .fromPath( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { reads_for_fastqc; fastq_files_for_size_est}    


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
    set pair_id, file("${pair_id}_tagged.fastq.gz") into tagged_files_for_alignment
    file("${pair_id}_tagged.fastq.gz") into tagged_files_for_fastqc
  
    """
        droptag -S -p ${task.cpus} -c ${configFile} ${reads}
        zcat *.tagged.*.gz >> ${pair_id}_tagged.fastq
        gzip ${pair_id}_tagged.fastq
        rm  *.fastq.gz.tagged.*.gz
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
    tag { fastq_file_for_size_est }

    input: 
    file(fastq_file_for_size_est) from fastq_files_for_size_est.first()
 
    output: 
    stdout into read_length

    script:
    def qc = new QualityChecker(input:fastq_file_for_size_est)
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
    set pair_id, file(tags) from STARmappedTags_for_est
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
    dropest -m -w -g ${annotationFile} -c ${configFile} -o ${pair_id}.est ${tags} 
    """
    

}

/*
*
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

*/

/*
*
*/
process dropReport {
    label 'indrop'
    publishDir rep_folder
    tag { pair_id }
    errorStrategy 'ignore'

    input:
    set pair_id, file(estimate) from estimates_rds

    output:
    set pair_id, file ("${pair_id}_report.html")  into outreport

    script:     
    """
    dropReport.Rsc -o ${pair_id}_report.html ${estimate}
    """
}



/*
 * Step 6. QualiMap QC. The default is using strand-specific-reverse. Should we try both directions? // better multiQC // we should try...

process qualimap {
    label 'big_mem_cpus'
    publishDir stat_folder, mode: 'copy'

    input:
    file annotationFile

    set pair_id, file(bamfolder) from STARmappedFolders_for_qualimap

    output:
    set pair_id, file("QUALIMAP_${pair_id}") into QualiMap

    script:
    //
    // Qualimap
    //
    """
    unset DISPLAY
    qualimap rnaseq --java-mem-size=16G -bam "${bamfolder}/${pair_id}Aligned.sortedByCoord.out.bam" -gtf ${annotationFile} -outdir "QUALIMAP_${pair_id}" -p ${qualiOption}
    """
}


 * Step 7. Multi QC.

    process multiQC_unfiltered {
        publishDir outputMultiQC

        input:
        file multiconfig
        file ribo_report
        file '*' from raw_fastqc_files.collect()
        file '*' from trimmed_fastqc_files.collect()        
        file '*' from STARmappedFolders_for_multiQC.collect()
        file '*' from logTrimming_for_QC.collect()
        file '*' from QualiMap.collect()

        output:
        file("multiqc_report.html") into multiQC 
    
        script:
         //
         // multiqc
         //
        """
        check_tool_version.pl -l fastqc,star,skewer,qualimap,ribopicker,bedtools,samtools > tools_mqc.txt
        multiqc -c ${multiconfig} .
        """
}

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




