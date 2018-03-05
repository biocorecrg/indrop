#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */
 

log.info "BIOCORE@CRG RNASeq - N F  ~  version 0.1"
log.info "========================================"
log.info "pairs               		 : ${params.pairs}"
log.info "genome               		 : ${params.genome}"
log.info "annotation           	     : ${params.annotation}"
log.info "config             	     : ${params.config}"
log.info "barcode_list               : ${params.barcode_list}"
log.info "output (output folder) 	 : ${params.output}"
log.info "email for notification 	 : ${params.email}"
log.info "\n"

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"


genomeFile          = file(params.genome)
annotationFile      = file(params.annotation) 
configFile        	= file(params.config) 
barcodeFile        	= file(params.barcode_list) 

outputfolder    = "${params.output}"
outputQC		= "${outputfolder}/QC"
outputMapping   = "${outputfolder}/Alignments";
filt_folder		= "${outputfolder}/Tagged_reads";
est_folder		= "${outputfolder}/Estimated_counts";
plots_folder	= "${outputfolder}/Plots";


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
    .into { read_pairs; reads_for_fastqc; fastq_files_for_size_est }


/*
 * Step 0. Run FastQC on raw data
*/
process QConRawReads {
	publishDir outputQC

	tag { read }

    input:
    set pair_id, file(read) from reads_for_fastqc

     output:
   	 file("*_fastqc.*") into raw_fastqc_files

     script:

    """
		fastqc ${read} 
    """
}

    
/*
 * Step 1. Launch droptag for tagging your files
 */
process dropTag {    
	publishDir filt_folder

	tag { pair_id }

    input:
    set pair_id, file(reads) from read_pairs
    file configFile
    
    output:
    set pair_id, "${pair_id}_tagged.fastq.gz" into tagged_files_for_alignment
    set pair_id, "${pair_id}_tagged.fastq.gz" into tagged_files_for_fastqc
  
    """
		droptag -S -p ${task.cpus} -c ${configFile} ${reads}
        zcat *.tagged.*.gz >> ${pair_id}_tagged.fastq
        gzip ${pair_id}_tagged.fastq
        rm 	*.fastq.gz.tagged.*.gz
    """
}   


/*
 * Step 2. FastQC of your trimmed files
 */

process QCFiltReads {
	publishDir outputQC

	tag { pair_id }

   	 input:
     set pair_id, file(filtered_read) from tagged_files_for_fastqc

     output:
   	 file("*_fastqc.*") into trimmed_fastqc_files

     script:

    """
		fastqc ${filtered_read}
    """
   }

/*
 * Step 3 extract read length of filtered reads?
*/

process getReadLength {   
	tag { fastq_files_for_size_est[1] }

    input: 
     set pair_id, file(fastq_files_for_size_est) from fastq_files_for_size_est.first()

 
	output: 
	stdout into read_length


    """
    estim_read_size.sh ${fastq_files_for_size_est[1]}
    """
} 


 /*
 * Step 4. Builds the genome index required by the mapping process
 */

    
process buildIndex {
	tag { genomeFile }

    input:
    file genomeFile
    file annotationFile
    val read_size from read_length.map {  it.trim().toInteger() } 

    output:
    file "STARgenome" into STARgenomeIndex

    echo true
    
    """
    echo read size is ${read_size}
	index_star.sh ${genomeFile} STARgenome STARgenome ${annotationFile} ${read_size-1} ${task.cpus}
    """
}

 /*
 * Step 5. Mapping with STAR
 */

process mapping {
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
	mappingPairs( pair_id, STARgenome, reads, task.cpus)  
}


process dropEst {
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

process dropReport {
	publishDir est_folder
	tag { pair_id }

	input:
	set pair_id, file(estimate) from estimates_rds

	output:
	set pair_id, file ("report.html")  into outreport

	script:		
    """
    dropReport.Rsc ${estimate}
    """
}

*/

/*
 * Step 6. QualiMap QC. The default is using strand-specific-reverse. Should we try both directions? // better multiQC // we should try...

process qualimap {
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
	    publishDir outputQC

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

/*
 * Mail notification

 
workflow.onComplete {
    def subject = 'indrop execution'
    def recipient = "${params.email}"
    def attachment = "${outputQC}/multiqc_report.html"

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
 */

/*
* ************** CUSTOM FUNCTIONS ****************
*/

def mappingPairs( pair_id, STARgenome, reads, cpus) { 
    """
		STAR --genomeDir ${STARgenome} \
				 --readFilesIn ${reads} \
				 --outSAMunmapped None \
				 --outSAMtype BAM SortedByCoordinate \
				 --runThreadN ${cpus} \
				 --quantMode GeneCounts \
				 --outFileNamePrefix ${pair_id} \
				 --readFilesCommand zcat			 	

			mkdir STAR_${pair_id}
			mv ${pair_id}Aligned* STAR_${pair_id}/.
			mv ${pair_id}SJ* STAR_${pair_id}/.
			mv ${pair_id}ReadsPerGene* STAR_${pair_id}/.
			mv ${pair_id}Log* STAR_${pair_id}/.   
    """
}
