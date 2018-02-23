#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */
 

log.info "BIOCORE@CRG RNASeq - N F  ~  version 0.1"
log.info "========================================"
log.info "pairs                  : ${params.pairs}"
log.info "genome                 : ${params.genome}"
log.info "annotation             : ${params.annotation}"
log.info "config                 : ${params.config}"
log.info "\n"

if (params.help) exit 1


genomeFile          = file(params.genome)
annotationFile      = file(params.annotation) 
barcodeFile        	= file(params.barcode) 
configFile        	= file(params.config) 

outputMapping   = "Alignments";
stat_folder		= "Statistics";
filt_folder		= "Filtered_reads";

if( !genomeFile.exists() ) exit 1, "Missing genome file: ${genomeFile}"
if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"

if (params.strand == "yes") qualiOption = "strand-specific-forward"
else if (params.strand != "no") qualiOption = "non-strand-specific"

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .set { read_pairs }
    
/*
 * Step 1. Launch droptag for tagging your files
 */
process dropTag {    
    input:
    set pair_id, file(reads) from read_pairs
    file configFile
    
    output:
    set pair_id, "*.tagged.*.fastq.gz" into tagged_files
  
    """
		droptag -S -p ${task.cpus} -c ${configFile} ${reads} 
    """
}   



/*
 * Step 2 extract read length
 */
process getReadLength {   
    input: 
    set pair_id, file(single_whole_filtered) from whole_filtered_files_for_size.first()
 
	output: 
	stdout into read_length


    """
    estim_read_size.sh ${single_whole_filtered}
    """
} 


 /*
 * Step 3. Builds the genome index required by the mapping process
 */
 
    
process buildIndex {
    input:
    file genomeFile
    file annotationFile
     val read_size from read_length.map {  it.trim().toInteger() } 

    output:
    file "STARgenome" into STARgenomeIndex
    
    """
	index_star.sh ${genomeFile} STARgenome STARgenome ${annotationFile} ${read_size-1} ${task.cpus}
    """
}

/*
process mapping {
	publishDir outputMapping

	input:
	file STARgenome from STARgenomeIndex
	set pair_id, file(reads) from whole_filtered_files_for_mapping

	output:
	set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_parsing
	set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_qualimap
        set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_multiQC

	script:		
	mappingPairs( pair_id, STARgenome, reads, task.cpus)  
}






 * Step 7. QualiMap QC. The default is using strand-specific-reverse. Should we try both directions? // better multiQC // we should try...

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



process multiQC {
    publishDir stat_folder, mode: 'copy'

    input:
    file 'STAR_*' from STARmappedFolders_for_multiQC.collect()
    file 'QUALIMAP_*' from QualiMap.collect()

    output:
    file("multiqc_report.html") into multiQC

    script:
    //
    // multiqc
    //
    """
    multiqc .
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
				 --outFileNamePrefix ${pair_id}
							 	

			mkdir STAR_${pair_id}
			mv ${pair_id}Aligned* STAR_${pair_id}/.
			mv ${pair_id}SJ* STAR_${pair_id}/.
			mv ${pair_id}ReadsPerGene* STAR_${pair_id}/.
			mv ${pair_id}Log* STAR_${pair_id}/.   
    """
}
