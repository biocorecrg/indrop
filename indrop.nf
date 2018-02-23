#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */
 

log.info "BIOCORE@CRG RNASeq - N F  ~  version 0.1"
log.info "========================================"
log.info "pairs                 	 : ${params.pairs}"
log.info "barcode        		 : ${params.barcode}"
log.info "chunkSize	            	 : ${params.chunkSize}"
log.info "mode (can be identity or 2mod) : ${params.mode}"
log.info "junction  			 : ${params.junction}"
log.info "UMI (lower limit)              : ${params.UMI}"
log.info "genome                         : ${params.genome}"
log.info "annotation                     : ${params.annotation}"
log.info "umisize                        : ${params.umisize}"
log.info "\n"

if (params.help) exit 1


genomeFile             = file(params.genome)
annotationFile         = file(params.annotation) 
barcodeFile        	= file(params.barcode) 

outputMapping	       = "Alignments";
stat_folder		= "Statistics";
filt_folder		= "Filtered_reads";

if( !genomeFile.exists() ) exit 1, "Missing genome file: ${genomeFile}"
if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"
if( !barcodeFile.exists() ) exit 1, "Missing annotation file: ${barcodeFile}"

if (params.strand == "yes") qualiOption = "strand-specific-forward"
else if (params.strand != "no") qualiOption = "non-strand-specific"
else if (params.strand != "reverse") qualiOption = "strand-specific-reverse"
else exit 1, "Please choose a valid strand option"

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .set { read_pairs }
    
/*
 * Step 1. It merges the two pairs and format them for being processed
 */
process mergePairs {    
    input:
    set pair_id, file(reads) from read_pairs
  
    output:
    set pair_id, "merged.txt" into merged_files
  
    """
    join_pairs.sh ${reads} > merged.txt
    """
}   

/*
 * Step 2. It filters chunks of merged paires
 */
process filterMerged {    
    input:
    set pair_id, file(mergeds) from merged_files.splitText(by:params.chunkSize, file:true)
 
    output:
    set pair_id, file("filtered.fq")  into filtered_files
    set pair_id, file("stats.txt")  into stats_files
    
    """
    indrops.pl -action preparereads -input ${mergeds} -junction ${params.junction} -size 8,${params.umisize}  -barcode_file ${params.barcode} -output filtered.fq -compare ${params.mode} > stats.txt
    """
} 

/*
 * Step 3. It cats the parsed fastq pieces in single files
 */
process catFiltered {    
	publishDir filt_folder, mode: 'copy'
    input:
    set pair_id, file('filtered_*') from filtered_files.groupTuple()
 
    output:
    set pair_id, file("${pair_id}_filtered.fq") into whole_filtered_files_for_stats
    set pair_id, file("${pair_id}_filtered.fq") into whole_filtered_files_for_size
    set pair_id, file("${pair_id}_filtered.fq") into whole_filtered_files_for_mapping
    
    """
    cat filtered_* > ${pair_id}_filtered.fq
    """
} 

/*
 * Step 4. It run statistics on whole parsed fastq files. WRONG???
 *
process makeStats {    
	publishDir stat_folder, mode: 'copy'
    input:
    set pair_id, file(whole_filtered) from whole_filtered_files_for_stats
 
    output:
    set pair_id, file("${pair_id}_cell_stats.txt") 
    set pair_id, file("${pair_id}_cell_stats.txt.dist") 
    
    """
    indrops.pl -action parsedstats -input ${whole_filtered} -umilimit ${params.UMI} -output ${pair_id}_cell_stats.txt
    """
} 
*/

/*
 * Step 5 
 */
process sumStats {    
	publishDir stat_folder, mode: 'copy'
    input:
    set pair_id, file('stats_*') from stats_files.groupTuple()
 
    output:
    set pair_id, file("${pair_id}_stats.txt") 

    """
    cat stats_* > temp_stats.txt
    indrops.pl -action joinstats -input temp_stats.txt -output ${pair_id}_stats.txt
    """
} 

/*
  * Step 6 extract read length
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
 * Step 7. Builds the genome index required by the mapping process
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

process parseMapping {
	publishDir outputMapping

	input:
	set pair_id, file(bamfolder) from STARmappedFolders_for_parsing
    	file annotationFile

	output:
	set pair_id, file("${pair_id}_parsed.aln")
	script:	
	
        """
	indrops.pl -action parsemapping -htseq_strand ${params.strand} -bam "${bamfolder}/${pair_id}Aligned.sortedByCoord.out.bam" -gtf ${annotationFile} -output "${pair_id}_parsed.aln" -umilimit ${params.UMI}
        """
}



/*
 * Step 7. QualiMap QC. The default is using strand-specific-reverse. Should we try both directions? // better multiQC // we should try...
*/
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

/*
 * Step 6. MutiQC QC.
*/
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
