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

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ╦┌┐┌┌┬┐┬─┐┌─┐┌─┐╔═╗╦  ╔═╗╦ ╦
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ║│││ ││├┬┘│ │├─┘╠╣ ║  ║ ║║║║
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ╩┘└┘─┴┘┴└─└─┘┴  ╚  ╩═╝╚═╝╚╩╝
                                                                                       
====================================================
BIOCORE@CRG indropSEQ - N F  ~  version ${version}
====================================================
pairs                         : ${params.pairs}
genome                        : ${params.genome}
annotation                    : ${params.annotation}
config                        : ${params.config}
barcode_list                  : ${params.barcode_list}
email                         : ${params.email}
mtgenes                       : ${params.mtgenes}
dbdir                         : ${params.dbdir}
version                       : ${params.version}
keepmulti                     : ${params.keepmulti}
library_tag                   : ${params.library_tag}
output (output folder)        : ${params.output}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

genomeFile          = file(params.genome)
annotationFile      = file(params.annotation) 
configFile          = file(params.config) 
barcodeFile         = file(params.barcode_list) 
if (params.mtgenes != "") mitocgenesFile  = file(params.mtgenes) 
//mitocgenesFile      = file(params.mtgenes)
db_folder		    = file(params.dbdir)
dropestScript       = file("$baseDir/docker/dropestr/dropReport.Rsc")

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
//if( !mitocgenesFile.exists() ) exit 1, "Missing mitocondrial genes file: ${mitocgenesFile}"

/*
* if (params.strand == "yes") qualiOption = "strand-specific-forward"
* else if (params.strand != "no") qualiOption = "non-strand-specific"
*/

if (params.version != "1-2" && params.version != "3_3" && params.version != "3_4") 
	exit 1, "Please define a valid version! It can be 1-2, 3_3, 3_4.\nRespectively version 1 or 2, version 3 with 3 files and version 4 with 4 files."

if (params.keepmulti != "NO" && params.keepmulti != "YES") 
	exit 1, "Please define a valid keepmulti value! It can YES or NO"

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.pairs, size: (params.version == "1-2") ? 2 : (params.version == "3-3") ? 3 : 4)                                          
    .ifEmpty { error "Cannot find any reads matching: ${params.pairs}" }  
    .into { read_pairs; fastq_files_for_size_est; ciccio}

Channel
    .fromPath( params.pairs )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { reads_for_fastqc}    


/*
 * Run FastQC on raw data
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

reads = Channel.empty()
fastq_file_for_size_est = Channel.empty()

if (params.version == "3_4") {
	read_pairs.map{
		[it[0], [it[1][1], it[1][3], it[1][0], it[1][2] ]]
	}.set{ reads}
	
	fastq_files_for_size_est.map{
		[it[0], [it[1][0]]]
	}.set{ fastq_file_for_size_est}
}    
else if (params.version == "3_3") {
	read_pairs.map{
		[it[0], [it[1][1], it[1][2], it[1][0] ]]
	}.set{ reads}
	fastq_files_for_size_est.map{
		[it[0], [it[1][0]]]
	}.set{ fastq_file_for_size_est}	
} else {
	read_pairs.set{reads}
	fastq_files_for_size_est.first().map{
		[it[0], [it[1][1]]]
	}.set{ fastq_file_for_size_est}
}



/*
 * Launch droptag for tagging your files
 */
process dropTag {    
    publishDir filt_folder
    label 'indrop_multi_cpus'

    tag { pair_id }

    input:
    set pair_id, file(seqs) from reads
    file configFile
    
    output:
    set pair_id, file("*.tagged.fastq.gz") into tagged_files_for_alignment
    file("*.tagged.fastq.gz") into tagged_files_for_fastqc
    set pair_id, file("*.tagged.params.gz") into params_files_for_estimation
    set pair_id, file("*.tagged.rds") into tagged_rds_for_report
  
    //zcat *.tagged.*.gz >> ${pair_id}_tagged.fastq
    //gzip ${pair_id}_tagged.fastq 
    script:
    def v3_params = ""
    if (params.library_tag != "") {
    	v3_params = "-t ${params.library_tag}"
    }
    """
        droptag -r 0 -S -s ${v3_params} -p ${task.cpus} -c ${configFile} ${seqs}
    """
    
}   


/*
 * FastQC of your trimmed files
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
 * Extract read length of filtered reads?
*/

process getReadLength {   
    tag { pairid }

    input: 
    
    set pairid, file(fastq_file) from fastq_file_for_size_est
 
    output: 
    stdout into read_length

    script:
    def qc = new QualityChecker(input:fastq_file)
    qc.getReadSize()
} 


 /*
 * Builds the genome index required by the mapping process
 */

    
process buildIndex {
    tag { genomeFile }
   // storeDir db_folder
    
    label 'index_mem_cpus'

    input:
    file genomeFile
    file annotationFile
    val read_size from read_length.map {  it.trim().toInteger() } 

    output:
    file ("STARgenome") into STARgenomeIndex
    
    script:
    def aligner = new NGSaligner(reference_file:genomeFile, index:"STARgenome", annotation_file:annotationFile, read_size:read_size-1, cpus:task.cpus)
    aligner.doIndexing("STAR")
}

/*
 * Mapping with STAR
 */

process mapping {
    label 'big_mem_cpus'
    tag { pair_id }

    input:
    file STARgenome from STARgenomeIndex
    set pair_id, file(reads) from tagged_files_for_alignment    

    output:
    set pair_id, file("STAR_${pair_id}/${pair_id}*.bam") into STARmappedTags_for_filter
    set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_qualimap
    set pair_id, file("STAR_${pair_id}") into STARmappedFolders_for_multiQC

// To add --outFilterMultimapNmax 1?
    script:
    """
    STAR   --runThreadN ${task.cpus} \
    --genomeDir ${STARgenome}   --readFilesIn ${reads}  --readFilesCommand zcat \
    --outFileNamePrefix ${pair_id}  --outSAMtype BAM Unsorted --outSAMunmapped None
        mkdir STAR_${pair_id} 
        mv ${pair_id}Aligned* STAR_${pair_id}/. 
        mv ${pair_id}SJ* STAR_${pair_id}/. 
        mv ${pair_id}Log* STAR_${pair_id}/. 
        if test -f "${pair_id}ReadsPerGene*"; 
           then mv ${pair_id}ReadsPerGene* STAR_${pair_id}/.
        fi       
    """
    //def aligner = new NGSaligner(id:pair_id, reads:reads, index:STARgenome, cpus:task.cpus,  output:"STAR_${pair_id}") 
    //aligner.doAlignment("STAR")  

}

if (params.keepmulti == "NO") {
	/*
	 * Removing multi-mapping
	 */
	 
	process removeMultimapping {
		publishDir outputMapping
		label 'big_mem_cpus'
		tag { pair_id }

		input:
		set pair_id, file(aln) from STARmappedTags_for_filter    

		output:
		set pair_id, file ("${pair_id}_univoc_s.bam") into STARmappedTags_for_est
	
		script:
		"""
		samtools view -H ${aln} > ${pair_id}_univoc_s.sam
		samtools view -@ ${task.cpus} ${aln} | grep \"\\<NH:i:1\\>\" >> ${pair_id}_univoc_s.sam
		samtools view -@ ${task.cpus} -Sb ${pair_id}_univoc_s.sam > ${pair_id}_univoc_s.bam 
		rm ${pair_id}_univoc_s.sam
		"""
 
	}
} else {
   STARmappedTags_for_est = STARmappedTags_for_filter
}


process dropEst {
    label 'indrop_one_cpu'
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
process dropReport {
    label 'dropreport'
    errorStrategy = 'ignore'
    publishDir rep_folder, mode: 'copy'
    tag { pair_id }

    input:
    set pair_id, file(estimate), file (droptag) from estimates_rds.join(tagged_rds_for_report)
    file (dropestScript)
    
    output:
    set pair_id, file ("${pair_id}_report.html")  into outreport

    script:
    def mitopar = ""
    def mitocmd = ""
    if (params.mtgenes != "") {
        mitopar = " -m mitoc.rds" 
        mitocmd = "gene_to_rds.r ${mitocgenesFile} mitoc.rds"
    }
    """
    ${mitocmd}
    Rscript --vanilla ${dropestScript} -t ${droptag} -o ${pair_id}_report.html ${mitopar} ${estimate} 
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


workflow.onComplete {
    println "Pipeline BIOCORE@indrop-Flow completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
* send mail
*/
workflow.onComplete {
    def subject = 'indrop-Flow execution'
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




