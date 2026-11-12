#!/usr/bin/env nextflow

/*
 * The following pipeline uses the markers positions of your current genome assembly
 * and finds the positions of these markers in a new genome version the user provides.
 */


// Defaults; override in params.config or on the command line (--multi_ratio 0.98)
params.paralog_ratio = 0.90   // competitor >= this * best bitscore -> Unique_paralog
params.multi_ratio   = 0.98   // competitor >= this * best bitscore -> Multi-mapped
params.merge_dist    = 1000   // HSPs closer than this (bp) on one sequence = one locus
params.blast_dust    = 'yes'  // 'no' disables low-complexity masking of queries


process Index_genome {
    tag "Indexing genome ${ref_genome.name}"
    publishDir "${params.reference}", mode: 'copy', overwrite: true

    input:
    path ref_genome

    output:
    path "${ref_genome}.fai", emit: fai

    script:
    """
    samtools faidx ${ref_genome}
    """
}


process Obtain_fasta {
    tag "Obtain FASTA sequences of markers ${ref_file.name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path old_chip
    path ref_file
    path ref_index

    output:
    path "seq_${params.v_current}.fasta", emit: seqs
    path "markers_na.txt",                emit: markers_na
    path "flank_start.txt",               emit: flank_start

    script:
    """
    python ${params.scripts}/flank_seq.py ${ref_file} ${old_chip} seq_${params.v_current}.fasta markers_na.txt flank_start.txt ${params.flank_size}
    """
}


process MakeBlastDB {
    tag "Make BLAST DB from ${ref_fasta.name}"

    input:
    path ref_fasta

    output:
    path "${ref_fasta.baseName}.*", emit: blast_db_files

    script:
    """
    makeblastdb -in ${ref_fasta} -parse_seqids -dbtype nucl -out ${ref_fasta.baseName}
    """
}


process BlastSearch {
    tag "BLAST sequences ${sequences_fasta.name}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path sequences_fasta
    path blast_db_files

    output:
    path "blast_${params.v_new}_results.tsv", emit: results

    script:
    """
    blastn -query ${sequences_fasta} -db ${blast_db_files[0].baseName} -out blast_${params.v_new}_results.tsv -outfmt '6 std qseq sseq' -dust ${params.blast_dust}
    """
}


process ExtractHitIds {
    tag "Extract hit IDs from BLAST results"

    input:
    path blast_results_file

    output:
    path "hit_fasta_${params.v_new}_ids.txt", emit: ids

    script:
    """
    # Extract query IDs (column 1), sort uniquely
    cut -f1 ${blast_results_file} | sort | uniq > hit_fasta_${params.v_new}_ids.txt
    """
}


process FindFailedBlastIds {
    tag "Find failed BLAST IDs"

    input:
    path full_fasta
    path hit_ids_file

    output:
    path "all_fasta_ids.txt"
    path "all_fasta_ids_sorted.txt"
    path "hit_fasta_ids_sorted.txt"
    path "failed_blast_${params.v_new}.txt", emit: failed

    script:
    """
    grep '^>' ${full_fasta} | sed 's/>//' > all_fasta_ids.txt
    sort all_fasta_ids.txt > all_fasta_ids_sorted.txt
    sort ${hit_ids_file} > hit_fasta_ids_sorted.txt

    # Output 4 columns to match the pipeline format
    comm -23 all_fasta_ids_sorted.txt hit_fasta_ids_sorted.txt | \
    awk '{print \$1"\tn.a.\tn.a.\tFailed_BLAST"}' > failed_blast_${params.v_new}.txt
    """
}


process BlastFilter {
    tag "Filter BLAST results and select position"
    //publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path blast_results
    path markers_file

    output:
    path "filtered_${params.v_new}.csv",  emit: filtered
    path "filtered_failed.tsv",           emit: failed
    path "positions_${params.v_new}.csv", emit: positions

    script:
    """
    python ${params.scripts}/blast_res_filter.py ${blast_results} filtered_${params.v_new}.csv filtered_failed.tsv ${params.flank_size}
    python ${params.scripts}/select_position.py filtered_${params.v_new}.csv positions_${params.v_new}.csv ${params.paralog_ratio} ${params.multi_ratio} ${params.merge_dist}
    """
}


process CreateChip {
    tag "Create final CHIP with positions in genome assembly version ${params.v_new}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path positions_file
    path length_file
    path markers_na
    path failed_blast
    path filtered_failed

    output:
    path "hits_${params.v_new}.txt"
    path "chip_${params.v_new}.txt"

    script:
    """
    python ${params.scripts}/find_position.py ${positions_file} ${length_file} hits_${params.v_new}.txt
    python ${params.scripts}/concatinate.py hits_${params.v_new}.txt ${markers_na} ${failed_blast} ${filtered_failed} chip_${params.v_new}.txt
    """
}


workflow {

    log.info """\
             ${workflow.manifest.name} v${workflow.manifest.version}
             Author : ${workflow.manifest.author}
             """.stripIndent()

    old_chip = file( params.markers )
    new_ref  = file( params.ref_new )
    old_ref  = file( params.ref_current )

    Index_genome( old_ref )

    Obtain_fasta( old_chip, old_ref, Index_genome.out.fai )

    MakeBlastDB( new_ref )

    BlastSearch( Obtain_fasta.out.seqs, MakeBlastDB.out.blast_db_files )

    ExtractHitIds( BlastSearch.out.results )

    FindFailedBlastIds( Obtain_fasta.out.seqs, ExtractHitIds.out.ids )

    BlastFilter( BlastSearch.out.results, old_chip )

    CreateChip(
        BlastFilter.out.positions,
        Obtain_fasta.out.flank_start,
        Obtain_fasta.out.markers_na,
        FindFailedBlastIds.out.failed,
        BlastFilter.out.failed
    )
}
