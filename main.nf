#!/usr/bin/env nextflow

/*
 * The following pipeline uses the markers positions of your current genome assembly
 * and finds the positions of these markers in a new genome version the user provides.
 * 
 */


log.info """\
         ${params.manifest.name} v${params.manifest.version}
         Author : ${params.manifest.author}
         ===================================================
         Pipeline to extract marker positions
         ----
         ref current from   : ${params.ref_current}
         ref new from      : ${params.ref_new}
         output to    : ${params.outdir}

         """
         .stripIndent()



scripts = file (params.scripts)

old_chip = file( params.markers )
new_ref = file( params.ref_new )
old_ref = file( params.ref_current )



process Index_genome {
    tag "Indexing genome ${ref_genome.name}"
    
    publishDir = [
            path: "${params.reference}",
            mode: 'copy', overwrite: true
        ]


    input:
    path ref_genome

    output:
    path "${ref_genome}.fai"

    script:
    """
    samtools faidx ${ref_genome}
    """
}



process Obtain_fasta {
    tag{"Obtain FASTA sequences of markers ${old_ref.name}"}
    publishDir = [
            path: "${params.outdir}",
            mode: 'copy', overwrite: true
        ]

    input:
    path old_chip
    path ref_file
    path ref_index

    output:
    path "seq_${params.v_current}.fasta"
    path "markers_na.txt"
    path "flank_start.txt"

    script:
    """

    # obtain sequence
    python ${scripts}/flank_seq.py $ref_file $old_chip seq_${params.v_current}.fasta markers_na.txt flank_start.txt ${params.flank_size}

    
    """
}

process MakeBlastDB {
    tag { "Make BLAST DB from ${ref_fasta.name}" }
    //publishDir path: "${params.outdir}/db", mode: 'copy', overwrite: true

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
    tag { "BLAST sequences ${sequences_fasta} againts db ${ref_fasta}" }
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path sequences_fasta
    path blast_db_files

    output:
    path "blast_${params.v_new}_results.tsv"

    script:
    """
    blastn -query ${sequences_fasta} -db ${blast_db_files[0].baseName} -out blast_${params.v_new}_results.tsv -outfmt 6
    """
}

process ExtractHitIds {
    tag "Extract hit IDs from BLAST results"

    //publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path blast_results_file

    output:
    path "hit_fasta_${params.v_new}_ids.txt"

    script:
    """
    # Extract subject IDs (column 2), sort uniquely
    cut -f1 ${blast_results_file} | sort | uniq > hit_fasta_${params.v_new}_ids.txt
    """
}


process FindFailedBlastIds {
    tag "Find failed BLAST IDs"

    //publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path full_fasta     // e.g., seq_ref2.1.fasta or all_seqv2.fasta
    path hit_ids_file   // e.g., hit_fasta_ids.txt (list of IDs with hits)

    output:
    path "all_fasta_ids.txt"
    path "all_fasta_ids_sorted.txt"
    path "hit_fasta_ids_sorted.txt"
    path "failed_blast_${params.v_new}.txt"

    script:
    """
    # Extract all sequence IDs from fasta headers (remove '>')
    grep '^>' ${full_fasta} | sed 's/>//' > all_fasta_ids.txt

    # Sort the full list and the hit list of IDs
    sort all_fasta_ids.txt > all_fasta_ids_sorted.txt
    sort ${hit_ids_file} > hit_fasta_ids_sorted.txt

    # Find IDs present in all but not in hits (failed)
    comm -23 all_fasta_ids_sorted.txt hit_fasta_ids_sorted.txt > failed_blast_${params.v_new}.txt
    """
}




process BlastFilter {
    tag { "Filter BLAST results and select position" }
    //publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path blast_results
    path markers_file

    output:
    path "filtered_${params.v_new}.csv"
    path "filtered_failed.tsv"
    path "positions_${params.v_new}.csv"

    script:
    """
    python ${scripts}/blast_res_filter.py ${blast_results} filtered_${params.v_new}.csv filtered_failed.tsv
    python ${scripts}/select_position.py filtered_${params.v_new}.csv ${markers_file} positions_${params.v_new}.csv
    """
}

process CreateChip {
    tag { "Create final CHIP with positions in genome assembly version ${params.v_new}  ${new_ref.baseName}" }
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

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
    python ${scripts}/find_position.py ${positions_file} ${length_file} hits_${params.v_new}.txt
    python ${scripts}/concatinate.py hits_${params.v_new}.txt  ${markers_na} ${failed_blast} ${filtered_failed} chip_${params.v_new}.txt
    """
}




workflow {

    idx_ref = Index_genome(old_ref)
    obtained_fasta = Obtain_fasta(old_chip, old_ref, idx_ref)
    blast_db_files = MakeBlastDB(new_ref)
    blast_results = BlastSearch(obtained_fasta[0], blast_db_files)
    hit_ids = ExtractHitIds(blast_results)
    failed_blast = FindFailedBlastIds(
    obtained_fasta[0],
    hit_ids
    )
    filter_results = BlastFilter(
        blast_results,
        old_chip
    )

    CreateChip(
        filter_results[2],
        obtained_fasta[2],
        obtained_fasta[1],
        failed_blast[3],
        filter_results[1]
    )
}
