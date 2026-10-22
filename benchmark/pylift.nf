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



process Alignment {
    tag "Create paf file ${ref_genome.name}"
    
    publishDir = [
            path: "${params.reference}",
            mode: 'copy', overwrite: true
        ]


    input:
    path ref_genome
    path target_genome

    output:
    path "ref_asmb1_target_asmb2.paf"

    script:
    """
    minimap2 -c ${ref_genome} ${target_genome} > ref_asmb1_target_asmb2.paf
    """
}

process Paf2chain {
    tag{"Obtain chain file from chip"}
    publishDir = [
            path: "${params.reference}",
            mode: 'copy', overwrite: true
        ]

    input:
    path paf

    output:
    path "ref_asmb1_target_asmb2.chain"

    script:
    """
    paf2chain -i ${paf} > ref_asmb1_target_asmb2.chain

    """
}


process Pyliftover {
    tag{"Obtain txt markers file"}
    publishDir = [
            path: "${params.outdir}",
            mode: 'copy', overwrite: true
        ]

    input:
    path chip
    path chain

    output:
    path "chip_pylift.tsv"

    script:
    """
    python ${params.scripts}/pylift_script.py -i ${chip} -o chip_pylift.tsv -c ${chain}

    """
}



workflow {

    paf = Alignment(old_ref, new_ref)
    chain= Paf2chain(paf)
    out = Pyliftover(old_chip, chain)
}