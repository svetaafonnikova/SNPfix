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

process Tsv2bed {
    tag{"Obtain bed file from chip"}
    publishDir = [
            path: "${params.outdir}",
            mode: 'copy', overwrite: true
        ]

    input:
    path chip

    output:
    path "markers_asm1_1based.bed"

    script:
    """
    awk -F'\t' 'BEGIN {OFS="\t"}
    NR==1 {
        for(i=1;i<=NF;i++) {
            if(\$i=="Chr") c_chr=i;
            if(\$i=="Position") c_pos=i;
            if(\$i=="MarkerID") c_id=i;
        }
        next
    }
    {
        # РџРµС‡Р°С‚СЊ: Chromosome, Start (Position - 1), End (Position), MarkerID
        print \$c_chr, (\$c_pos - 1), \$c_pos, \$c_id
    }' ${chip} > markers_asm1_1based.bed
    

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



process Crossmap {
    tag{"Obtain bed file with rustybam"}
    publishDir = [
            path: "${params.outdir}",
            mode: 'copy', overwrite: true
        ]

    input:
    path chain
    path bed

    output:
    path "chip_crossmap.bed"

    script:
    """
    CrossMap bed ${chain} ${bed}  chip_crossmap.bed

    """
}

process Bed2tsv {
    tag { "Make tsv chip from rustybam bed" }
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path bed

    output:
    path "chip_crossmap.tsv"

    script:
    """
    awk -F'\t' 'BEGIN {OFS="\t"; print "Name", "Chr", "Position"}
    {

        print \$4, \$1, \$3
    }' ${bed} > chip_crossmap.tsv
    """
}



workflow {

    paf = Alignment(old_ref, new_ref)
    chain=Paf2chain(paf)
    bedchip= Tsv2bed(old_chip)
    cross = Crossmap(chain, bedchip)
    newchip = Bed2tsv(cross)
}