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
        # Печать: Chromosome, Start (Position - 1), End (Position), MarkerID
        print \$c_chr, (\$c_pos - 1), \$c_pos, \$c_id
    }' ${chip} > markers_asm1_1based.bed
    

    """
}




process Rustybam {
    tag{"Obtain bed file with rustybam"}
    publishDir = [
            path: "${params.outdir}",
            mode: 'copy', overwrite: true
        ]

    input:
    path paf
    path bed

    output:
    path "asm2_coordinates_rustybam.bed"

    script:
    """
    rustybam liftover --bed  ${bed}  ${paf} > asm2_coordinates_rustybam.bed

    """
}

process Bed2tsv {
    tag { "Make tsv chip from rustybam bed" }
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path rusty

    output:
    path "chip_rustybam.tsv"

    script:
    """
    awk -F'\t' -v OFS='\t' '
BEGIN { print "Name", "Chr", "Position" }
/^#/ || /^[[:space:]]*\$/ { next }
{
    target_chr = \$1
    target_pos = \$4
    marker_name = ""
    for (i = 13; i <= NF; i++) {
        if (\$i ~ /^id:Z:/) {
            marker_name = substr(\$i, 6)
            break
        }
    }
    if (marker_name != "") print marker_name, target_chr, target_pos
}' ${rusty} > chip_rustybam.tsv
    """
}



workflow {

    paf = Alignment(old_ref, new_ref)
    bedchip= Tsv2bed(old_chip)
    bed = Rustybam(paf, bedchip)
    newchip = Bed2tsv(bed)
}