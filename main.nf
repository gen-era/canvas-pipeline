#!/usr/bin/env nextflow

import java.text.SimpleDateFormat
nextflow.enable.dsl=2

params.chip_id = ""
params.bpm = 'manifest-cluster/GSACyto_20044998_A1.bpm'
params.csv = 'gtcler/GSA-24v1-0_C1.csv'
params.egt = 'manifest-cluster/2003.egt'
params.fasta = '/staging/references/hg19/hg19.fa'
params.pfb = 'test/out.pfb'
params.hmm = '/path/to/default/hhall.hmm'
params.output_dir = 'outputs'
params.samplesheet = ''
params.band = "canvas-pipeline/"
params.tex_template = "asdadfa"
params.cnvs =  ""
params.protocol_id =  ""
params.institute = ""
params.sampleId = ""
params.position = ""
params.chip_type = ""
params.version = ""
params.include_plots = false

params.cnv_bed = ""
params.cnv_pk = ""
params.snap_probes = ""

idat_folder = "s3://canvas/chip_data/${params.chip_id}/idats"
output_dir = "s3://canvas/chip_data/${params.chip_id}"

workflow {
    if (params.cnv_bed) {
        sampleId = "${params.chip_id}_${params.position}"
        band = Channel.fromPath(params.band)
        bed = Channel.fromPath(params.cnv_bed)
        cnv_input = bed.map { bed -> [sampleId, params.cnv_pk, bed] }
        sample_summary = Channel.fromPath("s3://canvas/chip_data/${params.chip_id}/gtcs/gt_sample_summary.csv")

        plot_dir = Channel.fromPath("s3://canvas/chip_data/${params.chip_id}/plots/${sampleId}")
        baf = "s3://canvas/chip_data/${params.chip_id}/bedgraphs/${sampleId}.BAF.bedgraph.gz"
        lrr = "s3://canvas/chip_data/${params.chip_id}/bedgraphs/${sampleId}.LRR.bedgraph.gz"

        makesexfile(sample_summary)

        if (params.snap_probes.toBoolean()) {
            snap_probes(cnv_input, lrr)
            cnv_input = snap_probes.out.bed.map { bed -> [sampleId, params.cnv_pk, bed] }
        }

        cnv_addiscn(cnv_input, band)
        cnv_classification_bed(cnv_input, makesexfile.out.first())
        cnv_classification(cnv_classification_bed.out)


        bed_with_baf_lrr = cnv_classification_bed.out.map { bed -> [sampleId, params.cnv_pk, baf, lrr, bed[2]] }
        cnv_makeplots(bed_with_baf_lrr)

        cnv_make_bedgraphs(cnv_classification_bed.out)
    } else if (params.cnvs) {
        cnvs = Channel.fromPath(params.cnvs)
        println("cnvs")
        tex_template = Channel.fromPath(params.tex_template)
        sampleId = "${params.chip_id}_${params.position}"
        plot_dir = Channel.fromPath("s3://canvas/chip_data/${params.chip_id}/plots/${sampleId}")

        println(sampleId)
        maketemplatefromcnv(
            cnvs,
            sampleId,
            plot_dir,
            params.protocol_id,
            params.institute,
            params.chip_type,
            params.version,
            params.include_plots,
            tex_template
        )
        println(cnvs)
        template_dir = file(params.tex_template).getParent()

        makereport(maketemplatefromcnv.out, template_dir)
    } else {

    idat_folder = Channel.fromPath(idat_folder)
    bpm = Channel.fromPath(params.bpm)
    csv = Channel.fromPath(params.csv)
    egt = Channel.fromPath(params.egt)
    pfb = file(params.pfb)
    hmm = file(params.hmm)
    fasta = Channel.fromPath(params.fasta)
    band = Channel.fromPath(params.band)
    tex_template = Channel.fromPath(params.tex_template)

    fasta_dir = file(params.fasta).getParent()
    template_dir = file(params.tex_template).getParent()

    idat2gtc(idat_folder, bpm, egt)
    makesexfile(idat2gtc.out.sample_summary)
    fasta = Channel.fromPath(params.fasta)
    fasta_index = Channel.fromPath("${params.fasta}.fai")
    prepared_reference = prepare_reference(fasta, fasta_index)
    vcfs = gtc2vcf(idat2gtc.out.gtc_dir, bpm, csv, prepared_reference)

    vcfs .flatten() .map {
        vcf -> [vcf.getSimpleName(), vcf]
    }
    .set {named_vcfs}

    vcf2penncnv(named_vcfs)

    makebedgraphs(vcf2penncnv.out)

    penncnv_detect(vcf2penncnv.out, pfb, hmm, makesexfile.out.first())
    penncnv_clean_cnv(penncnv_detect.out, pfb)

    beds = penncnv2bed(penncnv_clean_cnv.out, makesexfile.out.first())
    classification(penncnv2bed.out)
    make_cnv_bedgraphs(beds)
    smooth_lrr(makebedgraphs.out)

    addiscn(penncnv_clean_cnv.out, band.first())

    bedgraphswithbeds = makebedgraphs.out.combine(beds, by:0)
    makeplots(bedgraphswithbeds)

    // cnvswithscores = classification.out.scoresheet
    //     .combine(addiscn.out, by:0)
    //     .combine(makeplots.out, by:0)


    // cnvswithscores.view()
    // if ( params.samplesheet ) {
    //     samplesheet = Channel.fromPath(params.samplesheet)
    //     samplesheet.view()

    //     samplesheet
    //         .splitCsv(header:["sample_id", "protocol_id", "institute"],
    //             sep:"\t",
    //             skip:1
    //         )
    //         | map { row ->
    //             [row.sample_id, row.protocol_id, row.institute]
    //         }
    //         | set {samplesheet}
    //     samplesheet.view()
    //     cnvswithscores = cnvswithscores.combine(samplesheet, by:0)
    // }
    // else {
    //     cnvswithscores = cnvswithscores.map { _ ->
    //         _ + [_[0], 'noinstitute']
    //     }
    // }
    // cnvswithscores.view()
    //maketemplate(cnvswithscores, tex_template)
    //makereport(maketemplate.out, template_dir)
    }
}

// Process to cache and prepare reference genome
process prepare_reference {
    storeDir "${workflow.workDir}/reference_cache"

    input:
    path fasta
    path fai

    output:
    tuple path("${fasta}"), path("${fasta}.fai")  // Output both fasta and its index

    script:
    """
    echo "staging reference"
    """
}


process makesexfile {
    memory "1 GB"
    cpus 1
    input:
    path sample_summary

    output:
    path "sex_file"

    shell:
    '''
    awk -F"," '{printf "%s\\t%s\\n", $1, $7}' !{sample_summary} > sex_file
    '''
}

process idat2gtc {
    publishDir "${output_dir}/", mode: "copy"

    input:
    path(idat_folder)
    path(bpm)
    path(egt)

    output:
    path "gtcs/", emit: gtc_dir
    path "gtcs/gt_sample_summary.csv", emit: sample_summary

    script:
    """
    mkdir gtcs/
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena \
    genotype call \
    --bpm-manifest ${bpm} \
    --cluster-file ${egt} \
    --idat-folder ${idat_folder} \
    --output-folder gtcs/
    """
}

process gtc2vcf {
    publishDir "${output_dir}/vcfs", mode: "copy"

    input:
    path gtcs
    path bpm
    path csv
    tuple path(fasta), path(fasta_idx)  // Take prepared reference as input

    output:
    path("*.vcf.gz")

    script:
    """
    /opt/dragena-linux-x64-DAv1.0.0/dragena/dragena genotype gtc-to-vcf \
        --gtc-folder ${gtcs} \
        --bpm-manifest ${bpm} \
        --csv-manifest ${csv} \
        --genome-fasta-file "${fasta}" \
        --output-folder ./
    """
}

process vcf2penncnv {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "staphb/bcftools:1.20"
    publishDir "${output_dir}/penncnv_inputs/", mode: "copy"

    input:
    tuple val(sampleId), path(vcf)

    output:
    tuple val(sampleId), path("${sampleId}.txt")

    script:
    """
    printf "Name\tChr\tPosition\tGType\tLog R Ratio\tB Allele Freq\n" > "${sampleId}.txt"
    bcftools query -f '%ID\t%CHROM\t%POS[\t%GT\t%LRR\t%BAF]\n' ${vcf} | sed 's/1\\/1/BB/g;s/0\\/0/AA/g;s/0\\/1/AB/g;s/.\\/./NC/g' >> "${sampleId}.txt"
    """
}

process makebedgraphs {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/bedgraphs", mode: "copy"

    input:
    tuple val(sampleId), path(txt)

    output:
    tuple val(sampleId), path("*BAF.bedgraph.gz"), path("*LRR.bedgraph.gz")


    shell:
    '''
    awk '{printf "%s\\t%s\\t%s\\t%s\\n", $2, $3, $3, $5}' !{txt} | sed 1d | sed '/\\.$/d' | gzip -c > !{sampleId}.LRR.bedgraph.gz
    awk '{printf "%s\\t%s\\t%s\\t%s\\n", $2, $3, $3, $6}' !{txt} | sed 1d | sed '/\\.$/d' | gzip -c > !{sampleId}.BAF.bedgraph.gz
    '''
}


process smooth_lrr {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/bedgraphs", mode: "copy"

    input:
    tuple val(sampleId), path(baf), path(lrr)

    output:
    tuple val(sampleId), path ("${sampleId}.LRR_smooth.bedgraph.gz")

    script:
    """
    zcat ${lrr} > lrr_file
    python3 ${projectDir}/smooth_lrr.py --input_file lrr_file --output_file ${sampleId}.LRR_smooth.bedgraph
    rm lrr_file
    gzip ${sampleId}.LRR_smooth.bedgraph
    """
}

process penncnv_detect {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "genomicslab/penncnv"
    publishDir "${output_dir}/cnvs", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path pfb
    path hmm
    path sex_file

    output:
    tuple val(sampleId), path("${sampleId}.cnv.txt")

    script:
    """
    grep "${sampleId}" "${sex_file}" > sample_sex
    sex=\$(cut -f2 sample_sex)
    /home/user/PennCNV/detect_cnv.pl -test -hmm ${hmm} -pfb ${pfb} ${txt} --confidence         -log ${sampleId}.log      -out ${sampleId}.autosomal.out.cnv
    
    if [[ "\$sex" == "M" || "\$sex" == "F" ]]; then
        /home/user/PennCNV/detect_cnv.pl -test -sexfile sample_sex -hmm ${hmm} -pfb ${pfb} ${txt} --confidence --chrx  -log ${sampleId}.chrx.log -out ${sampleId}.chrx.out.cnv
        /home/user/PennCNV/detect_cnv.pl -test -sexfile sample_sex -hmm ${hmm} -pfb ${pfb} ${txt} --confidence --chry  -log ${sampleId}.chry.log -out ${sampleId}.chry.out.cnv
    else
        /home/user/PennCNV/detect_cnv.pl -test -hmm ${hmm} -pfb ${pfb} ${txt} --confidence --chrx  -log ${sampleId}.chrx.log -out ${sampleId}.chrx.out.cnv
        /home/user/PennCNV/detect_cnv.pl -test -hmm ${hmm} -pfb ${pfb} ${txt} --confidence --chry  -log ${sampleId}.chry.log -out ${sampleId}.chry.out.cnv
    fi

    cat ${sampleId}.autosomal.out.cnv ${sampleId}.chrx.out.cnv ${sampleId}.chry.out.cnv > ${sampleId}.cnv.txt
    """
}

process penncnv_clean_cnv {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "genomicslab/penncnv"
    publishDir "${output_dir}/cnvs", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path pfb

    output:
    tuple val(sampleId), path("${sampleId}.cleaned.txt")


    script:
    """
    /home/user/PennCNV/clean_cnv.pl combineseg --fraction 0.5 --signalfile ${pfb} ${txt} > ${sampleId}.cleaned.txt
    """
}

process penncnv2bed {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/beds", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path sex_file

    output:
    tuple val(sampleId), path ("${sampleId}.bed")

    script:
    """
    python3 ${projectDir}/penncnv2bed.py --input_file ${txt} --output_file ${sampleId}.bed --sex_file ${sex_file}
    """
}


process make_cnv_bedgraphs {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/bedgraphs", mode: "copy"

    input:
    tuple val(sampleId), path(bed)

    output:
    tuple val(sampleId), path ("*.bedgraph.gz")

    script:
    """
    awk -F"\\t" '{printf "%s\\t%s\\t%s\\t1\\n", \$1, \$2, \$3}' "$bed" | gzip -c > "${sampleId}.CNV_pos.bedgraph.gz"
    awk -F"\\t" '{printf "%s\\t%s\\t%s\\t-1\\n", \$1, \$2, \$3}' "$bed" | gzip -c > "${sampleId}.CNV_neg.bedgraph.gz"
    """
}


process classification {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "ghcr.io/gen-era/classifycnv:3b6863daf18af577e3090f31fa24aad0403d5027"
    publishDir "${output_dir}/ClassifyCNV/", mode: "copy"

    input:
    tuple val(sampleId), path(bed)

    output:
    tuple val(sampleId), path("${sampleId}/Scoresheet.txt"), emit: scoresheet
    path "${sampleId}/Intermediate_files/*.bed"

    script:
    """
    python3 /ClassifyCNV/ClassifyCNV.py --infile ${bed} --GenomeBuild hg19 --precise --outdir "${sampleId}"
    """
}


process addiscn {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/cnvs/", mode: "copy"

    input:
    tuple val(sampleId), path(txt)
    path(band)

    output:
    tuple val(sampleId), path("${sampleId}_iscn.txt")

    script:
    """
    python3 ${projectDir}/getISCN.py --cnv ${txt} --band ${band} --output "${sampleId}_iscn.txt"
    """
}

process makeplots {
    tag "$sampleId"
    container "ghcr.io/gen-era/bedgraph-visualizer:70caf7a190d65bed88aeb8a9c552b37c7580a170"
    publishDir "${output_dir}/plots", mode: "copy"

    input:
    tuple val(sampleId), path(BAF), path(LRR), path(bed)

    output:
    tuple val(sampleId), path("${sampleId}/")

    script:
    """
    sed 's/^chr//' ${bed} > regions
    Rscript /usr/src/app/bedgraph-visualizer.R region_plot ${BAF} ${LRR} regions ${sampleId}
    Rscript /usr/src/app/bedgraph-visualizer.R genome_plot ${BAF} ${LRR} ${sampleId}
    """
}

process maketemplate {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/texs", mode: "copy"

    input:
    tuple val(sampleId), path(scoresheet), path(cnv), path(plot_dir), val(protocol_id), val(institute)
    path tex_template

    output:
    tuple val(sampleId), path("${sampleId}.tex"), path(plot_dir), val(protocol_id)

    script:
    """
    python3 ${projectDir}/tex_template_compile.py \
        --scoresheet_file ${scoresheet} \
        --cnv_file ${cnv} \
        --tex_template ${tex_template} \
        --output_file "${sampleId}.tex" \
        --plot_dir ${plot_dir} \
        --protocol_id "${protocol_id}" \
        --institute "${institute}"
    """
}

process maketemplatefromcnv {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/texs", mode: "copy"

    input:
    path cnvs 
    val sampleId 
    path plot_dir 
    val protocol_id 
    val institute 
    val chip_type 
    val version 
    val include_plots 
    path tex_template

    output:
    tuple val(sampleId), path("${sampleId}.tex"), path(plot_dir), val(protocol_id)

    script:
    def includePlotsArg = include_plots ? '--include_plots' : ''

    """
    python3 ${projectDir}/tex_template_compile.py \
        --tex_template ${tex_template} \
        --sample_id ${sampleId} \
        --output_file "${sampleId}.tex" \
        --plot_dir ${plot_dir} \
        --cnvs ${cnvs} \
        --protocol_id "${protocol_id}" \
        --institute "${institute}" \
        --chip_type "${chip_type}" \
        --version "${version}" \
        ${includePlotsArg}
    """
}

process makereport {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "texlive/texlive"
    publishDir "${output_dir}/pdfs", mode: "copy"

    input:
    tuple val(sampleId), path(tex), path(plot_dir), val(protocol_id)
    path template_dir

    output:
    tuple val(sampleId), path("${sampleId}_${protocol_id}_${timestamp}.pdf")

    script:
    Date now = new Date();
    SimpleDateFormat timestamp_formatter = new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
    timestamp = timestamp_formatter.format(now)

    """
    cp "${template_dir}"/* ./
    pdflatex ${tex} --jobname="${protocol_id}_${timestamp}.pdf"
    pdflatex ${tex} --jobname="${protocol_id}_${timestamp}.pdf"
    mv "${sampleId}.pdf" "${sampleId}_${protocol_id}_${timestamp}.pdf"
    """
}


process cnv_addiscn {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/cnvs/", mode: "copy"

    input:
    tuple val(sampleId), val(cnv_pk), path(txt)
    path(band)

    output:
    tuple val(sampleId), val(cnv_pk), path("${sampleId}_${cnv_pk}_iscn.txt")

    script:
    """
    python3 ${projectDir}/getISCN.py --bed ${txt} --band ${band} --output "${sampleId}_${cnv_pk}_iscn.txt"
    """
}

process cnv_classification_bed {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/beds", mode: "copy"

    input:
    tuple val(sampleId), val(cnv_pk), path(txt)
    path(sex_file)

    output:
    tuple val(sampleId), val(cnv_pk), path("${sampleId}_${cnv_pk}.bed")

    script:
    """
    python3 ${projectDir}/penncnv2bed.py --input_file ${txt} --sample_id ${sampleId} --format bed --output_file ${sampleId}_${cnv_pk}.bed --sex_file ${sex_file}
    """
}   

process cnv_classification {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    container "fauzul/classifycnv:1.0"
    publishDir "${output_dir}/ClassifyCNV/", mode: "copy"

    input:
    tuple val(sampleId), val(cnv_pk), path(bed)

    output:
    tuple val(sampleId), val(cnv_pk), path("${sampleId}_${cnv_pk}/Scoresheet.txt"), emit: scoresheet
    path "${sampleId}_${cnv_pk}/Intermediate_files/*.bed"

    script:
    """
    python3 /ClassifyCNV/ClassifyCNV.py --infile ${bed} --GenomeBuild hg19 --precise --outdir "${sampleId}_${cnv_pk}"
    """
}


process cnv_makeplots {
    tag "$sampleId"
    container "yserdem/bedgraph-visualizer"
    publishDir "${output_dir}/plots", mode: "copy"

    input:
        tuple val(sampleId), val(cnv_pk), path(BAF), path(LRR), path(bed)

    output:
     tuple val(sampleId), path("${sampleId}/")

    script:
    """
    sed 's/^chr//' ${bed} > regions
    Rscript /usr/src/app/bedgraph-visualizer.R region_plot ${BAF} ${LRR} regions ${sampleId}
    """
}

process cnv_make_bedgraphs {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    publishDir "${output_dir}/bedgraphs", mode: "copy"

    input:
    tuple val(sampleId), val(cnv_pk), path(bed)

    output:
    tuple val(sampleId), path ("${sampleId}_${cnv_pk}.CNV_pos.bedgraph.gz"), path ("${sampleId}_${cnv_pk}.CNV_neg.bedgraph.gz")

    script:
    """
    awk -F"\\t" '{printf "%s\\t%s\\t%s\\t1\\n", \$1, \$2, \$3}' "$bed" | gzip -c > "${sampleId}_${cnv_pk}.CNV_pos.bedgraph.gz"
    awk -F"\\t" '{printf "%s\\t%s\\t%s\\t-1\\n", \$1, \$2, \$3}' "$bed" | gzip -c > "${sampleId}_${cnv_pk}.CNV_neg.bedgraph.gz"
    """
}

process snap_probes {
    memory "1 GB"
    cpus 1
    tag "$sampleId"
    input:
    tuple val(sampleId), val(cnv_pk), path(bed)
    path(lrr)

    output:
    path("snapped.bed"), emit: bed

    script:
    """
    python3 ${projectDir}/snap_probes.py --bed ${bed} --lrr ${lrr} --output snapped.bed
    """
}
