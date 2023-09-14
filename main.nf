#!/usr/bin/env nextflow

params.sampleSheet = ""
params.gbzIndex = ""
params.minIndex = ""
params.distIndex = ""
params.refFasta = ""
params.refPath = ""

process GET_REF_PATHS {
    input:
        path("index.gbz")
        val(refPath)

    output:
        path("ref_paths.txt")

    """
    vg paths -x index.gbz -S $refPath -L > ref_paths.txt
    """
}

process SNARLS {
    input:
        path("index.gbz")

    output:
        path("index.snarls")

    """
    vg snarls \
        -t $task.cpus \
        index.gbz \
        > index.snarls
    """
}

process GIRAFFE {
    input:
        path(gbzIndex)
        path(minIndex)
        path(distIndex)
        tuple val(sample), path(r1), path(r2)

    output:
        tuple val(sample), path("${sample}.gaf.gz")

    """
    cp $distIndex local.dist
    vg giraffe -t $task.cpus -p \
        -Z $gbzIndex \
        -m $minIndex \
        -d local.dist \
        -f $r1 -f $r2 \
        -N $sample -R $sample \
        -o gaf \
        | gzip > ${sample}.gaf.gz
    """
}

process MERGE_GAFS {
    input:
        tuple val(sample), path(allGafs, stageAs: "?/*")

    output:
        tuple val(sample), path("${sample}_merged.gaf.gz")

    shell:
        def gafsString = allGafs instanceof List ? allGafs.join(" ") : allGafs
        """
        cat $gafsString > ${sample}_merged.gaf.gz
        """
}

process PACK {
    input:
        path(gbzIndex)
        tuple val(sample), path("${sample}.gaf.gz")

    output:
        tuple val(sample), path("${sample}.pack")

    """
    vg pack \
        -x $gbzIndex \
        -a ${sample}.gaf.gz \
        -o ${sample}.pack \
        -Q 5 \
        -s 5 \
        -t $task.cpus
    """
}

process CALL {
    input:
        path(gbzIndex)
        path(refPaths)
        path(snarls)
        tuple val(sample), path("${sample}.pack")

    output:
        tuple val(sample), path("${sample}.vcf")

    """
    ref_path_args=""
    while read path; do
        ref_paths_args+="-p \$path "
    done < $refPaths

    vg call \
        $gbzIndex \
        -a \
        -k ${sample}.pack \
        -r $snarls \
        -s ${sample} \
        -t ${task.cpus} \
        \$ref_paths_args \
        > ${sample}.vcf
    """
}

process NORM {
    input:
        path(refFasta)
        val(refPath)
        tuple val(sample), path(inVcf)

    output:
        path("${sample}.norm.vcf.gz")

    """
    bcftools head $inVcf | awk '{
        if (/^##contig/) {
            if (/$refPath/) {
                sub(/ID=${refPath}#0#/, "ID=", \$0)
                print
            }
        } else print}' > new_header.txt
    bcftools reheader -h new_header.txt $inVcf > reheadered.vcf
    bcftools norm -f $refFasta reheadered.vcf \
        | bcftools sort -T . -m ${task.memory * 0.6} \
        | bgzip > ${sample}.norm.vcf.gz
    """
}

process MERGE_VCFS {
    input:
        path(allVcfs)

    output:
        path("merged.vcf.gz")

    """
    for vcf in $allVcfs; do
        bcftools index \$vcf
    done
    bcftools merge $allVcfs | bgzip > merged.vcf.gz
    """
}

workflow {
    gbzIndex = file(params.gbzIndex)
    minIndex = file(params.minIndex)
    distIndex = file(params.distIndex)
    refFasta = file(params.refFasta)
    refPath = file(params.refPath)

    GET_REF_PATHS(gbzIndex, refPath)

    Channel.fromPath(params.sampleSheet)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.sample_id,
            file(row.r1),
            file(row.r2)
        ) }
        .set { reads }

    reads
        .map { tuple(it[0], 1) }
        .groupTuple()
        .map { tuple(it[0], it[1].sum()) }
        .set { sampleLibraryCounts }


    SNARLS(gbzIndex)

    GIRAFFE(gbzIndex, minIndex, distIndex, reads)

    sampleLibraryCounts
        .cross(GIRAFFE.out)
        .map { tuple(
            groupKey(it[0][0], it[0][1]),
            it[1][1]
        ) }
        .groupTuple() | MERGE_GAFS

    PACK(gbzIndex, MERGE_GAFS.out)
    CALL(gbzIndex, GET_REF_PATHS.out, SNARLS.out, PACK.out)
    NORM(refFasta, refPath, CALL.out)
    MERGE_VCFS(NORM.out.collect())
}
