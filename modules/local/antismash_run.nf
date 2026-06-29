process ANTISMASH_RUN {
    label 'antismash'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/antismash_local/**"), emit: results
    path 'versions.yml',                                    emit: versions

    script:
    def out = meta.id
    def gbk = "${params.target}/${out}/predict_results/${out}.gbk"
    """
    GBK="${gbk}"
    if [ ! -f "\$GBK" ] && [ -f "${gbk}.gz" ]; then
        zcat "${gbk}.gz" > ${out}.predict.gbk
        GBK=${out}.predict.gbk
    fi
    if [ ! -f "\$GBK" ]; then
        echo "ERROR: predict GBK not found: ${gbk}[.gz]" >&2
        exit 1
    fi
    source /etc/profile.d/modules.sh 2>/dev/null || true
    mkdir -p ${out}/antismash_local
    antismash --taxon ${params.antismash_taxon} \\
        --output-dir ${out}/antismash_local \\
        --genefinding-tool none \\
        --fullhmmer --clusterhmmer --cb-general --pfam2go \\
        -c ${task.cpus} \\
        \$GBK
    pigz ${out}/antismash_local/*.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version 2>&1 | grep -oP '(?<=antiSMASH )\\S+' || antismash --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/antismash_local
    touch ${out}/antismash_local/${out}.json.gz
    touch ${out}/antismash_local/index.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: 7.1.0
    END_VERSIONS
    """
}
