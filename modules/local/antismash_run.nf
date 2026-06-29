process ANTISMASH_RUN {
    label 'antismash'
    tag "${meta.id}"

    cpus   8
    memory '16 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/antismash_local/**")

    script:
    def out = meta.id
    def gbk = "${params.target}/${out}/predict_results/${out}.gbk"
    """
    # Accept a compressed prediction (.gbk.gz); antismash needs it uncompressed, so
    # inflate a local copy in the work dir when only the gzipped form is present.
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
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/antismash_local
    touch ${out}/antismash_local/${out}.json.gz
    touch ${out}/antismash_local/index.html
    """
}
