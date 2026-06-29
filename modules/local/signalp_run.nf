process SIGNALP_RUN {
    label 'signalp'
    label 'process_medium'
    tag "${meta.id}"

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/annotate_misc/signalp.results.txt"), emit: results
    path 'versions.yml',                                                    emit: versions

    script:
    def out      = meta.id
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    TMPDIR=\${SCRATCH:-/tmp}
    signalp6 -od \$TMPDIR/${out}_signalp \\
        -org euk --mode fast -format txt \\
        -fasta ${proteins} \\
        --write_procs ${task.cpus} -bs 16
    mkdir -p ${out}/annotate_misc
    cp \$TMPDIR/${out}_signalp/prediction_results.txt ${out}/annotate_misc/signalp.results.txt
    rm -rf \$TMPDIR/${out}_signalp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(signalp6 --version 2>&1 | grep -oP '\\d+\\.\\d+\\S*' | head -1)
    END_VERSIONS
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/signalp.results.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: 6.0g
    END_VERSIONS
    """
}
