process SIGNALP_RUN {
    label 'signalp'
    tag "${meta.id}"

    cpus   8
    memory '16 GB'
    time   '12h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/annotate_misc/signalp.results.txt")

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
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/signalp.results.txt
    """
}
