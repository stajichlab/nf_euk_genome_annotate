// IPRSCAN5
process INTERPROSCAN_RUN {
    label 'interproscan'
    tag "${meta.id}"

    cpus   8
    memory '32 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}/annotate_misc/iprscan.xml")

    script:
    def out      = meta.id
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    if [ ! -f "${proteins}" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}" >&2
        exit 1
    fi
    mkdir -p ${out}/annotate_misc
    interproscan.sh -i ${proteins} -f XML -o ${out}/annotate_misc/iprscan.xml \\
        -dp -goterms -pa -t p -cpu ${task.cpus}
    """

    stub:
    def out = meta.id
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    """
}
