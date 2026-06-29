process ASM_STATS {
    label 'setup'
    label 'process_low'

    storeDir { params.tables_dir }

    input:
    path samples
    path genome_dir

    output:
    path 'asm_stats.tsv.gz', emit: stats
    path 'versions.yml',     emit: versions

    script:
    """
    set -euo pipefail

    TMPFILE=\$(mktemp)
    trap 'rm -f \$TMPFILE' EXIT

    printf 'ASMID\\ttotal_length_bp\\tN50_bp\\tcontig_count\\n' > \$TMPFILE

    awk -F',' 'NR>1 {print \$2}' ${samples} | sort -u | while read asmid; do
        [ -z "\$asmid" ] && continue
        asmid="\$(echo "\$asmid" | xargs)"

        if [ -f "${genome_dir}/\${asmid}.fa.gz" ]; then
            genome="${genome_dir}/\${asmid}.fa.gz"
        elif [ -f "${genome_dir}/\${asmid}.fa" ]; then
            genome="${genome_dir}/\${asmid}.fa"
        elif [ -f "${genome_dir}/\${asmid}.masked.fasta.gz" ]; then
            genome="${genome_dir}/\${asmid}.masked.fasta.gz"
        elif [ -f "${genome_dir}/\${asmid}.masked.fasta" ]; then
            genome="${genome_dir}/\${asmid}.masked.fasta"
        else
            echo "[WARN] No genome file found for \${asmid} in ${genome_dir}" >&2
            continue
        fi

        total_bp=\$(seqkit stats -T "\$genome" 2>/dev/null | tail -n 1 | awk '{print \$4}')
        n50=\$(seqkit fx2tab -l "\$genome" 2>/dev/null | sort -rn -k2 | \\
            awk -v total="\$total_bp" 'BEGIN{sum=0} {sum+=\$2; if(sum >= total/2) {print \$2; exit}}')
        contigs=\$(seqkit stats -T "\$genome" 2>/dev/null | tail -n 1 | awk '{print \$3}')

        [ -z "\$total_bp" ] && total_bp="0"
        [ -z "\$n50" ] && n50="0"
        [ -z "\$contigs" ] && contigs="0"

        printf '%s\\t%s\\t%s\\t%s\\n' "\$asmid" "\$total_bp" "\$n50" "\$contigs" >> \$TMPFILE
    done

    pigz -c \$TMPFILE > asm_stats.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    printf 'ASMID\\ttotal_length_bp\\tN50_bp\\tcontig_count\\n' | pigz -c > asm_stats.tsv.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: 2.8.0
    END_VERSIONS
    """
}
