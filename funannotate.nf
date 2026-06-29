#!/usr/bin/env nextflow

/*
 * SOURCE: ../../../1KFG/common_annotate/pipeline/nextflow/funannotate.nf
 * Last synced: 2026-05-23
 * Changes vs source: removed nextflow.enable.dsl=2; params block moved to
 *                    conf/profile_annotate.config.
 *
 * Usage (from project root — a pipeline profile is REQUIRED; without it
 * params.taxondb / params.funannotate_db are null and parsing fails):
 *   sbatch nextflow/run_annotate.sh
 *   nextflow run nextflow/funannotate.nf -c nextflow/nextflow.config \
 *       -profile annotate,slurm,ucr_hpcc -resume
 */

// Data contract: every channel element is `tuple val(meta), val/path(genome)`.
// meta is a Map built by SampleUtils.makeMeta(row) — see lib/SampleUtils.groovy.
//   meta.id is the ONLY field used for tag{} and file naming.
//   meta.asmid, meta.species, meta.strain, meta.locustag, meta.busco,
//   meta.transl_table, meta.taxonid carry payload used inside process scripts.
//   header_length is NOT in meta — it comes from params.header_length (default 24).
//
// GENOME_CLEAN receives: tuple val(meta), path(genome_gz), val(taxondb)
//   → emits: tuple val(meta), path(genome_fa)   [storeDir writes input_clean_genomes/<asmid>.fa.gz]
// SRA_FETCH receives: val(species_tag), val(taxonid)   [only when --run_sra_fetch]
//   → emits: val(species_tag), path(norm_R1.fastq.gz), path(norm_R2.fastq.gz), path(se)
// RNASEQ_PREPARE receives: tuple val(species_tag), val(meta), val(genome_fa), path(r1), path(r2), path(se)
//   → emits: val(species_tag), path(trinity-GG.fasta)   [storeDir caches in rnaseq_data/]
// FUNANNOTATE_TRAIN receives: tuple val(meta), val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa)
//   → emits: tuple val(meta), val(genome_fa)
// FUNANNOTATE_PREDICT receives: tuple val(meta), val(genome_fa)
//   → emits metadata: tuple val(meta)

// Run funannotate train on the representative (first) assembly of each species, then
// archive the Trinity-GG transcripts (normalized reads are in rnaseq_reads)
// reads into rnaseq_data/ so all other strains can skip those expensive steps.
// storeDir skips this process entirely if all five output files already exist.
process RNASEQ_PREPARE {
    label 'funannotate'
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_data"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(species_tag), val(meta), val(genome_fa), path(r1), path(r2), path(se)

    output:
    tuple val(species_tag),
            path("${species_tag}.trinity-GG.fasta"), emit: shared

    script:
    def out           = meta.id
    def species       = meta.species
    def strain        = meta.strain
    def header_length = params.header_length
    """
    # ── Empty-reads sentinel: no RNA-seq found by SRA_FETCH / SRA_FETCH_SE ──
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ]; then
        echo "[INFO] No RNAseq reads for ${species_tag}; writing empty shared markers"
        touch ${species_tag}.trinity-GG.fasta
        exit 0
    fi

    # ── If representative was already trained, just extract shared files ──────
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    if [ -f "\$TRAIN_GFF3" ]; then
        echo "[INFO] Training already complete for ${out}; extracting shared files to rnaseq_data"
        TRAINDIR="${params.training_target}/${out}/training"
        TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
        if [ -n "\$TRINITY_FA" ]; then
            cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
        else
            touch ${species_tag}.trinity-GG.fasta
        fi
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    # ── Run full funannotate train on the representative genome ───────────────
    # Use SCRATCH for the funannotate output dir so Trinity/HISAT2/normalize
    # intermediates land on fast local storage and don't consume project quota.
    echo "[INFO] RNASEQ_PREPARE: running funannotate train for representative ${out} (species: ${species_tag})"

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    if [ -s "${r1}" ]; then
        funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    else
        echo "[INFO] RNASEQ_PREPARE: using single-end reads for ${out}"
        funannotate train -i "\$GENOME_IN" -o \$SCRATCH/${out} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            --stop_after_trinity --no_trimmomatic
    fi

    # ── Copy shared outputs to rnaseq_data/ ──────────────────────────────────
    TRAINDIR="\$SCRATCH/${out}/training"
    TRINITY_FA=\$(find \$TRAINDIR -maxdepth 1 -name "trinity.fasta" | head -1)
    if [ -n "\$TRINITY_FA" ]; then
        cp "\$TRINITY_FA" ${species_tag}.trinity-GG.fasta
    else
        echo "[WARN] No trinity.fasta found under \$TRAINDIR for ${out}"
        touch ${species_tag}.trinity-GG.fasta
    fi

    # ── Clean up scratch output dir (all intermediates were temporary) ────────
    rm -rf "\$SCRATCH/${out}"
    echo "[INFO] RNASEQ_PREPARE complete for ${species_tag}"
    """

    stub:
    def out = meta.id
    """
    echo ">stub_trinity_${species_tag}" > ${species_tag}.trinity-GG.fasta
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}

// For non-representative strains: funannotate train --trinity <shared_fasta> runs only
// PASA (skips Trimmomatic, normalization, HISAT2, and Trinity-GG assembly).
// Falls back to a full train when no shared Trinity is available (e.g. species with
// a single strain or when run_sra_fetch is false).
process FUNANNOTATE_TRAIN {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '96 GB'
    time   '120h'

    input:
    tuple val(meta), val(genome_fa), path(r1), path(r2), path(se), path(trinity_fa)

    output:
    tuple val(meta), val(genome_fa)

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def header_length = params.header_length
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no RNA-seq data at all ────────────────────────────────────────
    if [ ! -s "${r1}" ] && [ ! -s "${se}" ] && [ ! -s "${trinity_fa}" ]; then
        echo "[INFO] No RNAseq data for ${out}, skipping funannotate train"
        exit 0
    fi

    # ── Skip if training output already present and rnaseq is not newer than GBK ──
    # Accept a compressed prediction (.gbk.gz) as "done" so folders can be space-saved.
    TRAIN_GFF3="${params.training_target}/${out}/training/funannotate_train.pasa.gff3"
    PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk"
    [ -f "\$PREDICT_GBK" ] || PREDICT_GBK="${params.target}/${out}/predict_results/${out}.gbk.gz"
    if [ -f "\$TRAIN_GFF3" ]; then
        RETRAIN=0
        if [ -f "\$PREDICT_GBK" ]; then
            # Re-train if the rnaseq reads are newer than the existing prediction GBK.
            if [ -s "${r1}" ] && [ "${r1}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq R1 reads newer than predict GBK for ${out}; retraining"
#                rm -rf "${params.training_target}/${out}/training"
                RETRAIN=1
            elif [ -s "${se}" ] && [ "${se}" -nt "\$PREDICT_GBK" ]; then
                echo "[INFO] RNAseq SE reads newer than predict GBK for ${out}; retraining"
#                rm -rf "${params.training_target}/${out}/training"
                RETRAIN=1
            fi
        fi
        if [ \$RETRAIN -eq 0 ]; then
            echo "[INFO] Training already complete for ${out}; skipping"
            exit 0
        fi
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load funannotate/dev-1.8.18 
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.training_target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/mysql/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        # ──  may be unnecessary if overridden by -B option later? ──
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH/mysql_db
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Inflate a gzipped clean genome to a local uncompressed copy; funannotate cannot
    # read a gzipped FASTA via -i. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Use shared Trinity transcripts (PASA only) or run full train ──────────
    if [ -s "${trinity_fa}" ]; then
        if [ -s "${r1}" ]; then
            echo "[INFO] Running funannotate train (PASA+PE) for ${out} using shared Trinity"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        elif [ -s "${se}" ]; then
            echo "[INFO] Running funannotate train (PASA+SE) for ${out} using shared Trinity"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --single_norm ${se} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        else
            echo "[INFO] Running funannotate train (PASA only, no reads) for ${out} using shared Trinity"
            funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
                --trinity ${trinity_fa} --left_norm ${r1} --right_norm ${r2} \\
                --species "${species}" --strain "${strain}" \\
                --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
                --header_length ${header_length} \\
                --jaccard_clip --no-progress \\
                --max_intronlen ${params.max_intronlen} \\
                \$pasa_db_arg
        fi
    elif [ -s "${r1}" ]; then
        echo "[INFO] Running funannotate train (full PE, no shared Trinity) for ${out}"
        funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --left_norm ${r1} --right_norm ${r2} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --jaccard_clip --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    else
        echo "[INFO] Running funannotate train (full SE, no shared Trinity) for ${out}"
        funannotate train -i "\$GENOME_IN" -o ${params.training_target}/${out} \\
            --single_norm ${se} --aligners minimap2 \\
            --species "${species}" --strain "${strain}" \\
            --cpus ${task.cpus} --memory ${task.memory.toGiga()}G \\
            --header_length ${header_length} \\
            --no-progress --min_coverage 4 \\
            --max_intronlen ${params.max_intronlen} \\
            \$pasa_db_arg
    fi

    # ── Remove large intermediates not needed for predict or update ─────────────
    # Keeps: *.bam, *.bai, *.pasa.gff3, *.stringtie.gtf, *.transcripts.gff3
    TRAINDIR="${params.training_target}/${out}/training"
    echo "[INFO] Removing large training intermediates in \$TRAINDIR"
    rm -rf "\$TRAINDIR/hisat2"
    rm -rf "\$TRAINDIR/trinity_gg"
    echo "[INFO] Training cleanup complete for ${out}"
    echo "mysql is ${params.pasa_mysql}"
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] FUNANNOTATE_TRAIN stub for ${out}"
    mkdir -p ${params.training_target}/${out}/training
    touch ${params.training_target}/${out}/training/funannotate_train.pasa.gff3
    """
}

// Option B persistence model: funannotate predict computes DIRECTLY into the persistent
// per-genome dir (${params.target}/${out}), symmetric with FUNANNOTATE_TRAIN writing to
// training_target. funannotate checkpoints into predict_misc/, so a restart after an
// OOM/timeout/orchestrator death resumes completed steps in place rather than starting
// over. There is no publishDir copy and no work-dir<->target rsync: the durable output is
// written where downstream steps already read it. Large intermediates still go to the
// node-local --tmpdir. The Nextflow output is a small marker file (nothing consumes the
// predict dir as a channel; downstream rebuilds metadata from the CSV and gates on the
// on-disk GBK), so emitting a marker keeps the DAG edge without copying the result tree.
process FUNANNOTATE_PREDICT {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '32 GB'
    time   '32h'

    input:
    tuple val(meta), val(genome_fa)

    output:
    val meta, emit: metadata
    path("${meta.id}.predict.done"), emit: done

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def transl_table  = meta.transl_table
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load funannotate/dev-1.8.18
    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    PREDICTDIR="${params.target}/${out}"
    PREDICT_GBK="\$PREDICTDIR/predict_results/${out}.gbk"

    if [ "${params.debug.toBoolean()}" = "true" ]; then
        echo "[DEBUG] out=${out} asmid=${asmid} species=${species} strain=${strain}"
        echo "[DEBUG] locustag=${locustag} busco=${busco_lineage} transl_table=${transl_table}"
        echo "[DEBUG] proteins=${params.proteins} genome_fa=${genome_fa}"
        echo "[DEBUG] PREDICTDIR=\$PREDICTDIR TMPDIR=\$TMPDIR pwd=\$(pwd)"
    fi

    # ── Skip vs. refresh decision ─────────────────────────────────────────────
    # The workflow schedules this process when the GBK is missing OR stale (rnaseq/trinity
    # newer than the GBK, per staleRnaseq()). Re-derive staleness here from the same on-disk
    # timestamps so a current GBK short-circuits, but a stale one forces a clean re-predict.
    if [ -s "\$PREDICT_GBK" ]; then
        SPECIES_TAG=\$(printf '%s' "${species}" | sed -E 's/[[:space:]]+/_/g')
        STALE=0
        for f in "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_R1.fastq.gz" \\
                 "${launchDir}/rnaseq_reads/\${SPECIES_TAG}_norm_SE.fastq.gz" \\
                 "${launchDir}/rnaseq_data/\${SPECIES_TAG}.trinity-GG.fasta"; do
            if [ -s "\$f" ] && [ "\$f" -nt "\$PREDICT_GBK" ]; then STALE=1; fi
        done
        if [ "\$STALE" -eq 0 ]; then
            echo "[INFO] Prediction already complete and current for ${out}; nothing to do"
            touch ${out}.predict.done
            exit 0
        fi
        echo "[INFO] Stale prediction for ${out}: rnaseq/trinity newer than GBK — clearing predict outputs for a fresh run"
        rm -rf "\$PREDICTDIR/predict_results" "\$PREDICTDIR/predict_misc"
    fi

    mkdir -p "\$PREDICTDIR"

    # ── Guard against a corrupt partial from a previous attempt ───────────────
    # funannotate resumes from predict_misc/. If predict_results/ exists without a
    # predict_misc/ (a half-written tree with no checkpoints and no GBK), clear it so
    # predict starts the consensus/output step from a clean state instead of choking on it.
    if [ ! -d "\$PREDICTDIR/predict_misc" ] && [ -d "\$PREDICTDIR/predict_results" ]; then
        echo "[WARN] predict_results/ present without predict_misc/ for ${out}; clearing stale partial"
        rm -rf "\$PREDICTDIR/predict_results"
    fi

    # funannotate predict expects training data at <outdir>/training; point it at the
    # persistent training dir. The symlink lives in the persistent project tree (no
    # publishDir to recursively copy the target), so it is left in place.
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "\$PREDICTDIR/training"
    fi

    TBL2ASN_PARAMS="-l paired-ends"

    # Inflate a gzipped clean/masked genome to a local uncompressed copy; funannotate
    # cannot read a gzipped FASTA via -i, and the pre-flight awk below also needs plain
    # text. Plain (uncompressed) genomes pass through unchanged.
    GENOME_FA="${genome_fa}"
    case "\$GENOME_FA" in
        *.gz) echo "[INFO] Inflating compressed genome \$GENOME_FA"; pigz -dc "\$GENOME_FA" > genome_input.fa; GENOME_IN="\$(pwd)/genome_input.fa" ;;
        *)    GENOME_IN="\$GENOME_FA" ;;
    esac

    # ── Too-small-genome pre-flight guard ────────────────────────────────────
    # Assemblies that are both small AND fragmented cannot yield funannotate's
    # required 30 training models; predict would run for hours then abort with
    # "Not enough gene models N to train Augustus (30 required), exiting". Detect
    # that up front from cheap contig stats and skip cleanly (flag, no crash).
    # Requires BOTH gates so complete small genomes (e.g. Malassezia) are unaffected.
    # Disabled when predict_min_asm_bp=0.
    SKIP_REPORT="${params.target}/predict_skipped_too_small.tsv"
    if [ "${params.predict_min_asm_bp}" -gt 0 ]; then
        # Per-contig lengths -> sort descending -> N50 (portable; no gawk asort).
        read ASM_BP ASM_CTG ASM_N50 < <(
            awk '/^>/{if(len)print len;len=0;next}{len+=length(\$0)}END{if(len)print len}' "\$GENOME_IN" \\
            | sort -rn \\
            | awk '{L[NR]=\$1;tot+=\$1}END{half=tot/2;run=0;n50=0;for(i=1;i<=NR;i++){run+=L[i];if(run>=half){n50=L[i];break}}print tot, NR, n50}')
        echo "[INFO] Pre-flight assembly stats for ${out}: \${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}"
        SMALL=0; FRAG=0
        [ "\$ASM_BP" -lt "${params.predict_min_asm_bp}" ] && SMALL=1
        { [ "\$ASM_N50" -lt "${params.predict_frag_max_n50}" ] || [ "\$ASM_CTG" -gt "${params.predict_frag_max_contigs}" ]; } && FRAG=1
        if [ "\$SMALL" -eq 1 ] && [ "\$FRAG" -eq 1 ]; then
            echo "[WARN] ${out} is too small/fragmented for funannotate training (\${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_small_fragmented" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
    fi

    funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$PREDICTDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w codingquarry:0 glimmerhmm:0 \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table} --auto-skip-genemark || true

    # ── Post-predict catch ────────────────────────────────────────────────────
    # If predict produced no GBK, distinguish the known "too few training models"
    # outcome (an unfixable property of the assembly) from a genuine error. The
    # former is flagged and skipped so it does not abort the batch; anything else
    # still hard-fails so real problems surface.
    if [ ! -s "\$PREDICT_GBK" ]; then
        PLOG="\$PREDICTDIR/logfiles/funannotate-predict.log"
        if [ -f "\$PLOG" ] && grep -q "Not enough gene models .* to train Augustus" "\$PLOG"; then
            NMODELS=\$(grep -oE "Not enough gene models [0-9]+" "\$PLOG" | grep -oE "[0-9]+" | tail -1)
            echo "[WARN] ${out}: funannotate found only \${NMODELS:-<min} training models (needs 30); too small/fragmented to annotate — skipping" >&2
            mkdir -p "${params.target}"
            [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "funannotate_too_few_models:\${NMODELS:-NA}" "" "" "" >> "\$SKIP_REPORT"
            touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
            touch ${out}.predict.done
            exit 0
        fi
        echo "ERROR: funannotate predict did not produce expected GBK: \$PREDICT_GBK" >&2
        exit 1
    fi
    if [ -d "\$PREDICTDIR/predict_misc/ab_initio_parameters" ]; then
        mv "\$PREDICTDIR/predict_misc/ab_initio_parameters" "\$PREDICTDIR"
        mv "\$PREDICTDIR/predict_misc/trnascan.no-overlaps.gff3" "\$PREDICTDIR"
        rm -rf "\$PREDICTDIR/predict_misc"
        mkdir -p "\$PREDICTDIR/predict_misc"
        mv "\$PREDICTDIR/ab_initio_parameters" "\$PREDICTDIR/trnascan.no-overlaps.gff3" "\$PREDICTDIR/predict_misc"
    fi
    find "\$PREDICTDIR/predict_results/" -maxdepth 1 \\( -name "*.txt" -o -name "*.mrna-transcripts.fa" \\) -print0 \
        | xargs -0 --no-run-if-empty pigz
    sync
    touch ${out}.predict.done
    echo "[INFO] Prediction complete for ${out} at \$PREDICTDIR"
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] Would run funannotate predict for ${out} using ${genome_fa}"
    [ -f "${genome_fa}" ] || [ -f "${genome_fa}.gz" ] || { echo "ERROR: genome not found at ${genome_fa}[.gz]" >&2; exit 1; }
    mkdir -p ${params.target}/${out}/predict_results ${params.target}/${out}/predict_misc
    # non-empty so downstream size>0 gating (predict_ch / postpredict) is exercised
    echo "LOCUS stub_${out}" > ${params.target}/${out}/predict_results/${out}.gbk
    echo ">stub_${out}_p1" > ${params.target}/${out}/predict_results/${out}.proteins.fa
    touch ${out}.predict.done
    """
}

process FUNANNOTATE_ANNOTATE {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '32 GB'
    time   '48h'

    input:
    val(meta)

    output:
    tuple val(meta), path("${meta.id}.annotate.done"), emit: marker

    script:
    def out           = meta.id
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def antiSm    = file("${params.target}/${meta.id}/antismash_local/${meta.id}.gbk")
    def antiSmArg = antiSm.exists() ? "--antismash ${antiSm}" : ""
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}

    funannotate annotate -i ${params.target}/${out} -o ${params.target}/${out} \\
        --species "${species}" --strain "${strain}" \\
        --busco_db ${busco_lineage} --rename ${locustag} \\
        --sbt ${params.sbt_template} \\
        --header_length ${header_length} \\
        ${antiSmArg} \\
        --cpu ${task.cpus} --tmpdir \$TMPDIR

    EXPECTED_GBK="${params.target}/${out}/annotate_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate annotate did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    touch ${out}.annotate.done
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${params.target}/${out}/annotate_results ${params.target}/${out}/annotate_misc
    touch ${params.target}/${out}/annotate_results/${out}.gbk
    touch ${out}.annotate.done
    """
}

process FUNANNOTATE_UPDATE {
    label 'funannotate'
    tag "${meta.id}"

    cpus   16
    memory '96 GB'
    time   '48h'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    val meta

    script:
    def out           = meta.id
    def asmid         = meta.asmid
    def species       = meta.species
    def strain        = meta.strain
    def locustag      = meta.locustag
    def busco_lineage = meta.busco
    def header_length = params.header_length
    def pasa_db_arg = "--pasa_db sqlite"
    """
    # ── Skip if no reads (empty marker file from SRA_FETCH) ──────────────────
    if [ ! -s "${r1}" ]; then
        echo "[INFO] No RNAseq reads for ${out}, skipping funannotate update"
        exit 0
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    TMPDIR=\${SCRATCH:-/tmp}
    export PASACONF=""
    pasa_db_arg="--pasa_db sqlite"
    # ── Optional per-task MariaDB for PASA ────────────────────────────────────
    if [ "${params.pasa_mysql}" = "true" ]; then
        MYSQL_SCRATCH=${params.training_target}/${out}/training/mysql_db
        if [ ! -f \$MYSQL_SCRATCH/mysql/conf/my.cnf ]; then
            echo "[INFO] Setting up temporary MariaDB for PASA at \$MYSQL_SCRATCH"
            mkdir -p \$MYSQL_SCRATCH/db \$MYSQL_SCRATCH/conf
            rsync -a ${params.mysql_datadir}/mysql \$MYSQL_SCRATCH/db/ || \
                { echo "ERROR: Failed to copy mysql data from ${params.mysql_datadir}" >&2; exit 1; }
            cp ${params.pasa_conf_dir}/my.cnf \$MYSQL_SCRATCH/conf/my.cnf || \
                { echo "ERROR: Failed to copy my.cnf" >&2; exit 1; }
        fi
        MYHOSTNAME=\$(hostname -s)
        PORT=\$(shuf -i3000-4999 -n1)
        export PASACONF=\$MYSQL_SCRATCH/conf/pasa-local-\${MYHOSTNAME}.config.txt
        cp ${params.pasa_conf_dir}/conf.txt \$PASACONF
        sed -i "s/^MYSQLSERVER.*\$/MYSQLSERVER=\${MYHOSTNAME}:\${PORT}/" \$PASACONF
        perl -i -p -e "s/port = \\d+/port = \${PORT}/" \$MYSQL_SCRATCH/conf/my.cnf
        export SINGULARITY_BINDPATH=\$TMPDIR,\$MYSQL_SCRATCH/mysql_db
        stop_mysqldb() { singularity instance stop mysqldb_${asmid} 2>/dev/null || true; }
        trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM
        trap "stop_mysqldb" EXIT
        module load singularity
        singularity instance start --writable-tmpfs \\
            -B \$MYSQL_SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,\$MYSQL_SCRATCH/db/:/var/lib/mysql,\$MYSQL_SCRATCH/conf:/usr/conf \\
            ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
        pasa_db_arg="--pasa_db mysql"
        sleep 5
    fi

    # Link training data into work dir so funannotate update finds it at the relative path it expects.
    mkdir -p ${out}
    if [ -d "${params.training_target}/${out}/training" ]; then
        ln -sfn "${params.training_target}/${out}/training" "${out}/training"
    fi

    # r1/r2 are pre-normalized reads from SRA_FETCH (fastp-trimmed + bbnorm-normalized).
    # funannotate update will still run its internal alignment step against these.
    echo "[INFO] Running funannotate update for ${out}"
    funannotate update -i ${params.target}/${out} \\
        --left ${r1} --right ${r2} \\
        --cpus ${task.cpus} \\
        \$pasa_db_arg
    if [ "${params.pasa_mysql}" = "true" ]; then stop_mysqldb; fi
    echo "[INFO] stopped mysql"
    EXPECTED="${params.target}/${out}/update_results/${out}.gbk"
    if [ ! -f "\$EXPECTED" ]; then
        echo "ERROR: funannotate update did not produce expected GBK: \$EXPECTED" >&2
        exit 1
    fi
    """

    stub:
    def out = meta.id
    """
    echo "[STUB] FUNANNOTATE_UPDATE stub for ${out} (r1=${r1}, r2=${r2})"
    mkdir -p ${params.target}/${out}/update_results
    touch ${params.target}/${out}/update_results/${out}.tbl
    touch ${params.target}/${out}/update_results/${out}.gbk
    touch ${params.target}/${out}/update_results/${out}.gff3
    """
}

// A funannotate step's GenBank output may be stored uncompressed (.gbk) or
// gzip-compressed (.gbk.gz) so completed folders can be compressed to save space.
// Returns the existing non-empty file (preferring .gbk), or null if neither exists.
// Use this for completion/skip gating so a compressed result still counts as "done".
def gbkResult(String dir, String out) {
    def plain = file("${dir}/${out}.gbk")
    if (plain.exists() && plain.size() > 0) return plain
    def gz = file("${dir}/${out}.gbk.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return null
}

// Clean/masked genomes in input_clean_genomes may be stored gzip-compressed (.gz) to
// save space. Given the uncompressed base path (e.g. .../<asmid>.fa or
// .../<asmid>.masked.fasta), returns the existing non-empty file, preferring the
// compressed form. Falls back to the plain path object when neither exists, so callers'
// .exists() checks still report missing.
def genomeFile(String base) {
    def gz = file("${base}.gz")
    if (gz.exists() && gz.size() > 0) return gz
    return file(base)
}

def staleRnaseq(String out, String species) {
    def species_tag = species.replaceAll(/\s+/, '_')
    def gbk = gbkResult("${params.target}/${out}/predict_results", out)
    if (gbk == null) return false  // predict hasn't run yet; normal path handles it
    def r1      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_R1.fastq.gz")
    def se      = file("${launchDir}/rnaseq_reads/${species_tag}_norm_SE.fastq.gz")
    def trinity = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
    def r1_newer      = r1.exists()      && r1.size() > 0      && r1.lastModified()      > gbk.lastModified()
    def se_newer      = se.exists()      && se.size() > 0      && se.lastModified()      > gbk.lastModified()
    def trinity_newer = trinity.exists() && trinity.size() > 0 && trinity.lastModified() > gbk.lastModified()
    if (r1_newer || se_newer || trinity_newer) {
        log.info "stale prediction for ${out}: rnaseq/trinity newer than GBK — scheduling retrain+repredict"
        return true
    }
    return false
}

include { validateParameters; paramsSummaryLog; paramsHelp } from 'plugin/nf-schema'
include { ASM_STATS }        from './modules/local/asm_stats'
include { INPUT_CHECK }      from './subworkflows/local/input_check'
include { SETUP_DBS }        from './subworkflows/local/setup_dbs'
include { CLEAN_GENOMES }    from './subworkflows/local/clean_genomes'
include { MASK_GENOME }      from './subworkflows/local/mask_genome'
include { FETCH_RNASEQ }     from './subworkflows/local/fetch_rnaseq'
include { ANTISMASH_RUN }    from './modules/local/antismash_run'
include { INTERPROSCAN_RUN } from './modules/local/interproscan_run'
include { SIGNALP_RUN }      from './modules/local/signalp_run'

workflow {
    // `--help` prints schema-driven parameter help (grouped, with types/defaults) and exits.
    if (params.help) {
        log.info paramsHelp()
        exit 0
    }
    // Type-check params against nextflow_schema.json and log the resolved set.
    // (Unrecognised params warn rather than fail — see nextflow.config.)
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Fail fast with an actionable message when a pipeline profile was not selected
    // (these params come from conf/profile_annotate.config). Without it, downstream
    // file(params.funannotate_db) calls throw a cryptic "file() ... cannot be null".
    if( !params.taxondb || !params.funannotate_db )
        error "Missing params.taxondb / params.funannotate_db — add a pipeline profile, e.g. -profile annotate,slurm,module (or use: sbatch nextflow/run_annotate.sh)"

    // ── Samplesheet ingestion (INPUT_CHECK) ──────────────────────────────────
    // Parses samples CSV, applies taxon/asmid/suppress/n_test filters, builds
    // meta maps, and resolves genome paths. Two outputs:
    //   jobs        — tuple(meta, gz)  with genome existence filter (cleaning path)
    //   postpredict — meta only        no genome filter (annotate/update paths)
    INPUT_CHECK()
    def jobs = INPUT_CHECK.out.genomes

    def ch_versions = Channel.empty()
    if (params.debug.toBoolean()) {
        jobs.view { meta, gz -> "[CHANNEL] Submitting: out=${meta.id}, asmid=${meta.asmid}, transl_table=${meta.transl_table}, gz=${gz}" }
    }

    // Build/seed the three run-once databases. All use storeDir so they are no-ops
    // on any run where their target directories already exist.
    SETUP_DBS()
    def taxondb_ch = SETUP_DBS.out.taxondb

    CLEAN_GENOMES(jobs, taxondb_ch)

    if (!params.only_clean.toBoolean()) {
        def clean_genome_ch = CLEAN_GENOMES.out.genomes

        // ── Generate assembly statistics (for earlgrey_mask.nf SELECT_REPS) ────────
        // Generate asm_stats.tsv if --gen_asm_stats is true and the file doesn't exist.
        // This is used by earlgrey_mask.nf to select representative genomes per species.
        if (params.gen_asm_stats.toBoolean()) {
            def asm_stats_path = file(params.tables_dir).toAbsolutePath()
            def asm_stats_gz = file("${asm_stats_path}/asm_stats.tsv.gz")
            if (!asm_stats_gz.exists()) {
                log.info "Generating assembly statistics: ${asm_stats_gz}"
                ASM_STATS(
                    file(params.samples),
                    file(params.genome_dir)
                )
                ch_versions = ch_versions.mix(ASM_STATS.out.versions)
            } else {
                log.info "Assembly statistics already exist: ${asm_stats_gz}"
            }
        }

        // ── Repeat masking ────────────────────────────────────────────────────────
        // predict_genome_ch: tantan soft-masked (default) or clean/prior-masked genome.
        // MASK_GENOME handles the run_repeatmasker if/else and storeDir-cached masking.
        MASK_GENOME(clean_genome_ch)
        def predict_genome_ch = MASK_GENOME.out.genomes

        // Gate the predict chain on funannotate DB + augustus config being ready.
        // SETUP_DBS was already called above; its storeDir-cached outputs are free
        // on resumed runs. Gating here threads the dependency through the entire
        // downstream funannotate subgraph (train, predict, update, annotate).
        // (MASKREPEAT uses `funannotate mask`, which needs neither, so it is intentionally
        // left ungated and can run in parallel with these setup steps.)
        predict_genome_ch = predict_genome_ch
            .combine(SETUP_DBS.out.db)
            .combine(SETUP_DBS.out.config)
            .map { row -> row[0..-3] }

        // FUNANNOTATE_PREDICT input tuple drops taxonid (not needed after masking/clean).
        // When SRA is enabled: FETCH_RNASEQ fetches reads once per species; RNASEQ_PREPARE runs
        // funannotate train on the representative assembly and archives Trinity-GG transcripts
        // to rnaseq_data/; all other strains run FUNANNOTATE_TRAIN --trinity.
        def predict_input_ch
        def reads_ch = Channel.empty()
        if (params.run_sra_fetch.toBoolean()) {
            // Build per-species input: group assemblies, keep first taxonid per species.
            def sra_input = predict_genome_ch
                .map { meta, _genome_fa ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta.taxonid)
                }
                .groupTuple(by: 0)
                .map { species_tag, taxonids -> tuple(species_tag, taxonids[0]) }

            FETCH_RNASEQ(sra_input)
            reads_ch = FETCH_RNASEQ.out.reads

            if (!params.stop_after_sra_fetch.toBoolean() && !params.stop_after_sra_query.toBoolean()) {
            // Build per-assembly channel keyed by species_tag with SRA reads joined.
            // reads_ch is now a 4-tuple: (species_tag, r1, r2, se)
            def assembly_with_reads = predict_genome_ch
                .map { meta, genome_fa ->
                    def species_tag = meta.species.replaceAll(/\s+/, '_')
                    tuple(species_tag, meta, genome_fa)
                }
                .combine(reads_ch, by: 0)
            // assembly_with_reads: (species_tag, meta, genome_fa, r1, r2, se)

            // RNASEQ_PREPARE: run funannotate train --stop_after_trinity once per species on
            // the representative (first) assembly, then cache the Trinity-GG FASTA in rnaseq_data/
            // so all other strains share it. Normalized reads stay in rnaseq_reads/ (SRA_FETCH storeDir).
            // pasa.gff3 is NOT produced here (--stop_after_trinity stops before PASA);
            // it is produced by FUNANNOTATE_TRAIN for every strain including the representative.
            // Species whose representative r1 and se are both zero-length skip RNASEQ_PREPARE
            // entirely; an empty trinity FASTA is written locally without submitting a SLURM job.
            def repr_ch = assembly_with_reads
                .groupTuple(by: 0)
                .map { species_tag, metas, genomes, r1s, r2s, ses ->
                    tuple(species_tag, metas[0], genomes[0], r1s[0], r2s[0], ses[0])
                }

            def repr_branched = repr_ch.branch {
                has_reads: it[3].size() > 0 || it[5].size() > 0   // r1=[3] or se=[5]
                no_reads:  true
            }

            RNASEQ_PREPARE(repr_branched.has_reads)

            // For species with no RNA-seq reads, write an empty trinity FASTA to rnaseq_data/
            // in the driver process (no SLURM job) and emit it directly as a shared channel item.
            def empty_shared_ch = repr_branched.no_reads
                .map { species_tag, _meta, _gfa, _r1, _r2, _se ->
                    def empty_fa = file("${launchDir}/rnaseq_data/${species_tag}.trinity-GG.fasta")
                    if (!empty_fa.exists()) {
                        empty_fa.parent.mkdirs()
                        empty_fa.text = ''
                    }
                    tuple(species_tag, empty_fa)
                }

            def shared_ch = RNASEQ_PREPARE.out.shared.mix(empty_shared_ch)

            // Join shared Trinity from rnaseq_data back to every assembly for FUNANNOTATE_TRAIN.
            // Normalized reads (r1/r2/se) come from SRA_FETCH/SRA_FETCH_SE via assembly_with_reads.
            def train_input = assembly_with_reads
                .combine(shared_ch, by: 0)
                .map { species_tag, meta, genome_fa, r1, r2, se, trinity_fa ->
                    tuple(meta, genome_fa, r1, r2, se, trinity_fa)
                }
            // train_input: meta=0, genome_fa=1, r1=2, r2=3, se=4, trinity_fa=5

            // Branch on r1 (idx 2), se (idx 4), or trinity_fa (idx 5) sizes.
            // Assemblies with no RNA-seq bypass FUNANNOTATE_TRAIN entirely.
            def branched = train_input.branch {
                has_rnaseq: it[2].size() > 0 || it[4].size() > 0 || it[5].size() > 0
                no_rnaseq:  true
            }
            def predict_no_rnaseq = branched.no_rnaseq
                .map { meta, genome_fa, _r1, _r2, _se, _tf ->
                    tuple(meta, genome_fa)
                }

            // Skip TRAIN at the channel level when pasa.gff3 already exists and is non-empty,
            // UNLESS the rnaseq reads or trinity FASTA is newer than the existing prediction GBK
            // (staleRnaseq), in which case we re-run training so predict can be refreshed too.
            def train_todo = branched.has_rnaseq.filter { meta, _gfa, _r1, _r2, _se, _tf ->
                def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
                !gff3.exists() || gff3.size() == 0 || staleRnaseq(meta.id as String, meta.species as String)
            }
            def train_done = branched.has_rnaseq
                .filter { meta, _gfa, _r1, _r2, _se, _tf ->
                    def gff3 = file("${params.training_target}/${meta.id}/training/funannotate_train.pasa.gff3")
                    gff3.exists() && gff3.size() > 0 && !staleRnaseq(meta.id as String, meta.species as String)
                }
                .map { meta, genome_fa, _r1, _r2, _se, _tf ->
                    tuple(meta, genome_fa)
                }
            FUNANNOTATE_TRAIN(train_todo)
            predict_input_ch = FUNANNOTATE_TRAIN.out.mix(train_done).mix(predict_no_rnaseq)
            } // end if (!stop_after_sra_fetch && !stop_after_sra_query)
        } else {
            predict_input_ch = predict_genome_ch
        }

        if (!params.run_sra_fetch.toBoolean() || (!params.stop_after_sra_fetch.toBoolean() && !params.stop_after_sra_query.toBoolean())) {
        def predict_ch = predict_input_ch
            .filter { meta, _gfa ->
                gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) == null || staleRnaseq(meta.id as String, meta.species as String)
            }
        FUNANNOTATE_PREDICT(predict_ch)

        // ── Post-predict steps and annotation ────────────────────────────────────
        // postpredict: all samples with a completed predict_results/*.gbk, whether
        // produced in this run or a prior one. This is the source for all optional
        // pre-annotate steps and for FUNANNOTATE_ANNOTATE itself.
        // INPUT_CHECK.out.samples already has taxon/asmid/suppress/n_test filters applied;
        // we just add the predict-results existence check on top.
        def postpredict = INPUT_CHECK.out.samples
            // Only genomes whose prediction was already complete AND current in a PRIOR run.
            // This is the exact logical complement of the predict_ch filter, so this set is
            // disjoint from the genomes (re)predicted in THIS run (which arrive via
            // FUNANNOTATE_PREDICT.out.metadata below). Keeping them disjoint means no genome
            // is fed downstream twice and stale genomes correctly wait for the fresh predict.
            .filter { meta ->
                gbkResult("${params.target}/${meta.id}/predict_results", meta.id as String) != null && !staleRnaseq(meta.id as String, meta.species as String)
            }

        // annotate_ready_ch threads through optional pre-annotate steps. Each optional
        // step splits the channel into "needs to run" vs "already done", processes the
        // former, then mixes the freshly-completed items back. FUNANNOTATE_ANNOTATE only
        // fires once all requested optional steps are complete for a given sample.
        // Joining ANTISMASH/INTERPRO/SIGNALP output back through predict_meta reconstructs
        // the metadata tuple while encoding the dependency edge in the channel DAG.
        //
        // Same-run completion gate: genomes predicted in THIS run flow in via
        // FUNANNOTATE_PREDICT.out.metadata (a real channel edge, so downstream waits for
        // predict to finish), while prior-run genomes flow in via postpredict (available
        // immediately). The two sets are disjoint by the filters above, so a plain mix
        // needs no dedup. (The optional steps below are still each gated behind their
        // params — run_antismash/interpro/signalp/update/annotate, all default false.)
        def predict_meta = postpredict.mix(FUNANNOTATE_PREDICT.out.metadata)
        def annotate_ready_ch = predict_meta

        if (params.run_antismash.toBoolean()) {
            def as_todo = annotate_ready_ch.filter { meta ->
                def asDir = file("${params.target}/${meta.id}/antismash_local")
                !(asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') })
            }
            def as_done = annotate_ready_ch.filter { meta ->
                def asDir = file("${params.target}/${meta.id}/antismash_local")
                asDir.isDirectory() && asDir.list()?.any { it.endsWith('.json') || it.endsWith('.json.gz') }
            }
            ANTISMASH_RUN(as_todo)
            ch_versions = ch_versions.mix(ANTISMASH_RUN.out.versions)
            def as_completed = ANTISMASH_RUN.out.results
                .map { meta, _files -> meta }
            annotate_ready_ch = as_completed.mix(as_done)
        }

        if (params.run_interpro.toBoolean()) {
            def ipr_todo = annotate_ready_ch.filter { meta ->
                !file("${params.target}/${meta.id}/annotate_misc/iprscan.xml").exists()
            }
            def ipr_done = annotate_ready_ch.filter { meta ->
                file("${params.target}/${meta.id}/annotate_misc/iprscan.xml").exists()
            }
            INTERPROSCAN_RUN(ipr_todo)
            ch_versions = ch_versions.mix(INTERPROSCAN_RUN.out.versions)
            def ipr_completed = INTERPROSCAN_RUN.out.results
                .map { meta, _xml -> meta }
            annotate_ready_ch = ipr_completed.mix(ipr_done)
        }

        if (params.run_signalp.toBoolean()) {
            def sp_todo = annotate_ready_ch.filter { meta ->
                !file("${params.target}/${meta.id}/annotate_misc/signalp.results.txt").exists()
            }
            def sp_done = annotate_ready_ch.filter { meta ->
                file("${params.target}/${meta.id}/annotate_misc/signalp.results.txt").exists()
            }
            SIGNALP_RUN(sp_todo)
            ch_versions = ch_versions.mix(SIGNALP_RUN.out.versions)
            def sp_completed = SIGNALP_RUN.out.results
                .map { meta, _txt -> meta }
            annotate_ready_ch = sp_completed.mix(sp_done)
        }

        if (params.run_update.toBoolean()) {
            if (params.run_sra_fetch.toBoolean()) {
                // UPDATE runs from predict results in parallel with antismash/interpro/signalp.
                // Reads are joined from SRA_FETCH (storeDir-cached, so prior-run reads are reused).
                // The join on upd_signal gates annotate_ready_ch so ANNOTATE waits for UPDATE.
                def upd_input = predict_meta
                    .map { meta ->
                        def species_tag = meta.species.replaceAll(/\s+/, '_')
                        tuple(species_tag, meta)
                    }
                    .combine(reads_ch, by: 0)
                    .map { _st, meta, r1, r2 ->
                        tuple(meta, r1, r2)
                    }
                def upd_todo = upd_input.filter { meta, _r1, _r2 ->
                    gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) == null
                }
                def upd_done_signal = upd_input
                    .filter { meta, _r1, _r2 ->
                        gbkResult("${params.target}/${meta.id}/update_results", meta.id as String) != null
                    }
                    .map { meta, _r1, _r2 -> tuple(meta.id, 'upd') }
                FUNANNOTATE_UPDATE(upd_todo)
                def upd_signal = FUNANNOTATE_UPDATE.out
                    .map { meta -> tuple(meta.id, 'upd') }
                    .mix(upd_done_signal)
                annotate_ready_ch = annotate_ready_ch
                    .map { meta -> tuple(meta.id, meta) }
                    .join(upd_signal)
                    .map { _id, meta, _flag -> meta }
            } else {
                log.warn "run_update=true but run_sra_fetch=false; funannotate update skipped (no reads available)"
            }
        }

        if (params.run_annotate.toBoolean()) {
            FUNANNOTATE_ANNOTATE(annotate_ready_ch.filter { meta ->
                gbkResult("${params.target}/${meta.id}/annotate_results", meta.id as String) == null
            })
        }
        } // end if (!params.stop_after_sra_fetch || !params.run_sra_fetch)
    }

    // Collect software versions from all processes that emit versions.yml.
    // Written to logs/software_versions.yml alongside the trace file.
    ch_versions
        .unique()
        .collectFile(
            name:     'software_versions.yml',
            storeDir: "${launchDir}/logs/nextflow",
            newLine:  true
        )
}

