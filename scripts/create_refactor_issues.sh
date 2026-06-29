#!/usr/bin/env bash
#
# create_refactor_issues.sh — create the modularization backlog on GitHub.
#
# Idempotent-ish: skips any issue whose exact title already exists (open OR
# closed). Creates the `refactor`/`epic` labels and a milestone if missing.
# Requires: gh (authenticated), run from inside the repo.
#
#   ./scripts/create_refactor_issues.sh            # create everything
#   DRY_RUN=1 ./scripts/create_refactor_issues.sh  # print what would happen
#
# Mirrors .github/REFACTOR_ISSUES.md — keep them in sync.

set -euo pipefail

REPO="$(gh repo view --json nameWithOwner -q .nameWithOwner)"
MILESTONE="DSL2 modularization"
DRY_RUN="${DRY_RUN:-0}"

say()  { printf '%s\n' "$*"; }
run()  { if [ "$DRY_RUN" = 1 ]; then say "DRY: $*"; else "$@"; fi; }

ensure_label() {  # name color description
    if ! gh label list --repo "$REPO" --limit 200 | grep -qiE "^$1[[:space:]]"; then
        run gh label create "$1" --repo "$REPO" --color "$2" --description "$3" || true
    fi
}

ensure_milestone() {
    if ! gh api "repos/$REPO/milestones?state=all" -q '.[].title' | grep -qxF "$MILESTONE"; then
        run gh api "repos/$REPO/milestones" -f title="$MILESTONE" \
            -f description="Break funannotate.nf monolith into modules + subworkflows" >/dev/null || true
    fi
}

issue_exists() {  # title
    gh issue list --repo "$REPO" --state all --limit 500 --json title -q '.[].title' \
        | grep -qxF "$1"
}

make_issue() {  # title labels body
    local title="$1" labels="$2" body="$3"
    if issue_exists "$title"; then
        say "skip (exists): $title"
        return
    fi
    run gh issue create --repo "$REPO" --title "$title" --label "$labels" \
        --milestone "$MILESTONE" --body "$body"
}

ensure_label refactor 1d76db "DSL2 modularization work"
ensure_label epic     5319e7 "Tracking epic"
ensure_milestone

make_issue "EPIC: Modularize nf_funannotate1 into DSL2 modules + subworkflows" "refactor,epic" \
"Break the 2,470-line funannotate.nf monolith into one-process-per-file modules composed by subworkflows, on a meta-map data contract. See REFACTORING_PLAN.md and .github/REFACTOR_ISSUES.md. Child issues are ordered; the meta-map and INPUT_CHECK issues block the rest."

make_issue "Adopt meta-map data contract (Phase 0, BLOCKER)" "refactor" \
"Replace the positional 10-tuple with \`tuple val(meta), path(genome)\`. meta.id is the only naming field; header_length becomes params.header_length. No further extraction until this lands. See REFACTORING_PLAN.md Principle 0."

make_issue "Shared INPUT_CHECK subworkflow (dedupe earlgrey)" "refactor" \
"Extract samplesheet parse + taxon/asmid/suppress filtering (duplicated in funannotate.nf and earlgrey_mask.nf) into subworkflows/local/input_check.nf emitting the meta channel. Both entrypoints consume it."

make_issue "Repo skeleton + relocate existing modules (Phase 1)" "refactor" \
"Create main.nf, workflows/, subworkflows/local/, modules/local/. Move asm_stats.nf and split annotation_tools.nf into antismash/signalp/interproscan single-process modules."

make_issue "versions.yml + conf/base.config + conf/modules.config (Phase 2)" "refactor" \
"Each process emits versions.yml; resources move to label-based conf/base.config; per-process publishDir/ext.args move to conf/modules.config. The pattern every later module copies."

make_issue "Extract setup modules (Phase 3)" "refactor,good first issue" \
"Extract SETUP_TAXONDB / SETUP_FUNANNOTATE_DB / SETUP_AUGUSTUS_CONFIG into modules/local + subworkflows/local/setup_dbs.nf. Preserve storeDir run-at-most-once caching."

make_issue "Genome clean + prepare_genome subworkflow (Phase 4)" "refactor" \
"Extract GENOME_CLEAN / GENOME_CLEAN_BATCH + subworkflows/local/prepare_genome.nf (clean -> asm_stats -> mask). Preserve FCS-GX /dev/shm staging and the skip-already-cleaned batch gating."

make_issue "Masking subworkflow + per-tool modules (Phase 5)" "refactor" \
"Replace the single masking mega-process: one module per masker (tantan now; repeatmodeler/repeatmasker/earlgrey follow-ups), selection in subworkflows/local/mask.nf keyed on params.mask_tool. Support NONE."

make_issue "RNA-seq fetch subworkflow (Phase 6, hardest)" "refactor" \
"Extract SRA_QUERY/_BATCH, COLLECT_SRA_QUERY, WRITE_EMPTY_READS, SRA_FETCH/_SE, RNASEQ_PREPARE into modules + subworkflows/local/rnaseq.nf. Done after pattern proven. Preserve shared Trinity-GG and maxForks."

make_issue "Funannotate predict subworkflow (Phase 7)" "refactor" \
"Extract FUNANNOTATE_TRAIN/PREDICT/UPDATE into modules + subworkflows/local/predict.nf. Keep pre-flight size/fragmentation checks and post-flight not-enough-models guard. update is optional/param-gated."

make_issue "Annotation subworkflow (Phase 8)" "refactor" \
"subworkflows/local/annotate.nf composing optional antismash/signalp/interproscan modules -> funannotate_annotate.nf, each independently param-gated."

make_issue "Consolidate ucr_hpcc institutional profile + portable container path (Phase 9)" "refactor" \
"module->ucr_hpcc rename is done. Fold UCR SLURM partitions/clusterOptions into the institutional profile; ensure a fully portable container/conda run works without Lmod. Document the new-site path."

make_issue "nf-core hygiene (Phase 9, stretch)" "refactor,documentation" \
"docs/usage.md + docs/output.md, assets/schema_input.json, naming decision (nf-core forbids underscores/digits), optional MultiQC + nf-test. Record nf-core submission vs nf-core-inspired decision in REFACTORING_PLAN.md."

say "Done."
