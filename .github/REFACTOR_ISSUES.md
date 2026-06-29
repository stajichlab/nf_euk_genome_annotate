# Refactor issue tracker

Source-of-truth backlog for the DSL2 modularization (see `REFACTORING_PLAN.md`).
Create these on GitHub with `scripts/create_refactor_issues.sh` (uses `gh`).
Each non-epic issue closes one `## Phase N` of the plan and must pass the
stub-run gate before merge:

```
nextflow config . -profile test && nextflow run . -profile test -stub-run && nextflow run . --help
```

---

## EPIC: Modularize nf_funannotate1 into DSL2 modules + subworkflows
labels: refactor, epic

Break the 2,470-line `funannotate.nf` monolith into one-process-per-file modules
composed by subworkflows, on a `meta`-map data contract. nf-core-*inspired*, not
nf-core-submitted (see `REFACTORING_PLAN.md` for the verdict). Child issues #1–#12
below; they are ordered — #1 and #2 block the rest.

Definition of done: monolith replaced by `main.nf` + `workflows/` +
`subworkflows/local/` + `modules/local/`; every process emits `versions.yml`;
earlgrey/funannotate share `INPUT_CHECK`; CI green at every step.

---

## 1. Adopt `meta`-map data contract  (Phase 0 — BLOCKER)
labels: refactor

Replace the positional 10-tuple
`tuple(out, asmid, species, strain, locustag, busco, header_length, transl_table, gz, taxonid)`
with `tuple val(meta), path(genome)` where `meta` is the map defined in
`REFACTORING_PLAN.md` (Principle 0). `meta.id` is the only naming field;
`header_length` becomes `params.header_length`. No more extraction until this lands.

- [ ] Build `meta` in the workflow channel construction
- [ ] Update every call site to consume `meta`
- [ ] `params.header_length` added to schema with default 24
- [ ] Stub-run gate green

## 2. Shared `INPUT_CHECK` subworkflow (dedupe earlgrey)  (Phase 0)
labels: refactor

Extract the samplesheet `splitCsv` + taxon/asmid/suppress filtering (duplicated
in `funannotate.nf` and `earlgrey_mask.nf`) into `subworkflows/local/input_check.nf`
emitting the `meta` channel. Both entrypoints call it.

- [ ] `subworkflows/local/input_check.nf` emits `ch_genomes` (meta + genome)
- [ ] `funannotate.nf` and `earlgrey_mask.nf` both consume it; no duplicated parse
- [ ] Stub-run gate green for both entrypoints

## 3. Repo skeleton + relocate existing modules  (Phase 1)
labels: refactor

Create `main.nf`, `workflows/funannotate.nf`, `subworkflows/local/`,
`modules/local/`. Move `modules/asm_stats.nf` and `modules/annotation_tools.nf`
to one-process-per-file under `modules/local/` (split `annotation_tools.nf` into
`antismash.nf` / `signalp.nf` / `interproscan.nf`).

- [ ] Directory skeleton in place; `main.nf` is a thin entrypoint
- [ ] `annotation_tools.nf` split into 3 single-process modules
- [ ] Stub-run gate green

## 4. `versions.yml` + `conf/base.config` + `conf/modules.config`  (Phase 2)
labels: refactor

Establish the conventions every later module copies: each process emits
`versions.yml`; resources move to label-based `conf/base.config`
(`process_low/medium/high`); per-process `publishDir`/`ext.args` move to
`conf/modules.config`.

- [ ] `conf/base.config` with resource labels
- [ ] `conf/modules.config` with publishDir + ext.args
- [ ] At least one module emits and the workflow collects `versions.yml`
- [ ] Stub-run gate green

## 5. Setup modules  (Phase 3)
labels: refactor, good first issue

Extract `SETUP_TAXONDB`, `SETUP_FUNANNOTATE_DB`, `SETUP_AUGUSTUS_CONFIG` into
`modules/local/` + `subworkflows/local/setup_dbs.nf`. Preserve `storeDir`
(run-at-most-once) caching. Good first real extraction — validates the gate.

- [ ] 3 modules + `setup_dbs.nf` subworkflow; storeDir preserved
- [ ] Stub-run gate green

## 6. Genome clean + `prepare_genome` subworkflow  (Phase 4)
labels: refactor

Extract `GENOME_CLEAN` / `GENOME_CLEAN_BATCH` into modules and a
`subworkflows/local/prepare_genome.nf` (clean → asm_stats → mask).
**Preserve the FCS-GX `/dev/shm` staging** and the "skip already-cleaned" batch
gating so a fully-cleaned batch never pays the ~30-min staging cost.

- [ ] Modules + `prepare_genome.nf`; FCS-GX /dev/shm staging preserved
- [ ] Batch padding/skip behavior unchanged
- [ ] Stub-run gate green

## 7. Masking subworkflow + per-tool modules  (Phase 5)
labels: refactor

Replace the planned single masking mega-process. One module per masker
(`mask_tantan.nf` now; `mask_repeatmodeler.nf`, `mask_repeatmasker.nf`,
`mask_earlgrey.nf` as stubs/follow-ups); selection logic in
`subworkflows/local/mask.nf` keyed on `params.mask_tool`.

- [ ] `mask.nf` selects one masker by param; `NONE` path supported
- [ ] `mask_tantan.nf` extracted; others stubbed with clear TODO
- [ ] Stub-run gate green

## 8. RNA-seq fetch subworkflow  (Phase 6 — hardest)
labels: refactor

Extract `SRA_QUERY` / `SRA_QUERY_BATCH` / `COLLECT_SRA_QUERY` /
`WRITE_EMPTY_READS` / `SRA_FETCH` / `SRA_FETCH_SE` / `RNASEQ_PREPARE` into modules
+ `subworkflows/local/rnaseq.nf`. Done **after** the pattern is proven on easier
processes. Preserve per-species shared Trinity-GG output and `maxForks` limits.

- [ ] 7 modules + `rnaseq.nf`; shared Trinity-GG semantics preserved
- [ ] maxForks / rate limits preserved
- [ ] Stub-run gate green

## 9. Funannotate predict subworkflow  (Phase 7)
labels: refactor

Extract `FUNANNOTATE_TRAIN` / `FUNANNOTATE_PREDICT` / `FUNANNOTATE_UPDATE` into
modules + `subworkflows/local/predict.nf`. Keep the pre-flight assembly
size/fragmentation validation and post-flight "not enough models" guard.

- [ ] 3 modules + `predict.nf`; pre/post-flight checks preserved
- [ ] `update` is optional (param-gated)
- [ ] Stub-run gate green

## 10. Annotation subworkflow  (Phase 8)
labels: refactor

`subworkflows/local/annotate.nf` composing optional `antismash` / `signalp` /
`interproscan` modules → `funannotate_annotate.nf`. Each optional tool
independently param-gated.

- [ ] `annotate.nf` with per-tool gating; merges results into funannotate annotate
- [ ] Stub-run gate green

## 11. Consolidate `ucr_hpcc` institutional profile + portable container path  (Phase 9)
labels: refactor

The `module`→`ucr_hpcc` rename is done. Finish the repivot: fold UCR SLURM
partitions / `clusterOptions` into the institutional profile, and ensure a fully
portable run works via per-module conda/biocontainer directives (no Lmod).
Document the "copy to `conf/provision_<site>.config`" path for new sites.

- [ ] UCR partition config consolidated under the institutional profile
- [ ] Portable container/conda path runs without any UCR modules
- [ ] `docs/` note for adding a new institution

## 12. nf-core hygiene (stretch)  (Phase 9)
labels: refactor, documentation

`docs/usage.md` + `docs/output.md`, `assets/schema_input.json` (samplesheet
schema), pipeline naming decision (nf-core forbids underscores/digits), optional
MultiQC + `nf-test`. Decide explicitly whether to pursue nf-core submission or
stay nf-core-inspired.

- [ ] `docs/usage.md`, `docs/output.md`
- [ ] `assets/schema_input.json`
- [ ] naming + submission decision recorded in `REFACTORING_PLAN.md`
