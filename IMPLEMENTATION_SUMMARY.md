# nf_funannotate1 Modularization: Implementation Summary & Review Request

**Date:** 2026-06-28  
**Status:** Plan Complete, Ready for Expert Review  
**Branch:** main  

---

## Executive Summary

A comprehensive modularization plan has been developed to refactor the monolithic `funannotate.nf` pipeline into reusable, independently-testable modules. The plan organizes processes into 4 phases reflecting the actual workflow structure:

1. **Genome Preprocessing** (AAFTF + Repeat Masking)
2. **Gene Prediction** (RNA-seq + Funannotate)
3. **Annotation** (Post-prediction tools)
4. **Setup & Utilities** (Database initialization)

Three working implementations are production-ready, with a detailed rollout plan for the remaining components.

---

## What's Complete ✅

### 1. ASM_STATS Module
- **File:** `modules/asm_stats.nf`
- **Status:** ✅ Working, production-ready
- **Function:** Generate assembly statistics (ASMID, total_length_bp, N50_bp, contig_count)
- **Note:** Will be relocated to `modules/AAFTF/asm_stats.nf` in Phase 1
- **Integration:** Used by both `funannotate.nf` and `earlgrey_mask.nf`

### 2. Optional SELECT_REPS
- **File:** `earlgrey_mask.nf` (updated)
- **Status:** ✅ Working, production-ready
- **Feature:** `--skip_select_reps` flag
- **Function:** Process all genomes for EarlGrey without size filtering
- **Parameter:** `--gen_asm_stats` (default: true) auto-generates assembly stats

### 3. Annotation Tools Module
- **File:** `modules/annotation_tools.nf`
- **Status:** ✅ Working, production-ready
- **Processes:**
  - ANTISMASH_RUN (secondary metabolite detection)
  - INTERPROSCAN_RUN (protein domain annotation)
  - SIGNALP_RUN (signal peptide prediction)
- **Note:** Will be relocated to `modules/annotate/annotation_tools.nf` in Phase 3

### 4. Documentation
- **REFACTORING_PLAN.md** - Complete phase breakdown with workflow diagram
- **modules/README.md** - Development guidelines, usage examples, testing procedures
- **MODULE_STRUCTURE.txt** - ASCII diagram showing module hierarchy and data flow

---

## Proposed Modularization Plan

### Phase 1: Genome Preprocessing (AAFTF & Repeat Masking)

**Directory Structure:**
```
modules/
├── AAFTF/
│   ├── asm_stats.nf (move from root)
│   ├── FCS_GX.nf (extract from GENOME_CLEAN)
│   ├── sourpurge.nf (extract from GENOME_CLEAN)
│   └── vecscreen.nf (new)
└── repeatmasking/
    └── masking.nf (strategy selector: TANTAN, REPEATMODELER, REPEATMASKER, EARLGREY, NONE)
```

**Key Decisions:**
- Consolidate all genome cleaning tools in modules/AAFTF/
- Create flexible masking strategy selector to support multiple methods
- Preserve FCS-GX /dev/shm caching pattern during extraction

---

### Phase 2: Gene Prediction (RNA-seq & Funannotate)

**Directory Structure:**
```
modules/
├── rnaseq_fetch/
│   ├── sra_query.nf (SRA_QUERY, SRA_QUERY_BATCH, COLLECT_SRA_QUERY)
│   ├── sra_fetch.nf (SRA_FETCH, SRA_FETCH_SE, WRITE_EMPTY_READS)
│   └── prepare.nf (RNASEQ_PREPARE)
└── funannotate/
    ├── train.nf (FUNANNOTATE_TRAIN)
    ├── predict.nf (FUNANNOTATE_PREDICT)
    └── update.nf (FUNANNOTATE_UPDATE - optional)
```

**Key Decisions:**
- Separate RNA-seq tools into dedicated module for reusability
- Keep funannotate core (train/predict/update) together
- Preserve complex channel workflows and skip-caching patterns

---

### Phase 3: Annotation Workflow

**Directory Structure:**
```
modules/annotate/
├── annotation_tools.nf (move from root - ANTISMASH, SIGNALP, INTERPROSCAN)
└── funannotate.nf (FUNANNOTATE_ANNOTATE)
```

**Key Decisions:**
- Group annotation tools together for easy composition
- Enable independent tool selection (run any/all of them)
- Support conditional execution per tool

---

### Phase 4: Setup & Utilities

**Directory Structure:**
```
modules/setup/
└── databases.nf (SETUP_TAXONDB, SETUP_FUNANNOTATE_DB, SETUP_AUGUSTUS_CONFIG)
```

**Key Decisions:**
- Foundational utilities used across all phases
- Preserve storeDir caching behavior (run at most once)
- Maintain current parameter organization

---

## Implementation Priority

1. **Phase 2.1 (RNA-seq Fetch)** ← START HERE
   - Least interdependent
   - High reusability potential
   - Enables flexible training workflows

2. **Phase 2.2 (Funannotate Core)**
   - Enables flexible prediction pipelines
   - Supports optional update step
   - Well-tested existing code

3. **Phase 1 (Genome Preprocessing)**
   - More complex conditional branching
   - Foundation for rest of pipeline
   - Larger refactoring effort

4. **Phase 3 (Annotation)**
   - Depends on Phase 2
   - Enables composition of optional tools
   - Most flexible phase

5. **Phase 4 (Setup)**
   - Last; foundational only
   - No blocking dependencies

---

## Benefits & Use Cases

### Immediate Benefits
| Benefit | Details |
|---------|---------|
| **Reusability** | Modules can be used in different pipelines (e.g., funannotate-only vs. funannotate+earlgrey) |
| **Flexibility** | Easy to swap masking strategies or annotation tools without modifying main workflow |
| **Testability** | Each module testable in isolation with `-stub-run` |
| **Maintainability** | Smaller, focused files easier to understand and modify |
| **Documentation** | Clear module contracts (inputs/outputs/parameters) |

### Future Use Cases
- Simplified parameter composition for specific workflows
- Easier integration of new tools (e.g., new masking strategies)
- Support for different organism pipelines (fungi → eukaryotes → prokaryotes)
- CI/CD pipeline testing per module
- Community contributions of alternative modules

---

## Documentation Files

All files are committed to git and available in the repository:

1. **REFACTORING_PLAN.md** (177 lines)
   - Detailed phase breakdown
   - Complete workflow diagram
   - Module-by-module responsibility assignment
   - Rollout plan with priority order

2. **modules/README.md** (176 lines)
   - Development guidelines for module contributors
   - Module documentation template
   - Testing procedures (stub-run, unit tests)
   - Parameter handling best practices
   - Implementation checklist

3. **MODULE_STRUCTURE.txt** (280 lines)
   - ASCII diagram of module hierarchy
   - Data flow visualization for each phase
   - Input/output contracts per process
   - Current vs. planned status matrix

4. **IMPLEMENTATION_SUMMARY.md** (this file)
   - Executive summary
   - Complete artifacts list
   - Review focus areas
   - Specific recommendations for reviewers

---

## Items Requiring Expert Review

### 1. Architectural Soundness
- [ ] Are phase dependencies correctly ordered?
- [ ] Are module boundaries well-defined?
- [ ] Will module reusability work as designed?
- [ ] Any missing interdependencies?

### 2. Implementation Feasibility
- [ ] Are extraction patterns consistent?
- [ ] Will existing parameters need major reorganization?
- [ ] Hidden interdependencies violating modularity?
- [ ] Can each phase be developed independently?

### 3. Performance Impact
- [ ] Will include statements add overhead?
- [ ] Are storeDir/publishDir patterns preserved?
- [ ] Resume/checkpoint functionality unaffected?
- [ ] Any workflow DAG complexity changes?

### 4. Backward Compatibility
- [ ] Can existing wrapper scripts continue working?
- [ ] Will parameter changes be transparent?
- [ ] How to handle user-provided nextflow.config?
- [ ] Deprecation strategy for old parameter names?

### 5. Testing Strategy
- [ ] Are stub-run tests sufficient for module validation?
- [ ] Should we add integration tests?
- [ ] How to validate module independence?
- [ ] CI/CD pipeline changes needed?

### 6. Documentation Adequacy
- [ ] Are module templates clear and complete?
- [ ] Need more implementation examples?
- [ ] Data contracts documented sufficiently?
- [ ] Guidelines for when/where to create new modules?

---

## Specific Items for Review

### Critical Decision Points

1. **RNA-seq Fetch Extraction (Phase 2.1)**
   - Complex interdependencies: SRA_QUERY → SRA_FETCH → RNASEQ_PREPARE
   - Are channel flows properly documented?
   - Can parameters be cleanly separated per module?
   - **Question:** Should this be extracted in smaller sub-phases?

2. **Repeat Masking Strategy Selector (Phase 1)**
   - Multiple conditional branches (NONE, TANTAN, REPEATMODELER, REPEATMASKER, EARLGREY)
   - How should strategy selection work?
   - Should each be a separate module or combined?
   - **Question:** Is `modules/repeatmasking/masking.nf` the right structure?

3. **AAFTF Contamination Tools (Phase 1)**
   - FCS_GX currently in GENOME_CLEAN_BATCH with complex /dev/shm staging
   - Can this be cleanly extracted without breaking performance?
   - How to preserve one-time 30-min staging cost amortization?
   - **Question:** Should FCS_GX remain batched or become per-genome?

4. **Setup Database Ordering (Phase 4)**
   - Current: main workflow calls SETUP_*, then branches
   - Can setup be fully parallelized from main?
   - Do modules need to call setup themselves or assume it's done?
   - **Question:** Should setup be implicit or explicit in each module?

---

## Current Git History

```
80d05fe Add comprehensive module structure documentation
9867d96 Update modularization plan with detailed workflow structure
ec925c5 Add annotation_tools module and refactoring plan
4c50fd8 Add modular ASM_STATS and make SELECT_REPS optional
b13cd84 Merge pull request #1 from stajichlab/copilot/fix-parse-config-stub-run
```

---

## Next Steps (Pending Approval)

### Immediate (This Week)
- [ ] Schedule expert code review
- [ ] Collect feedback on critical decision points
- [ ] Refine Phase 2.1 extraction strategy based on review

### Short-term (This Sprint)
- [ ] Start Phase 2.1 implementation (RNA-seq modules)
- [ ] Create feature branch for Phase 2.1
- [ ] Extract and test individual modules

### Medium-term (Next 2-3 Sprints)
- [ ] Complete remaining phases
- [ ] Update main funannotate.nf to include modules progressively
- [ ] Document module usage patterns

### Long-term (Post-implementation)
- [ ] Update CI/CD to test modules independently
- [ ] Create community contribution guidelines for new modules
- [ ] Document use cases for different organism types

---

## Appendix: Files Ready for Review

All documentation is production-ready and committed to the main branch:

1. `REFACTORING_PLAN.md` - Complete architecture plan
2. `modules/README.md` - Developer guidelines
3. `MODULE_STRUCTURE.txt` - Visual reference
4. `modules/asm_stats.nf` - Working implementation
5. `modules/annotation_tools.nf` - Working implementation
6. `earlgrey_mask.nf` - Updated with optional SELECT_REPS
7. `conf/profile_annotate.config` - Updated parameters

---

**Prepared by:** Claude Code (AI)  
**Date:** 2026-06-28  
**Review Status:** ⏳ Awaiting Expert Review
