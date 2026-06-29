/*
 * INPUT_CHECK — canonical samplesheet ingestion
 *
 * Parses params.samples (CSV with header), applies taxon / ASMID / suppress /
 * n_test filters, builds SampleUtils.makeMeta(row), and resolves the genome
 * file from the GENOME column or the NCBI_ASM source directory.
 *
 * Emits two channels so callers can choose the right level of filtering:
 *
 *   samples  — val(meta)              filtered by taxon/asmid/suppress/n_test
 *                                     but NOT by genome existence (use for any
 *                                     workflow that just needs sample metadata,
 *                                     e.g. postpredict or earlgrey select path)
 *
 *   genomes  — tuple(val(meta), path(gz))
 *                                     additionally filtered to samples whose
 *                                     raw genome file exists on disk (use as
 *                                     the primary jobs channel in funannotate)
 *
 * All filtering logic lives here so it does not have to be duplicated across
 * funannotate.nf, earlgrey_mask.nf, or any future entrypoint.
 */

workflow INPUT_CHECK {

    main:

    // ── Suppress list ────────────────────────────────────────────────────────
    def suppressSet = (params.suppress && file(params.suppress, glob: false).exists())
        ? file(params.suppress, glob: false).readLines()
              .collect { it.trim().split(',')[0].trim() }
              .findAll { it && !it.startsWith('#') }
              .toSet()
        : ([] as Set)
    if (suppressSet) {
        log.info "Suppress list loaded: ${suppressSet.size()} ASMIDs will be skipped"
    }

    // ── Taxonomy filter ──────────────────────────────────────────────────────
    // Parse --taxon RANK:VALUE (e.g. --taxon PHYLUM:Ascomycota).
    def taxonFilter
    if (params.taxon) {
        def parts = (params.taxon as String).split(':', 2)
        if (parts.size() != 2 || !parts[0] || !parts[1]) {
            error "--taxon must be in RANK:VALUE format, e.g. --taxon PHYLUM:Ascomycota"
        }
        def taxRank  = parts[0].toUpperCase()
        def taxValue = parts[1]
        log.info "Taxonomy filter: ${taxRank} = '${taxValue}'"
        taxonFilter = { row -> row[taxRank]?.trim() == taxValue }
    } else {
        taxonFilter = { row -> true }
    }

    // ── ASMID filter ─────────────────────────────────────────────────────────
    def asmidFilter = params.asmid
        ? { row -> row.ASMID?.trim() == (params.asmid as String).trim() }
        : { row -> true }
    if (params.asmid) {
        log.info "ASMID filter: processing only '${params.asmid}'"
    }

    // ── Build channel with meta + resolved genome path ───────────────────────
    // ch_raw passes taxon/asmid/suppress/n_test filters; gz resolution is
    // done here so both emit channels share the same genome-path logic.
    def ch_raw = channel.fromPath(params.samples)
        .splitCsv(header: true)
        .filter(taxonFilter)
        .filter(asmidFilter)
        .map { row ->
            def meta       = SampleUtils.makeMeta(row)
            def genome_col = row.GENOME?.trim()
            def gz = genome_col
                ? (genome_col.startsWith('/') ? file(genome_col) : file("${launchDir}/${genome_col}"))
                : file("${params.source}/${meta.asmid}/${meta.asmid}_genomic.fna.gz")
            tuple(meta, gz)
        }
        .filter { meta, gz -> meta.id && meta.asmid }
        .take((params.n_test as int) > 0 ? params.n_test as int : -1)
        .filter { meta, gz ->
            if (suppressSet.contains(meta.asmid)) {
                log.info "Suppressing ${meta.id} (asmid=${meta.asmid})"
                return false
            }
            return true
        }

    // ── genomes: additionally require the genome file to exist ────────────────
    def ch_genomes = ch_raw
        .filter { meta, gz ->
            if (!gz.exists()) {
                log.warn "Missing genome for ${meta.id} (asmid=${meta.asmid}): ${gz}"
                return false
            }
            if (params.debug.toBoolean()) {
                log.info "Queuing ${meta.id}: genome=${gz} (${gz.size()} bytes)"
            }
            return true
        }

    emit:
    samples = ch_raw.map { meta, _gz -> meta }
    genomes = ch_genomes
}
