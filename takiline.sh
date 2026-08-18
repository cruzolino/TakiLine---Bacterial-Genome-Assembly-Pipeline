#!/usr/bin/env bash
# ==============================================================================
# TakiLine — Bacterial Genome Assembly Pipeline
# Modes: Illumina (PE) | Nanopore | Hybrid (Illumina + Nanopore) | PacBio HiFi
# ==============================================================================
set -Eeuo pipefail
IFS=$'\n\t'

# ── Defaults ───────────────────────────────────────────────────────────────────
THREADS=$(nproc --all 2>/dev/null || echo 8)
MEMORY="32G"
OUTDIR="takiline"
ASSEMBLER="spades"
READ_TYPE="illumina"
MIN_CONTIG_LENGTH=500
GENOME_SIZE="5m"
SAMPLE="isolate"
PLASMID_DETECTION=false   # opt-in via --plasmid; see run_quality_assessment()
PLASMID_DB="enterobacteriaceae"
PLASMID_DB_PATH=""   # empty = auto-detect at runtime; see run_quality_assessment()
MEDAKA_MODEL=""       # empty = auto-select via --bacteria; see run_polishing()
BUSCO_LINEAGE=""       # empty = auto-select via --auto-lineage-prok; see run_quality_assessment()
AMR_DETECTION=false        # opt-in via --amr (AMRFinderPlus)
ABRICATE_DETECTION=false   # opt-in via --abricate (secondary AMR screen)
ABRICATE_DB="card"         # not "ncbi" — that mirrors AMRFinderPlus's own db; see README
MLST_DETECTION=false       # opt-in via --mlst
SPECIES=""             # -S "Genus species"; empty = auto-detect via Kraken2
KRAKEN2_DB=""          # -K DIR, falls back to $KRAKEN2_DB_PATH
QC_ONLY=false
HIFI_MODE=false
RESUME=false
CLEANUP=true

# ── Colors ─────────────────────────────────────────────────────────────────────
readonly RED='\033[0;31m' GREEN='\033[0;32m' YELLOW='\033[1;33m' \
         BLUE='\033[0;34m' CYAN='\033[0;36m' NC='\033[0m'

# ── Logging ────────────────────────────────────────────────────────────────────
log()  { echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $*"; }
ok()   { echo -e "${GREEN}[✓ $(date +%H:%M:%S)]${NC} $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
err()  { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

# ── Derived-path helpers ───────────────────────────────────────────────────────
# Centralized so a resumed run still has these set even when run_qc is skipped.
_set_derived_paths() {
    TRIMMED_R1="${OUTDIR}/qc/${SAMPLE}_R1.trimmed.fastq.gz"
    TRIMMED_R2="${OUTDIR}/qc/${SAMPLE}_R2.trimmed.fastq.gz"
    FILTERED_LONG="${OUTDIR}/qc/${SAMPLE}_filtered.fastq.gz"
}

# Input filename, minus .fastq/.fq(.gz) — used to name the final assembly. Prefers
# LONG (the primary input in long-read/hybrid mode) and falls back to R1.
_input_basename() {
    local b; b=$(basename "${LONG:-${R1}}")
    b="${b%.gz}"; b="${b%.fastq}"; b="${b%.fq}"
    echo "${b}"
}

# ── Sentinel helpers (resume logic) ───────────────────────────────────────────
sentinel_path() { echo "${OUTDIR}/logs/.done_${1}"; }

stage_done() {
    local stage="$1"
    touch "$(sentinel_path "${stage}")"
}

skip_if_done() {
    local stage="$1"
    if ${RESUME} && [[ -f "$(sentinel_path "${stage}")" ]]; then
        ok "Skipping ${stage} (already completed; use without -r to re-run)"
        return 0
    fi
    return 1
}

# ── Resume-consistency fingerprint ────────────────────────────────────────────
# Detects -r resuming with different -s/-a/-g/-1/-2/-l/-S than the run that
# produced the sentinels in OUTDIR. Read back positionally, not sourced/eval'd.
resume_fingerprint_path() { echo "${OUTDIR}/logs/.run_fingerprint"; }

write_resume_fingerprint() {
    printf '%s\n' "${SAMPLE}" "${ASSEMBLER}" "${READ_TYPE}" "${GENOME_SIZE}" "${HIFI_MODE}" \
        "$(basename "${R1:-}")" "$(basename "${R2:-}")" "$(basename "${LONG:-}")" "${SPECIES}" \
        > "$(resume_fingerprint_path)"
}

check_resume_fingerprint() {
    local fp; fp="$(resume_fingerprint_path)"
    [[ -f "${fp}" ]] || return 0

    local prev_sample prev_assembler prev_read_type prev_genome_size prev_hifi_mode prev_r1 prev_r2 prev_long prev_species
    { read -r prev_sample; read -r prev_assembler; read -r prev_read_type; \
      read -r prev_genome_size; read -r prev_hifi_mode; \
      read -r prev_r1; read -r prev_r2; read -r prev_long; read -r prev_species; \
    } < "${fp}"

    local cur_r1 cur_r2 cur_long
    cur_r1="$(basename "${R1:-}")"; cur_r2="$(basename "${R2:-}")"; cur_long="$(basename "${LONG:-}")"

    if [[ "${prev_sample}" != "${SAMPLE}" || "${prev_assembler}" != "${ASSEMBLER}" || \
          "${prev_read_type}" != "${READ_TYPE}" || "${prev_genome_size}" != "${GENOME_SIZE}" || \
          "${prev_hifi_mode}" != "${HIFI_MODE}" || "${prev_species}" != "${SPECIES}" || \
          "${prev_r1}" != "${cur_r1}" || "${prev_r2}" != "${cur_r2}" || "${prev_long}" != "${cur_long}" ]]; then
        err "Resume (-r) flags don't match the run recorded in ${OUTDIR}: sample='${prev_sample}'→'${SAMPLE}', assembler='${prev_assembler}'→'${ASSEMBLER}', read-type='${prev_read_type}'→'${READ_TYPE}', genome-size='${prev_genome_size}'→'${GENOME_SIZE}', hifi='${prev_hifi_mode}'→'${HIFI_MODE}', species='${prev_species}'→'${SPECIES}'. Use a different -o for a new run, or match the original flags to resume this one."
    fi
}

# ── Parallel gzip ──────────────────────────────────────────────────────────────
# pigz when available, else gzip. Accepts the same flags as gzip.
_gzip() {
    if command -v pigz &>/dev/null; then
        pigz -p "${THREADS}" "$@"
    else
        gzip "$@"
    fi
}

# ── Genome-size conversion ──────────────────────────────────────────────────────
# "5m"/"0.5g" -> raw base count. Returns 1 on bad input; caller must check status.
genome_size_to_bases() {
    local gs="${1,,}"
    if [[ "${gs}" =~ ^([0-9]+(\.[0-9]+)?)m$ ]]; then
        awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]} * 1000000}"
    elif [[ "${gs}" =~ ^([0-9]+(\.[0-9]+)?)g$ ]]; then
        awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]} * 1000000000}"
    else
        return 1
    fi
}

# ── Memory-limit conversion ─────────────────────────────────────────────────────
# "-m" shorthand (e.g. "32G") -> raw integer GB. Same contract as genome_size_to_bases().
memory_to_gb_int() {
    local mem="${1,,}"
    if [[ "${mem}" =~ ^([0-9]+(\.[0-9]+)?)g$ ]]; then
        awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]}}"
    else
        return 1
    fi
}

# ── Tool requirement check ─────────────────────────────────────────────────────
check_tools() {
    local missing=()
    for t in "$@"; do
        command -v "${t}" &>/dev/null || missing+=("${t}")
    done
    if (( ${#missing[@]} > 0 )); then
        err "Missing required tools: ${missing[*]}. Install them and retry."
    fi
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
${GREEN}TakiLine${NC}

${YELLOW}Usage:${NC}  $0 [OPTIONS]

${YELLOW}Input (at least one required):${NC}
  -1 FILE   Illumina forward reads (R1.fastq.gz)
  -2 FILE   Illumina reverse reads (R2.fastq.gz)
  -l FILE   Long reads – Nanopore or PacBio HiFi (.fastq[.gz])
  --hifi    Treat -l reads as PacBio HiFi (activates flye --pacbio-hifi)

${YELLOW}Assembly:${NC}
  -a STR    Assembler: spades|skesa|unicycler (Illumina/hybrid)
                       flye|canu|raven|trycycler (long-read) [default: spades]
            trycycler: 12-subsample consensus; auto-drops clusters that fail
            reconcile (review ${OUTDIR}/02_assembly/trycycler/cluster/).

${YELLOW}General:${NC}
  -o DIR    Output directory          [default: ${OUTDIR}]
  -t INT    Threads                   [default: auto = ${THREADS}]
  -m STR    Memory limit              [default: ${MEMORY}]
  -s STR    Sample name               [default: ${SAMPLE}]
  -g STR    Expected genome size      [default: ${GENOME_SIZE}]
  -c INT    Min contig length (bp)    [default: ${MIN_CONTIG_LENGTH}]
  -S STR    Expected species, e.g. "Escherichia coli" [default: auto-detect via Kraken2]
            Recommended — skips Kraken2. If omitted, Kraken2 is REQUIRED.
  -K DIR    Kraken2 database path     [default: \$KRAKEN2_DB_PATH env var]
  -P STR    PlasmidFinder DB          [default: ${PLASMID_DB}]
  -D DIR    PlasmidFinder DB path     [default: auto-detect]
  -M STR    Medaka model (Nanopore polishing) [default: auto-select via --bacteria]
            e.g. -M r941_min_sup_g507 for legacy/R9.4.1 reads.
  -L STR    BUSCO lineage odb10 [default: auto-select via --auto-lineage-prok]
            e.g. -L enterobacterales_odb10 to force a specific lineage.
  -q        QC only (skip assembly)
  -r        Resume from last checkpoint
  -k        Keep intermediates (raw-data symlinks, trimmed/filtered FASTQs,
            assembler working dirs) instead of cleaning up after a full run
  -h|--help Show this help

${YELLOW}Optional typing/resistance analyses (off by default, hard requirement once passed):${NC}
  --plasmid   PlasmidFinder replicon typing (uses -P/-D above)
  --amr       AMRFinderPlus (AMR genes + point mutations, --plus panel)
  --abricate  abricate secondary AMR screen (db: ${ABRICATE_DB})
  --mlst      mlst sequence typing (PubMLST schemes)

${YELLOW}Examples:${NC}
  # Illumina PE, species informed (skips Kraken2)
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -S "Escherichia coli" -t 16

  # Illumina PE, species unknown — Kraken2 required (-K or \$KRAKEN2_DB_PATH)
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -K /path/to/kraken2_db

  # Nanopore
  $0 -l ont.fq.gz -s Salmonella -a flye -g 4.8m -t 16

  # Hybrid
  $0 -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz -s Klebsiella -a unicycler

  # PacBio HiFi
  $0 -l hifi.fastq.gz -s Pseudomonas -a flye --hifi -g 6.5m

  # Resume interrupted run
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r

  # Full typing workup: assembly + plasmids + AMR (both tools) + MLST
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 --plasmid --amr --abricate --mlst
EOF
    exit 0
}

# ── Argument parsing ───────────────────────────────────────────────────────────
ARGS=()
for arg in "$@"; do
    case "${arg}" in
        --hifi)     HIFI_MODE=true ;;
        --plasmid)  PLASMID_DETECTION=true ;;
        --amr)      AMR_DETECTION=true ;;
        --abricate) ABRICATE_DETECTION=true ;;
        --mlst)     MLST_DETECTION=true ;;
        --help)     usage ;;
        *)          ARGS+=("${arg}") ;;
    esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

while getopts "1:2:l:o:t:s:g:a:m:c:S:K:P:D:M:L:qrkh" opt; do
    case "${opt}" in
        1) R1="${OPTARG}" ;;
        2) R2="${OPTARG}" ;;
        l) LONG="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        s) SAMPLE="${OPTARG}" ;;
        g) GENOME_SIZE="${OPTARG}" ;;
        a) ASSEMBLER="${OPTARG}" ;;
        m) MEMORY="${OPTARG}" ;;
        c) MIN_CONTIG_LENGTH="${OPTARG}" ;;
        S) SPECIES="${OPTARG}" ;;
        K) KRAKEN2_DB="${OPTARG}" ;;
        P) PLASMID_DB="${OPTARG}" ;;
        D) PLASMID_DB_PATH="${OPTARG}" ;;
        M) MEDAKA_MODEL="${OPTARG}" ;;
        L) BUSCO_LINEAGE="${OPTARG}" ;;
        q) QC_ONLY=true ;;
        r) RESUME=true ;;
        k) CLEANUP=false ;;
        h) usage ;;
        *) err "Invalid option: -${OPTARG-}" ;;
    esac
done

# ── Validate inputs ────────────────────────────────────────────────────────────
validate_inputs() {
    log "Validating inputs and environment..."

    [[ -z "${R1:-}" && -z "${LONG:-}" ]] && err "No input reads provided. Use -h for help."

    if { [[ -n "${R1:-}" ]] && [[ -z "${R2:-}" ]]; } ||
       { [[ -z "${R1:-}" ]] && [[ -n "${R2:-}" ]]; }; then
        err "Both -1 and -2 are required for paired-end Illumina reads."
    fi

    if [[ -n "${LONG:-}" && -n "${R1:-}" ]]; then
        READ_TYPE="hybrid"
        log "Mode: ${CYAN}Hybrid${NC} (Illumina + Long-read)"
    elif [[ -n "${LONG:-}" ]]; then
        READ_TYPE="long"
        log "Mode: ${CYAN}Long-read${NC}$(${HIFI_MODE} && echo ' [PacBio HiFi]' || echo ' [Nanopore]')"
    else
        READ_TYPE="illumina"
        log "Mode: ${CYAN}Illumina paired-end${NC}"
    fi

    [[ -n "${R1:-}"   && ! -f "${R1}"   ]] && err "R1 file not found: ${R1}"
    [[ -n "${R2:-}"   && ! -f "${R2}"   ]] && err "R2 file not found: ${R2}"
    [[ -n "${LONG:-}" && ! -f "${LONG}" ]] && err "Long-read file not found: ${LONG}"

    case "${ASSEMBLER}" in
        spades|skesa|unicycler)
            [[ "${READ_TYPE}" == "long" ]] && \
                err "${ASSEMBLER} cannot assemble long-reads only. Use flye, canu, or raven."
            if [[ "${READ_TYPE}" == "hybrid" && "${ASSEMBLER}" != "unicycler" ]]; then
                err "${ASSEMBLER} does not use long reads — the -l input would be silently discarded in hybrid mode. Use -a unicycler for hybrid assembly, or drop -l to run ${ASSEMBLER} on Illumina reads only."
            fi
            ;;
        flye|canu|raven|trycycler)
            [[ "${READ_TYPE}" == "illumina" ]] && \
                err "${ASSEMBLER} is a long-read assembler. For Illumina use spades, skesa, or unicycler."
            [[ "${READ_TYPE}" == "hybrid" ]] && \
                warn "${ASSEMBLER} assembles only the long reads; Illumina reads will be used for post-assembly Pilon polishing, not integrated hybrid assembly. Use -a unicycler for true hybrid assembly."
            ;;
        *) err "Unknown assembler: ${ASSEMBLER}" ;;
    esac

    [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || err "Threads must be a positive integer."
    [[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]] || err "Min contig length must be a positive integer."
    [[ "${GENOME_SIZE,,}" =~ ^[0-9]+(\.[0-9]+)?[mg]$ ]] || \
        err "Genome size must be a number followed by m or g (e.g. 5m, 4.8m, 0.5g)."
    memory_to_gb_int "${MEMORY}" &>/dev/null || \
        err "Memory limit must be a number in gigabytes followed by G (e.g. 32G, 1.5G)."

    # ── Preflight tool checks (all at once) ───────────────────────────────────
    local required_tools=(fastqc fastp quast realpath)
    [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]] && \
        required_tools+=(filtlong)
    case "${ASSEMBLER}" in
        spades)    required_tools+=(spades.py) ;;
        skesa)     required_tools+=(skesa) ;;
        unicycler) required_tools+=(unicycler) ;;
        flye)      required_tools+=(flye) ;;
        canu)      required_tools+=(canu) ;;
        raven)     required_tools+=(raven) ;;
        trycycler) required_tools+=(trycycler flye miniasm minipolish raven minimap2) ;;
    esac
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then
        required_tools+=(bowtie2 samtools pilon)
    fi
    if [[ "${READ_TYPE}" == "long" && "${HIFI_MODE}" == false ]]; then
        required_tools+=(medaka_consensus)
    fi
    # Species: informing it (-S) skips detection entirely; otherwise Kraken2 is a
    # hard requirement (not best-effort) — an unidentified/uninformed sample isn't
    # a case the pipeline should silently proceed past.
    [[ -z "${SPECIES}" ]] && required_tools+=(kraken2)
    # Opt-in typing/resistance tools: requested via flag => hard requirement,
    # not a best-effort skip (unlike CheckM2/BUSCO, which stay opportunistic).
    ${PLASMID_DETECTION}  && required_tools+=(plasmidfinder.py)
    ${AMR_DETECTION}      && required_tools+=(amrfinder)
    ${ABRICATE_DETECTION} && required_tools+=(abricate)
    ${MLST_DETECTION}     && required_tools+=(mlst)
    check_tools "${required_tools[@]}"

    # Kraken2 needs a database, not just the binary — only checked when it'll actually run.
    if [[ -z "${SPECIES}" ]]; then
        [[ -z "${KRAKEN2_DB}" ]] && KRAKEN2_DB="${KRAKEN2_DB_PATH:-}"
        [[ -n "${KRAKEN2_DB}" && -d "${KRAKEN2_DB}" ]] || \
            err "No species given (-S) and no Kraken2 database found. Pass -K <db_dir> (or set \$KRAKEN2_DB_PATH), or supply -S \"Genus species\" to skip auto-detection."
    fi

    # ── Directory scaffold ────────────────────────────────────────────────────
    mkdir -p "${OUTDIR}"/{00_raw_data,qc/post_trim,02_assembly,\
03_polishing,reports/{quast,checkm,busco,plasmidfinder,amrfinderplus,abricate,mlst,kraken2},logs}

    [[ -n "${R1:-}"   ]] && ln -sf "$(realpath "${R1}")"   "${OUTDIR}/00_raw_data/$(basename "${R1}")"
    [[ -n "${R2:-}"   ]] && ln -sf "$(realpath "${R2}")"   "${OUTDIR}/00_raw_data/$(basename "${R2}")"
    [[ -n "${LONG:-}" ]] && ln -sf "$(realpath "${LONG}")" "${OUTDIR}/00_raw_data/$(basename "${LONG}")"

    # ── Resume-consistency check ──────────────────────────────────────────────
    ${RESUME} && check_resume_fingerprint
    write_resume_fingerprint

    _set_derived_paths

    ok "Input validation passed"
}

# ── QC ─────────────────────────────────────────────────────────────────────────
run_qc() {
    skip_if_done "qc" && return

    log "[1/5] Quality control..."

    # ── Illumina QC ──────────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then

        # Start FastQC on raw reads in background while fastp runs in parallel.
        log "FastQC on raw reads (background) + fastp trimming (foreground)..."
        fastqc "${R1}" "${R2}" \
            --outdir "${OUTDIR}/qc/" \
            --threads "${THREADS}" \
            --quiet \
            > "${OUTDIR}/logs/fastqc_raw.log" 2>&1 &
        local fastqc_raw_pid=$!

        fastp \
            -i "${R1}" -I "${R2}" \
            -o "${TRIMMED_R1}" -O "${TRIMMED_R2}" \
            --thread "${THREADS}" \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --length_required 50 \
            --detect_adapter_for_pe \
            --correction \
            --overrepresentation_analysis \
            --html "${OUTDIR}/qc/${SAMPLE}_fastp.html" \
            --json "${OUTDIR}/qc/${SAMPLE}_fastp.json" \
            --report_title "${SAMPLE}" \
            > "${OUTDIR}/logs/fastp.log" 2>&1

        # Wait for raw FastQC before running FastQC on trimmed reads
        wait "${fastqc_raw_pid}" || warn "FastQC (raw) finished with warnings; check fastqc_raw.log."

        log "FastQC on trimmed reads..."
        fastqc "${TRIMMED_R1}" "${TRIMMED_R2}" \
            --outdir "${OUTDIR}/qc/post_trim/" \
            --threads "${THREADS}" \
            --quiet \
            > "${OUTDIR}/logs/fastqc_post.log" 2>&1

        if command -v multiqc &>/dev/null; then
            multiqc "${OUTDIR}/qc/" \
                --outdir "${OUTDIR}/reports/" \
                --filename "${SAMPLE}_multiqc" \
                --no-megaqc-upload \
                --quiet \
                > "${OUTDIR}/logs/multiqc.log" 2>&1
        else
            warn "multiqc not found; skipping MultiQC report."
        fi
    fi

    # ── Long-read QC ──────────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]]; then

        # Trycycler needs a deeper pool (~100x) than a single assembler (~20x) for independent subsamples.
        local _gs_bases
        _gs_bases=$(genome_size_to_bases "${GENOME_SIZE}") || \
            err "Cannot parse genome size '${GENOME_SIZE}' for TARGET_BASES calculation."
        local _target_cov=20
        [[ "${ASSEMBLER}" == "trycycler" ]] && _target_cov=100
        local TARGET_BASES=$(( _gs_bases * _target_cov ))

        # Start NanoPlot in background while filtlong runs in foreground.
        local nanoplot_pid=""
        if command -v NanoPlot &>/dev/null; then
            log "NanoPlot QC (background) + Filtlong filtering (foreground)..."
            NanoPlot \
                --fastq "${LONG}" \
                --outdir "${OUTDIR}/qc/nanoplot/" \
                --threads "${THREADS}" \
                --plots dot kde \
                > "${OUTDIR}/logs/nanoplot.log" 2>&1 &
            nanoplot_pid=$!
        else
            warn "NanoPlot not found; skipping long-read QC plots."
            log "Filtering long reads with Filtlong (target: ${TARGET_BASES} bases)..."
        fi

        # Quality threshold differs by platform: Nanopore Q7, HiFi Q20.
        local _filtlong_q=7
        ${HIFI_MODE} && _filtlong_q=20

        filtlong \
            --min_length 1000 \
            --min_mean_q "${_filtlong_q}" \
            --target_bases "${TARGET_BASES}" \
            "${LONG}" \
            2> "${OUTDIR}/logs/filtlong.log" \
        | _gzip -1 -c > "${FILTERED_LONG}"

        if [[ -n "${nanoplot_pid}" ]]; then
            wait "${nanoplot_pid}" || warn "NanoPlot finished with warnings; check nanoplot.log."
        fi
    fi

    ok "QC complete"
    stage_done "qc"
}

# ── Species identification ────────────────────────────────────────────────────
# -S skips this. Otherwise Kraken2 is a hard requirement (checked in validate_inputs()).
run_species_id() {
    skip_if_done "species_id" && {
        _restore_species
        return
    }

    log "[2/5] Species identification..."

    if [[ -n "${SPECIES}" ]]; then
        ok "Species informed: ${SPECIES} (Kraken2 skipped)"
    else
        log "No species informed (-S); running Kraken2 to auto-detect..."
        local kraken_report="${OUTDIR}/reports/kraken2/${SAMPLE}_kraken2_report.tsv"
        local kraken_reads_args=()
        if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then
            kraken_reads_args=(--paired "${TRIMMED_R1}" "${TRIMMED_R2}")
        else
            kraken_reads_args=("${FILTERED_LONG}")
        fi

        kraken2 \
            --db "${KRAKEN2_DB}" \
            --threads "${THREADS}" \
            --report "${kraken_report}" \
            --output "${OUTDIR}/reports/kraken2/${SAMPLE}_kraken2_output.tsv" \
            "${kraken_reads_args[@]}" \
            > "${OUTDIR}/logs/kraken2.log" 2>&1

        # Report columns: pct, clade_reads, direct_reads, rank_code, taxid, name
        # (name is indented by tree depth; rank_code "S" = species level).
        local top_line
        top_line=$(awk -F'\t' '$4=="S"' "${kraken_report}" | sort -t $'\t' -k1,1 -rn | head -1)
        [[ -n "${top_line}" ]] || \
            err "Kraken2 produced no species-level classification for ${SAMPLE}. Check ${OUTDIR}/logs/kraken2.log and ${kraken_report} — or pass -S \"Genus species\" to skip auto-detection."

        local top_pct top_name
        top_pct=$(awk -F'\t' '{gsub(/^ +| +$/,"",$1); print $1}' <<< "${top_line}")
        top_name=$(awk -F'\t' '{gsub(/^ +| +$/,"",$6); print $6}' <<< "${top_line}")
        SPECIES="${top_name}"

        ok "Kraken2 species call: ${SPECIES} (${top_pct}% of classified reads)"
        if awk "BEGIN{exit !(${top_pct} < 70)}"; then
            warn "Top species (${SPECIES}) accounts for only ${top_pct}% of classified reads — possible contamination or a mixed sample. Review ${kraken_report}."
        fi
    fi

    echo "${SPECIES}" > "${OUTDIR}/logs/.species"
    stage_done "species_id"
}

# Restore SPECIES when resuming past the species-id stage.
_restore_species() {
    local species_file="${OUTDIR}/logs/.species"
    [[ -f "${species_file}" ]] && SPECIES=$(cat "${species_file}")
    log "Resumed with species: ${SPECIES:-not determined}"
}

# Drops contigs >10% off the cluster's median length (mirrors reconcile's own check).
_trycycler_prune_length_outliers() {
    local contigs_dir=$1
    local f len lengths=() sorted median ratio
    for f in "${contigs_dir}"/*.fasta; do
        len=$(seqkit stats -T "${f}" | awk 'NR==2{print $5}')
        lengths+=("${len}")
    done
    sorted=($(printf '%s\n' "${lengths[@]}" | sort -n))
    median=${sorted[$(( ${#sorted[@]} / 2 ))]}
    for f in "${contigs_dir}"/*.fasta; do
        len=$(seqkit stats -T "${f}" | awk 'NR==2{print $5}')
        ratio=$(LC_NUMERIC=C awk -v l="${len}" -v m="${median}" 'BEGIN{d=l/m; if(d<1)d=1/d; printf "%.3f", d}')
        if awk "BEGIN{exit !(${ratio} > 1.1)}"; then
            warn "  dropping outlier contig $(basename "${f}") (length ${len} vs cluster median ${median}, ratio ${ratio})"
            rm -f "${f}"
        fi
    done
}

# Drops contigs reconcile names as failing (e.g. circularisation), then retries.
_trycycler_drop_named_contigs() {
    local contigs_dir=$1; shift
    local name f
    for name in "$@"; do
        f="${contigs_dir}/${name}.fasta"
        [[ -f "${f}" ]] || continue
        warn "  dropping contig ${name} (named by Trycycler reconcile as failing)"
        rm -f "${f}"
    done
}

# ── Assembly ───────────────────────────────────────────────────────────────────
run_assembly() {
    skip_if_done "assembly" && {
        _restore_assembly_path
        return
    }

    log "[3/5] Assembly with ${ASSEMBLER}..."

    local MEM_INT
    MEM_INT=$(memory_to_gb_int "${MEMORY}") || err "Cannot parse memory limit '${MEMORY}'."

    case "${ASSEMBLER}" in

        spades)
            spades.py \
                -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}" \
                -o "${OUTDIR}/02_assembly/spades/" \
                -t "${THREADS}" -m "${MEM_INT}" \
                --isolate \
                --cov-cutoff auto \
                > "${OUTDIR}/logs/spades.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/spades/contigs.fasta"
            ;;

        skesa)
            mkdir -p "${OUTDIR}/02_assembly/skesa"
            skesa \
                --fastq "${TRIMMED_R1},${TRIMMED_R2}" \
                --cores "${THREADS}" \
                --memory "${MEM_INT}" \
                --contigs_out "${OUTDIR}/02_assembly/skesa/contigs.fasta" \
                > "${OUTDIR}/logs/skesa.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/skesa/contigs.fasta"
            ;;

        unicycler)
            local uni_args=(-1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}"
                            -o "${OUTDIR}/02_assembly/unicycler/"
                            -t "${THREADS}")
            [[ "${READ_TYPE}" == "hybrid" ]] && \
                uni_args+=(-l "${FILTERED_LONG}" --mode normal) || \
                uni_args+=(--mode conservative)
            unicycler "${uni_args[@]}" > "${OUTDIR}/logs/unicycler.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/unicycler/assembly.fasta"
            ;;

        flye)
            # --nano-hq requires Flye >= 2.9; auto-detect and fall back to --nano-raw.
            local flye_read_flag
            if ${HIFI_MODE}; then
                flye_read_flag="--pacbio-hifi"
            else
                local flye_version=""
                flye_version=$(LC_ALL=C flye --version 2>&1 | grep -Eo '[0-9]+\.[0-9]+' | head -1) || true
                if awk "BEGIN{exit !(${flye_version:-0} >= 2.9)}"; then
                    flye_read_flag="--nano-hq"
                else
                    warn "Flye ${flye_version} detected (< 2.9); using --nano-raw. Upgrade to Flye >= 2.9 for --nano-hq."
                    flye_read_flag="--nano-raw"
                fi
            fi
            flye \
                ${flye_read_flag} "${FILTERED_LONG}" \
                --out-dir "${OUTDIR}/02_assembly/flye/" \
                --threads "${THREADS}" \
                --genome-size "${GENOME_SIZE}" \
                --iterations 3 \
                > "${OUTDIR}/logs/flye.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/flye/assembly.fasta"
            ;;

        canu)
            local canu_read_flag="-nanopore"
            ${HIFI_MODE} && canu_read_flag="-pacbio-hifi"
            canu \
                -p "${SAMPLE}" \
                -d "${OUTDIR}/02_assembly/canu/" \
                genomeSize="${GENOME_SIZE}" \
                useGrid=false \
                maxThreads="${THREADS}" \
                maxMemory="${MEMORY}" \
                ${canu_read_flag} "${FILTERED_LONG}" \
                > "${OUTDIR}/logs/canu.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/canu/${SAMPLE}.contigs.fasta"
            ;;

        raven)
            mkdir -p "${OUTDIR}/02_assembly/raven"
            raven \
                --threads "${THREADS}" \
                --graphical-fragment-assembly "${OUTDIR}/02_assembly/raven/assembly.gfa" \
                "${FILTERED_LONG}" \
                > "${OUTDIR}/02_assembly/raven/assembly.fasta" \
                2> "${OUTDIR}/logs/raven.log"
            ASSEMBLY="${OUTDIR}/02_assembly/raven/assembly.fasta"
            ;;

        trycycler)
            local tc_dir="${OUTDIR}/02_assembly/trycycler"
            mkdir -p "${tc_dir}"

            log "Trycycler [1/6]: subsampling reads into independent read sets..."
            local tc_genome_bp
            tc_genome_bp=$(genome_size_to_bases "${GENOME_SIZE}") || \
                err "Cannot parse genome size '${GENOME_SIZE}' for Trycycler."
            trycycler subsample \
                -r "${FILTERED_LONG}" \
                -o "${tc_dir}/subsampled" \
                --genome_size "${tc_genome_bp}" \
                -t "${THREADS}" \
                > "${OUTDIR}/logs/trycycler_subsample.log" 2>&1
            local sub_fastqs=("${tc_dir}/subsampled"/sample_*.fastq)
            [[ -e "${sub_fastqs[0]}" ]] || \
                err "Trycycler subsample produced no read subsets. Check ${OUTDIR}/logs/trycycler_subsample.log"

            log "Trycycler [2/6]: assembling ${#sub_fastqs[@]} subsamples (Flye | Miniasm+Minipolish | Raven, rotated)..."
            mkdir -p "${tc_dir}/assemblies"
            local i sub_assembler sub_reads sub_id flye_flag
            for i in "${!sub_fastqs[@]}"; do
                sub_reads="${sub_fastqs[$i]}"
                sub_id=$(printf "%02d" $((i + 1)))
                case $(( i % 3 )) in
                    0) sub_assembler="flye" ;;
                    1) sub_assembler="miniasm" ;;
                    2) sub_assembler="raven" ;;
                esac
                if [[ -s "${tc_dir}/assemblies/assembly_${sub_id}.fasta" ]]; then
                    log "  subsample ${sub_id}/${#sub_fastqs[@]}: ${sub_assembler} — already assembled, skipping (resume)"
                    continue
                fi
                log "  subsample ${sub_id}/${#sub_fastqs[@]}: ${sub_assembler}"
                case "${sub_assembler}" in
                    flye)
                        flye_flag="--nano-hq"
                        ${HIFI_MODE} && flye_flag="--pacbio-hifi"
                        flye ${flye_flag} "${sub_reads}" \
                            --out-dir "${tc_dir}/assemblies/tmp_${sub_id}" \
                            --threads "${THREADS}" --genome-size "${GENOME_SIZE}" \
                            > "${OUTDIR}/logs/trycycler_assemble_${sub_id}.log" 2>&1 \
                        && cp "${tc_dir}/assemblies/tmp_${sub_id}/assembly.fasta" \
                              "${tc_dir}/assemblies/assembly_${sub_id}.fasta"
                        ;;
                    miniasm)
                        local miniasm_ava_flag="ava-ont"
                        ${HIFI_MODE} && miniasm_ava_flag="ava-pb"
                        minimap2 -x "${miniasm_ava_flag}" -t "${THREADS}" "${sub_reads}" "${sub_reads}" \
                            2> "${OUTDIR}/logs/trycycler_assemble_${sub_id}.log" \
                            > "${tc_dir}/assemblies/tmp_${sub_id}.paf" \
                        && miniasm -f "${sub_reads}" "${tc_dir}/assemblies/tmp_${sub_id}.paf" \
                            > "${tc_dir}/assemblies/tmp_${sub_id}.gfa" \
                            2>> "${OUTDIR}/logs/trycycler_assemble_${sub_id}.log" \
                        && minipolish --threads "${THREADS}" "${sub_reads}" "${tc_dir}/assemblies/tmp_${sub_id}.gfa" \
                            > "${tc_dir}/assemblies/tmp_${sub_id}_polished.gfa" \
                            2>> "${OUTDIR}/logs/trycycler_assemble_${sub_id}.log" \
                        && awk '/^S/{print ">"$2"\n"$3}' "${tc_dir}/assemblies/tmp_${sub_id}_polished.gfa" \
                            > "${tc_dir}/assemblies/assembly_${sub_id}.fasta"
                        ;;
                    raven)
                        raven --threads "${THREADS}" "${sub_reads}" \
                            > "${tc_dir}/assemblies/assembly_${sub_id}.fasta" \
                            2> "${OUTDIR}/logs/trycycler_assemble_${sub_id}.log"
                        ;;
                esac
                [[ -s "${tc_dir}/assemblies/assembly_${sub_id}.fasta" ]] || {
                    warn "Subsample ${sub_id} (${sub_assembler}) produced no assembly; skipping it. Check ${OUTDIR}/logs/trycycler_assemble_${sub_id}.log"
                    rm -f "${tc_dir}/assemblies/assembly_${sub_id}.fasta"
                }
            done

            local n_sub_ok
            n_sub_ok=$(find "${tc_dir}/assemblies" -maxdepth 1 -name 'assembly_*.fasta' | wc -l)
            (( n_sub_ok >= 3 )) || \
                err "Only ${n_sub_ok} of ${#sub_fastqs[@]} subsample assemblies succeeded; Trycycler needs multiple independent assemblies to reconcile. Check ${OUTDIR}/logs/trycycler_assemble_*.log"

            log "Trycycler [3/6]: clustering contigs across ${n_sub_ok} assemblies..."
            trycycler cluster \
                -a "${tc_dir}/assemblies"/assembly_*.fasta \
                -r "${FILTERED_LONG}" \
                -o "${tc_dir}/cluster" \
                -t "${THREADS}" \
                > "${OUTDIR}/logs/trycycler_cluster.log" 2>&1
            [[ -d "${tc_dir}/cluster" ]] || \
                err "Trycycler cluster failed. Check ${OUTDIR}/logs/trycycler_cluster.log"

            log "Trycycler [4/6]: reconciling clusters (failing ones auto-dropped; review ${tc_dir}/cluster/ afterward)..."
            local cluster_dirs=() cdir cname
            for cdir in "${tc_dir}/cluster"/cluster_*; do
                [[ -d "${cdir}" ]] || continue
                cname=$(basename "${cdir}")
                # reconcile's exit code is unreliable (0 on some failures, 1 on others); check its output file instead.
                local rc_log="${OUTDIR}/logs/trycycler_reconcile_${cname}.log"
                trycycler reconcile -c "${cdir}" -r "${FILTERED_LONG}" -t "${THREADS}" > "${rc_log}" 2>&1 || true
                local reconciled=false attempt=1 pruned n_remaining
                [[ -s "${cdir}/2_all_seqs.fasta" ]] && reconciled=true
                while ! ${reconciled} && (( attempt <= 3 )); do
                    pruned=false
                    if grep -q "too much length difference" "${rc_log}"; then
                        _trycycler_prune_length_outliers "${cdir}/1_contigs"
                        pruned=true
                    fi
                    if grep -q "failed to circularise" "${rc_log}"; then
                        local bad_names_arr=()
                        mapfile -t bad_names_arr < <(grep -oP 'failed to circularise: \K[^.]*' "${rc_log}" | tr ',' '\n' | sed 's/^ *//;s/ *$//')
                        _trycycler_drop_named_contigs "${cdir}/1_contigs" "${bad_names_arr[@]}"
                        pruned=true
                    fi
                    ${pruned} || break
                    n_remaining=$(find "${cdir}/1_contigs" -maxdepth 1 -name '*.fasta' | wc -l)
                    (( n_remaining >= 2 )) || break
                    attempt=$(( attempt + 1 ))
                    rc_log="${OUTDIR}/logs/trycycler_reconcile_${cname}_attempt${attempt}.log"
                    warn "Retrying Trycycler reconcile for ${cname} (attempt ${attempt}) after automatic contig pruning..."
                    trycycler reconcile -c "${cdir}" -r "${FILTERED_LONG}" -t "${THREADS}" > "${rc_log}" 2>&1 || true
                    [[ -s "${cdir}/2_all_seqs.fasta" ]] && reconciled=true
                done
                if ${reconciled}; then
                    (( attempt > 1 )) && warn "Trycycler reconcile succeeded for ${cname} after ${attempt} attempt(s) with automatic contig pruning."
                    cluster_dirs+=("${cdir}")
                else
                    warn "Trycycler reconcile failed for ${cname} after automatic pruning attempts; dropping this cluster. Check ${OUTDIR}/logs/trycycler_reconcile_${cname}*.log"
                fi
            done
            (( ${#cluster_dirs[@]} > 0 )) || \
                err "All clusters failed reconciliation; Trycycler cannot produce a consensus for this sample. Check ${OUTDIR}/logs/trycycler_reconcile_*.log"

            log "Trycycler [5/6]: multiple sequence alignment + read partitioning across ${#cluster_dirs[@]} surviving cluster(s)..."
            for cdir in "${cluster_dirs[@]}"; do
                cname=$(basename "${cdir}")
                trycycler msa -c "${cdir}" -t "${THREADS}" \
                    > "${OUTDIR}/logs/trycycler_msa_${cname}.log" 2>&1 || \
                    err "Trycycler msa failed for ${cname}. Check ${OUTDIR}/logs/trycycler_msa_${cname}.log"
            done
            trycycler partition -r "${FILTERED_LONG}" -c "${cluster_dirs[@]}" -t "${THREADS}" \
                > "${OUTDIR}/logs/trycycler_partition.log" 2>&1 || \
                err "Trycycler partition failed. Check ${OUTDIR}/logs/trycycler_partition.log"

            log "Trycycler [6/6]: generating consensus per cluster..."
            : > "${tc_dir}/consensus_combined.fasta"
            for cdir in "${cluster_dirs[@]}"; do
                cname=$(basename "${cdir}")
                trycycler consensus -c "${cdir}" -t "${THREADS}" \
                    > "${OUTDIR}/logs/trycycler_consensus_${cname}.log" 2>&1
                if [[ -s "${cdir}/7_final_consensus.fasta" ]]; then
                    cat "${cdir}/7_final_consensus.fasta" >> "${tc_dir}/consensus_combined.fasta"
                else
                    warn "Trycycler consensus failed for ${cname}; its replicon will be missing from the final assembly. Check ${OUTDIR}/logs/trycycler_consensus_${cname}.log"
                fi
            done

            ASSEMBLY="${tc_dir}/consensus_combined.fasta"
            ;;
    esac

    [[ -f "${ASSEMBLY}" && -s "${ASSEMBLY}" ]] || \
        err "Assembly output not found or empty. Check ${OUTDIR}/logs/${ASSEMBLER}.log"

    # ── Contig filtering ─────────────────────────────────────────────────────
    log "Filtering contigs shorter than ${MIN_CONTIG_LENGTH} bp..."
    local filtered="${ASSEMBLY%.fasta}_filtered.fasta"

    if command -v seqkit &>/dev/null; then
        seqkit seq --min-len "${MIN_CONTIG_LENGTH}" --quiet "${ASSEMBLY}" \
            > "${filtered}"
    else
        warn "seqkit not found; using awk fallback for contig filtering."
        awk -v min="${MIN_CONTIG_LENGTH}" '
            /^>/ {
                if (hdr != "" && length(seq) >= min) print hdr "\n" seq
                hdr = $0; seq = ""
                next
            }
            { seq = seq $0 }
            END { if (hdr != "" && length(seq) >= min) print hdr "\n" seq }
        ' "${ASSEMBLY}" > "${filtered}"
    fi

    ASSEMBLY="${filtered}"
    local ctg_count grep_rc=0
    # Distinguish a real zero-contig count (grep_rc=1) from an I/O error (>1).
    ctg_count=$(grep -c '^>' "${ASSEMBLY}") || grep_rc=$?
    (( grep_rc > 1 )) && \
        err "Could not count contigs in ${ASSEMBLY} (grep exited ${grep_rc}) — check the file exists and is readable."
    [[ "${ctg_count}" -eq 0 ]] && \
        err "No contigs remain after filtering to ≥ ${MIN_CONTIG_LENGTH} bp. Assembly failed, or -c (min contig length) is too strict for this dataset."
    ok "Assembly complete: ${ctg_count} contigs ≥ ${MIN_CONTIG_LENGTH} bp"

    echo "${ASSEMBLY}" > "${OUTDIR}/logs/.assembly_path"
    stage_done "assembly"
}

# Restore ASSEMBLY variable when resuming past the assembly stage.
_restore_assembly_path() {
    local path_file="${OUTDIR}/logs/.assembly_path"
    [[ -f "${path_file}" ]] || err "Cannot resume: assembly path file missing. Re-run without -r."
    ASSEMBLY=$(cat "${path_file}")
    [[ -f "${ASSEMBLY}" && -s "${ASSEMBLY}" ]] || \
        err "Cannot resume: assembly file missing or empty: ${ASSEMBLY}"
    log "Resumed with assembly: ${ASSEMBLY}"
}

# ── Polishing ──────────────────────────────────────────────────────────────────
run_polishing() {
    skip_if_done "polishing" && {
        _restore_assembly_path
        return
    }

    log "[4/5] Polishing assembly..."

    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then

        local avail_ram
        avail_ram=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo 2>/dev/null || echo 0)
        local MEM_INT
        MEM_INT=$(memory_to_gb_int "${MEMORY}") || err "Cannot parse memory limit '${MEMORY}'."
        if (( avail_ram > 0 && MEM_INT > avail_ram )); then
            warn "Pilon memory limit (${MEMORY}) exceeds available RAM (${avail_ram}G). " \
                 "Consider reducing -m or closing other processes."
        fi

        log "Mapping Illumina reads for Pilon polishing..."
        bowtie2-build --threads "${THREADS}" --quiet \
            "${ASSEMBLY}" "${OUTDIR}/03_polishing/bt2_idx" \
            > "${OUTDIR}/logs/bowtie2_build.log" 2>&1

        bowtie2 \
            -x "${OUTDIR}/03_polishing/bt2_idx" \
            -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}" \
            --threads "${THREADS}" \
            --very-sensitive-local \
            --no-unal \
            2> "${OUTDIR}/logs/bowtie2.log" \
        | samtools sort -@ "${THREADS}" \
            -o "${OUTDIR}/03_polishing/mapped.bam"

        samtools index -@ "${THREADS}" "${OUTDIR}/03_polishing/mapped.bam"

        log "Running Pilon..."
        pilon \
            -Xmx"${MEM_INT}"g \
            --genome "${ASSEMBLY}" \
            --frags "${OUTDIR}/03_polishing/mapped.bam" \
            --output "${SAMPLE}_pilon" \
            --outdir "${OUTDIR}/03_polishing/" \
            --threads "${THREADS}" \
            --changes --vcf \
            --fix all \
            > "${OUTDIR}/logs/pilon.log" 2>&1

        ASSEMBLY="${OUTDIR}/03_polishing/${SAMPLE}_pilon.fasta"
        ok "Pilon polishing done"
    elif [[ "${READ_TYPE}" == "long" && "${HIFI_MODE}" == false ]]; then
        log "Running Medaka polishing (Nanopore consensus)..."
        local medaka_model_args=(--bacteria)
        [[ -n "${MEDAKA_MODEL}" ]] && medaka_model_args=(-m "${MEDAKA_MODEL}")
        medaka_consensus \
            -i "${FILTERED_LONG}" \
            -d "${ASSEMBLY}" \
            -o "${OUTDIR}/03_polishing/medaka" \
            -t "${THREADS}" \
            "${medaka_model_args[@]}" \
            > "${OUTDIR}/logs/medaka.log" 2>&1 || true

        if [[ -s "${OUTDIR}/03_polishing/medaka/consensus.fasta" ]]; then
            ASSEMBLY="${OUTDIR}/03_polishing/medaka/consensus.fasta"
            ok "Medaka polishing done"
        else
            warn "Medaka did not produce consensus.fasta; keeping unpolished assembly. Check ${OUTDIR}/logs/medaka.log — older/legacy reads often need an explicit model via -M (e.g. -M r941_min_sup_g507), since --bacteria auto-selection requires basecaller metadata in the read headers."
        fi
    else
        log "PacBio HiFi assembly — skipping external polishing (assembler handles consensus)."
    fi

    ok "Polishing complete: ${ASSEMBLY}"
    echo "${ASSEMBLY}" > "${OUTDIR}/logs/.assembly_path"
    stage_done "polishing"
}

# ── Quality assessment ─────────────────────────────────────────────────────────
run_quality_assessment() {
    skip_if_done "qa" && return
    log "[5/5] Quality assessment..."

    local _gs_bases
    _gs_bases=$(genome_size_to_bases "${GENOME_SIZE}") || \
        err "Cannot parse genome size '${GENOME_SIZE}' for QUAST --est-ref-size."

    quast "${ASSEMBLY}" \
        -o "${OUTDIR}/reports/quast/" \
        -t "${THREADS}" \
        --min-contig "${MIN_CONTIG_LENGTH}" \
        --gene-finding \
        --rna-finding \
        --est-ref-size "${_gs_bases}" \
        > "${OUTDIR}/logs/quast.log" 2>&1

    if command -v checkm2 &>/dev/null; then
        log "Running CheckM2..."
        mkdir -p "${OUTDIR}/reports/checkm2"
        checkm2 predict \
            --input "${ASSEMBLY}" \
            --output-directory "${OUTDIR}/reports/checkm2/" \
            --threads "${THREADS}" \
            --force \
            > "${OUTDIR}/logs/checkm2.log" 2>&1
    elif command -v checkm &>/dev/null; then
        log "Running CheckM (legacy)..."
        local checkm_input="${OUTDIR}/reports/checkm/input"
        mkdir -p "${checkm_input}"
        cp "${ASSEMBLY}" "${checkm_input}/${SAMPLE}.fasta"
        checkm lineage_wf \
            -x fasta -t "${THREADS}" \
            "${checkm_input}" \
            "${OUTDIR}/reports/checkm/" \
            > "${OUTDIR}/logs/checkm.log" 2>&1
        checkm qa \
            "${OUTDIR}/reports/checkm/lineage.ms" \
            "${OUTDIR}/reports/checkm/" \
            -o 2 --tab_table \
            > "${OUTDIR}/reports/checkm/checkm_results.tsv"
    else
        warn "Neither checkm2 nor checkm found; skipping completeness assessment."
    fi

    if command -v busco &>/dev/null; then
        local busco_lineage_args=(--auto-lineage-prok)
        [[ -n "${BUSCO_LINEAGE}" ]] && busco_lineage_args=(-l "${BUSCO_LINEAGE}")

        log "Running BUSCO (lineage: ${BUSCO_LINEAGE:-auto-detect via --auto-lineage-prok})..."
        if ! busco \
            -i "${ASSEMBLY}" \
            -o "${SAMPLE}_busco" \
            "${busco_lineage_args[@]}" \
            -m genome \
            -c "${THREADS}" \
            --out_path "${OUTDIR}/reports/busco/" \
            --quiet \
            > "${OUTDIR}/logs/busco.log" 2>&1
        then
            if [[ -z "${BUSCO_LINEAGE}" ]]; then
                warn "BUSCO auto-lineage placement failed; retrying with generic bacteria_odb10. Check ${OUTDIR}/logs/busco.log — pass -L <lineage_odb10> to force a specific lineage instead (e.g. -L enterobacterales_odb10)."
                rm -rf "${OUTDIR}/reports/busco/${SAMPLE}_busco"
                busco \
                    -i "${ASSEMBLY}" \
                    -o "${SAMPLE}_busco" \
                    -l bacteria_odb10 \
                    -m genome \
                    -c "${THREADS}" \
                    --out_path "${OUTDIR}/reports/busco/" \
                    --quiet \
                    >> "${OUTDIR}/logs/busco.log" 2>&1 || warn "BUSCO fallback run also failed; skipping completeness-by-lineage assessment. Check ${OUTDIR}/logs/busco.log."
            else
                warn "BUSCO failed with lineage ${BUSCO_LINEAGE}. Check ${OUTDIR}/logs/busco.log."
            fi
        fi
    else
        warn "BUSCO not found; skipping."
    fi

    # Opt-in typing/resistance analyses; tool presence already checked in validate_inputs().
    if ${PLASMID_DETECTION}; then
        # -p needs the DB directory (auto-detected below if -D not given); -d takes the DB name.
        local pf_db_path="${PLASMID_DB_PATH}"
        if [[ -z "${pf_db_path}" ]]; then
            pf_db_path=$(find "${CONDA_PREFIX:-/usr}/share" -maxdepth 2 -type d \
                -path '*/plasmidfinder-*/database' 2>/dev/null | head -1)
        fi
        [[ -n "${pf_db_path}" && -f "${pf_db_path}/config" ]] || \
            err "PlasmidFinder database not found (looked under \${CONDA_PREFIX}/share/plasmidfinder-*/database). Pass -D <path>, or run the tool's download-db.sh."

        log "Running PlasmidFinder (db: ${PLASMID_DB}, path: ${pf_db_path})..."
        plasmidfinder.py \
            -i "${ASSEMBLY}" \
            -o "${OUTDIR}/reports/plasmidfinder/" \
            -p "${pf_db_path}" \
            -d "${PLASMID_DB}" \
            > "${OUTDIR}/logs/plasmidfinder.log" 2>&1
        ok "PlasmidFinder complete"
    fi

    if ${AMR_DETECTION}; then
        log "Running AMRFinderPlus (--plus panel: AMR + stress/virulence genes)..."
        amrfinder \
            -n "${ASSEMBLY}" \
            -o "${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv" \
            --threads "${THREADS}" \
            --plus \
            > "${OUTDIR}/logs/amrfinder.log" 2>&1
        ok "AMRFinderPlus complete"
    fi

    if ${ABRICATE_DETECTION}; then
        log "Running abricate (db: ${ABRICATE_DB})..."
        abricate \
            --db "${ABRICATE_DB}" \
            --threads "${THREADS}" \
            "${ASSEMBLY}" \
            > "${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv" \
            2> "${OUTDIR}/logs/abricate.log"
        ok "abricate complete"
    fi

    if ${MLST_DETECTION}; then
        log "Running mlst..."
        mlst \
            --threads "${THREADS}" \
            "${ASSEMBLY}" \
            > "${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv" \
            2> "${OUTDIR}/logs/mlst.log"
        ok "mlst complete"
    fi

    ok "Quality assessment complete"
    stage_done "qa"
}

# ── Cleanup ────────────────────────────────────────────────────────────────────
# Removes heavy intermediates after a full run, keeping only reports/logs/genome. -k disables this.
run_cleanup() {
    log "Cleaning up intermediate files (use -k to keep them)..."

    local final_assembly="${OUTDIR}/$(_input_basename).fasta"
    mv "${ASSEMBLY}" "${final_assembly}"
    ASSEMBLY="${final_assembly}"
    echo "${ASSEMBLY}" > "${OUTDIR}/logs/.assembly_path"

    rm -rf "${OUTDIR}/00_raw_data" "${OUTDIR}/02_assembly" "${OUTDIR}/03_polishing"
    rm -f "${TRIMMED_R1:-}" "${TRIMMED_R2:-}" "${FILTERED_LONG:-}"
    # Force qc/assembly/polishing/qa to redo on a later -r resume of this cleaned OUTDIR.
    rm -f "${OUTDIR}/logs/.done_qc" "${OUTDIR}/logs/.done_assembly" \
          "${OUTDIR}/logs/.done_polishing" "${OUTDIR}/logs/.done_qa"

    ok "Cleanup complete — kept qc/, reports/, logs/, and ${final_assembly}"
}

# ── Final report ───────────────────────────────────────────────────────────────
generate_report() {
    log "Generating final report..."

    local n_contigs total_len n50 largest
    n_contigs=$(grep -c '^>' "${ASSEMBLY}" 2>/dev/null || echo "N/A")

    if command -v assembly-stats &>/dev/null; then
        local _stats
        _stats=$(assembly-stats "${ASSEMBLY}" 2>/dev/null)
        # Fields are comma-terminated (e.g. "sum = 123, n = 1"); strip the comma before matching.
        total_len=$(awk '/^sum[ =]|sum[_ ]len|total[_ ]len/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}' <<< "${_stats}")
        n50=$(      awk '/^N50[ =]/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}'         <<< "${_stats}")
        # largest shares a line with sum; scan forward from the "largest" field, not from line start.
        largest=$(  awk '{for(i=1;i<=NF;i++){if($i~/largest/){for(j=i;j<=NF;j++){v=$j; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}}}' <<< "${_stats}")
        total_len="${total_len:-N/A}"
        n50="${n50:-N/A}"
        largest="${largest:-N/A}"
    else
        total_len="N/A (install assembly-stats)"
        n50="N/A"; largest="N/A"
    fi

    # Pull one headline number out of each opt-in typing/resistance report, if it ran.
    local mlst_st="N/A" amr_gene_count="N/A" abricate_gene_count="N/A" pf_hit_count="N/A"
    if ${MLST_DETECTION} && [[ -s "${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv" ]]; then
        mlst_st=$(awk -F'\t' 'NR==1{print $3}' "${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv")
    fi
    if ${AMR_DETECTION} && [[ -s "${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv" ]]; then
        amr_gene_count=$(( $(wc -l < "${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv") - 1 ))
    fi
    if ${ABRICATE_DETECTION} && [[ -s "${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv" ]]; then
        abricate_gene_count=$(( $(wc -l < "${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv") - 1 ))
    fi
    if ${PLASMID_DETECTION}; then
        local pf_dir="${OUTDIR}/reports/plasmidfinder"
        local pf_tab
        pf_tab=$(find "${pf_dir}" -maxdepth 1 -iname 'results_tab.tsv' 2>/dev/null | head -1)
        if [[ -n "${pf_tab}" && -s "${pf_tab}" ]]; then
            pf_hit_count=$(( $(wc -l < "${pf_tab}") - 1 ))
        elif [[ -s "${pf_dir}/data.json" ]] && grep -q '"No hit found"' "${pf_dir}/data.json"; then
            # PlasmidFinder only writes results_tab.tsv when there's at least one hit;
            # a confirmed-zero result leaves just data.json.
            pf_hit_count=0
        fi
    fi

    local species_report_note=""
    [[ -f "${OUTDIR}/reports/kraken2/${SAMPLE}_kraken2_report.tsv" ]] && \
        species_report_note=" (Kraken2-detected — review \`${OUTDIR}/reports/kraken2/${SAMPLE}_kraken2_report.tsv\` if unexpected)"

    # Empty unless at least one --plasmid/--amr/--abricate/--mlst flag was passed.
    # Built in both Markdown (SUMMARY.md) and HTML (SUMMARY.html) forms together
    # from the same values, so the two reports can't drift out of sync.
    local typing_section="" typing_section_html=""
    if ${MLST_DETECTION} || ${AMR_DETECTION} || ${ABRICATE_DETECTION} || ${PLASMID_DETECTION}; then
        typing_section+="## Typing & Resistance"$'\n\n'
        typing_section_html="<h2>Typing &amp; Resistance</h2>"$'\n'"<ul>"$'\n'
        if ${MLST_DETECTION}; then
            typing_section+="- MLST sequence type: **${mlst_st}** (\`${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv\`)"$'\n'
            typing_section_html+="<li>MLST sequence type: <strong>${mlst_st}</strong> (<code>${OUTDIR}/reports/mlst/${SAMPLE}_mlst.tsv</code>)</li>"$'\n'
        fi
        if ${AMR_DETECTION}; then
            typing_section+="- AMRFinderPlus: **${amr_gene_count}** gene(s)/mutation(s) (\`${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv\`)"$'\n'
            typing_section_html+="<li>AMRFinderPlus: <strong>${amr_gene_count}</strong> gene(s)/mutation(s) (<code>${OUTDIR}/reports/amrfinderplus/${SAMPLE}_amrfinder.tsv</code>)</li>"$'\n'
        fi
        if ${ABRICATE_DETECTION}; then
            typing_section+="- abricate (${ABRICATE_DB}): **${abricate_gene_count}** hit(s) (\`${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv\`)"$'\n'
            typing_section_html+="<li>abricate (${ABRICATE_DB}): <strong>${abricate_gene_count}</strong> hit(s) (<code>${OUTDIR}/reports/abricate/${SAMPLE}_abricate_${ABRICATE_DB}.tsv</code>)</li>"$'\n'
        fi
        if ${PLASMID_DETECTION}; then
            typing_section+="- PlasmidFinder: **${pf_hit_count}** replicon hit(s) (\`${OUTDIR}/reports/plasmidfinder/\`)"$'\n'
            typing_section_html+="<li>PlasmidFinder: <strong>${pf_hit_count}</strong> replicon hit(s) (<code>${OUTDIR}/reports/plasmidfinder/</code>)</li>"$'\n'
        fi
        typing_section+=$'\n'
        typing_section_html+="</ul>"$'\n'
    fi

    local input_files_html=""
    [[ -n "${R1:-}"   ]] && input_files_html+="<li>Illumina R1: <code>${R1}</code></li>"$'\n'
    [[ -n "${R2:-}"   ]] && input_files_html+="<li>Illumina R2: <code>${R2}</code></li>"$'\n'
    [[ -n "${LONG:-}" ]] && input_files_html+="<li>Long reads: <code>${LONG}</code></li>"$'\n'

    cat > "${OUTDIR}/reports/SUMMARY.md" <<EOF
# TakiLine Assembly Report

| Field           | Value                       |
|:----------------|:----------------------------|
| Sample          | ${SAMPLE}                   |
| Species         | ${SPECIES:-not determined}${species_report_note} |
| Date            | $(date '+%Y-%m-%d %H:%M')   |
| Read type       | ${READ_TYPE}                |
| Assembler       | ${ASSEMBLER}                |
| HiFi mode       | ${HIFI_MODE}                |

## Input Files
$(  [[ -n "${R1:-}"   ]] && echo "- Illumina R1: \`${R1}\`")
$(  [[ -n "${R2:-}"   ]] && echo "- Illumina R2: \`${R2}\`")
$(  [[ -n "${LONG:-}" ]] && echo "- Long reads:  \`${LONG}\`")

## Assembly Statistics

| Metric          | Value        |
|:----------------|:-------------|
| Contigs         | ${n_contigs} |
| Total length    | ${total_len} |
| Largest contig  | ${largest}   |
| N50             | ${n50}       |
| Min contig len  | ${MIN_CONTIG_LENGTH} bp |

${typing_section}## Output Files

| Step            | Path                                        |
|:----------------|:--------------------------------------------|
| Final assembly  | \`${ASSEMBLY}\`                             |
| QC reports      | \`${OUTDIR}/qc/\`                        |
| Quality reports | \`${OUTDIR}/reports/\`                   |
| Logs            | \`${OUTDIR}/logs/\`                         |

## Suggested Next Steps
1. **Annotation**     – \`Bakta\` (INSDC-ready) or \`PGAP\` for NCBI submission
2. **AMR screening**  – built in via \`--amr\`/\`--abricate\` (see above if run); or \`ResFinder\` for a third opinion
3. **MLST typing**    – built in via \`--mlst\` (see above if run)
4. **Phylogenetics**  – \`IQ-TREE2\` or \`FastTree\`
5. **Pan-genome**     – \`Panaroo\` (preferred) or \`Roary\`
EOF

    cat > "${OUTDIR}/reports/SUMMARY.html" <<EOF
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TakiLine Report — ${SAMPLE}</title>
<style>
  body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; max-width: 900px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; line-height: 1.5; }
  h1 { border-bottom: 3px solid #2c7a4b; padding-bottom: .3rem; }
  h2 { border-bottom: 1px solid #ddd; padding-bottom: .2rem; margin-top: 2rem; }
  table { border-collapse: collapse; width: 100%; margin: 1rem 0; }
  th, td { border: 1px solid #ddd; padding: .5rem .75rem; text-align: left; }
  th { background: #2c7a4b; color: #fff; }
  tr:nth-child(even) { background: #f6f6f6; }
  code { background: #eee; padding: .1rem .3rem; border-radius: 3px; font-size: .9em; }
  ul { padding-left: 1.3rem; }
</style>
</head>
<body>
<h1>TakiLine Assembly Report</h1>
<table>
<tr><th>Field</th><th>Value</th></tr>
<tr><td>Sample</td><td>${SAMPLE}</td></tr>
<tr><td>Species</td><td>${SPECIES:-not determined}${species_report_note}</td></tr>
<tr><td>Date</td><td>$(date '+%Y-%m-%d %H:%M')</td></tr>
<tr><td>Read type</td><td>${READ_TYPE}</td></tr>
<tr><td>Assembler</td><td>${ASSEMBLER}</td></tr>
<tr><td>HiFi mode</td><td>${HIFI_MODE}</td></tr>
</table>

<h2>Input Files</h2>
<ul>
${input_files_html}</ul>

<h2>Assembly Statistics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Contigs</td><td>${n_contigs}</td></tr>
<tr><td>Total length</td><td>${total_len}</td></tr>
<tr><td>Largest contig</td><td>${largest}</td></tr>
<tr><td>N50</td><td>${n50}</td></tr>
<tr><td>Min contig len</td><td>${MIN_CONTIG_LENGTH} bp</td></tr>
</table>

${typing_section_html}
<h2>Output Files</h2>
<table>
<tr><th>Step</th><th>Path</th></tr>
<tr><td>Final assembly</td><td><code>${ASSEMBLY}</code></td></tr>
<tr><td>QC reports</td><td><code>${OUTDIR}/qc/</code></td></tr>
<tr><td>Quality reports</td><td><code>${OUTDIR}/reports/</code></td></tr>
<tr><td>Logs</td><td><code>${OUTDIR}/logs/</code></td></tr>
</table>

<h2>Suggested Next Steps</h2>
<ol>
<li><strong>Annotation</strong> – <code>Bakta</code> (INSDC-ready) or <code>PGAP</code> for NCBI submission</li>
<li><strong>AMR screening</strong> – built in via <code>--amr</code>/<code>--abricate</code> (see above if run); or <code>ResFinder</code> for a third opinion</li>
<li><strong>MLST typing</strong> – built in via <code>--mlst</code> (see above if run)</li>
<li><strong>Phylogenetics</strong> – <code>IQ-TREE2</code> or <code>FastTree</code></li>
<li><strong>Pan-genome</strong> – <code>Panaroo</code> (preferred) or <code>Roary</code></li>
</ol>
</body>
</html>
EOF

    echo -e "\n${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}       PIPELINE COMPLETE ✓            ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "  Sample    : ${SAMPLE}"
    echo -e "  Assembly  : ${ASSEMBLY}"
    echo -e "  Contigs   : ${n_contigs}"
    echo -e "  N50       : ${n50}"
    echo -e "  Report    : ${OUTDIR}/reports/SUMMARY.md (${OUTDIR}/reports/SUMMARY.html)"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"

    # Disabled last, after every command that could still fail here.
    trap - ERR
}

# ── Error trap ─────────────────────────────────────────────────────────────────
trap 'err "Pipeline failed at line ${LINENO}. Check logs in ${OUTDIR}/logs/"' ERR

# ── Main ───────────────────────────────────────────────────────────────────────
main() {
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}               TakiLine               ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"

    validate_inputs
    run_qc
    run_species_id

    if ${QC_ONLY}; then
        ok "QC-only mode. Exiting."
        trap - ERR
        exit 0
    fi

    run_assembly
    run_polishing
    run_quality_assessment
    ${CLEANUP} && run_cleanup
    generate_report
}

main "$@"
