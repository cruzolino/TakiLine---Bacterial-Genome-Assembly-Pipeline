#!/usr/bin/env bash
# ==============================================================================
# Bacterial Genome Assembly Pipeline — v6.1.1
# Modes: Illumina (PE) | Nanopore | Hybrid (Illumina + Nanopore) | PacBio HiFi
#
# Changes from v6.1.0 (bug-fix release):
#   - BUG FIX (critical): `set -o errtrace` added — the ERR trap previously
#     never fired inside functions (i.e. during the entire pipeline), so
#     failures exited silently with no "Pipeline failed at line X" message.
#   - BUG FIX: Flye version detection no longer crashes the whole run when
#     `grep -oP` is unavailable/fails (non-GNU grep, non-UTF8 locale); now
#     uses portable `grep -E` under a forced C locale and never aborts the
#     script on failure to detect a version.
#   - BUG FIX: SUMMARY.md assembly-stats parsing — real assembly-stats output
#     terminates numeric fields with a comma (e.g. "N50 = 4655090, n = 1"),
#     which the old strict `^[0-9]+$` match never matched, silently picking up
#     the wrong field. Total length was never populated at all because the
#     regex looked for "sum_len"/"total_len" instead of the real "sum" field.
#   - BUG FIX: genome-size validation regex tightened to reject malformed
#     values (e.g. "4.8.2m"), which previously passed validation and silently
#     corrupted the Filtlong TARGET_BASES calculation.
#   - BUG FIX: hybrid-mode assembler validation — spades/skesa now hard-error
#     in hybrid mode instead of silently discarding the long reads (they were
#     QC'd but never assembled with); flye/canu/raven now warn that Illumina
#     reads are used for Pilon polishing only, not integrated hybrid assembly.
#   - Default output directory renamed from `bacterial_assembly` to `takiline`.
#   - All v6.1.0 changes retained.
# ==============================================================================
set -Eeuo pipefail
IFS=$'\n\t'

# ── Version ────────────────────────────────────────────────────────────────────
readonly VERSION="6.1.1"

# ── Defaults ───────────────────────────────────────────────────────────────────
THREADS=$(nproc --all 2>/dev/null || echo 8)
MEMORY="32G"
OUTDIR="takiline"
ASSEMBLER="spades"
READ_TYPE="illumina"
MIN_CONTIG_LENGTH=500
GENOME_SIZE="5m"
SAMPLE="isolate"
PLASMID_DETECTION=true
PLASMID_DB="enterobacteriaceae"
QC_ONLY=false
HIFI_MODE=false
RESUME=false

# ── Colors ─────────────────────────────────────────────────────────────────────
readonly RED='\033[0;31m' GREEN='\033[0;32m' YELLOW='\033[1;33m' \
         BLUE='\033[0;34m' CYAN='\033[0;36m' NC='\033[0m'

# ── Logging ────────────────────────────────────────────────────────────────────
log()  { echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $*"; }
ok()   { echo -e "${GREEN}[✓ $(date +%H:%M:%S)]${NC} $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
err()  { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

# ── Derived-path helpers ───────────────────────────────────────────────────────
# Paths for QC outputs are deterministic; centralising them here fixes resume:
# when a run resumes past the QC stage, run_assembly/run_polishing need these
# variables, but run_qc was never called to set them.
_set_derived_paths() {
    TRIMMED_R1="${OUTDIR}/01_qc/${SAMPLE}_R1.trimmed.fastq.gz"
    TRIMMED_R2="${OUTDIR}/01_qc/${SAMPLE}_R2.trimmed.fastq.gz"
    FILTERED_LONG="${OUTDIR}/01_qc/${SAMPLE}_filtered.fastq.gz"
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

# ── Parallel gzip ──────────────────────────────────────────────────────────────
# Use pigz (parallel gzip) when available; transparently fall back to gzip.
# Accepts the same flags as gzip.
_gzip() {
    if command -v pigz &>/dev/null; then
        pigz -p "${THREADS}" "$@"
    else
        gzip "$@"
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
${GREEN}Bacterial Genome Assembly Pipeline v${VERSION}${NC}

${YELLOW}Usage:${NC}  $0 [OPTIONS]

${YELLOW}Input (at least one required):${NC}
  -1 FILE   Illumina forward reads (R1.fastq.gz)
  -2 FILE   Illumina reverse reads (R2.fastq.gz)
  -l FILE   Long reads – Nanopore or PacBio HiFi (.fastq[.gz])
  --hifi    Treat -l reads as PacBio HiFi (activates flye --pacbio-hifi)

${YELLOW}Assembly:${NC}
  -a STR    Assembler: spades|skesa|unicycler (Illumina/hybrid)
                       flye|canu|raven (long-read)        [default: spades]

${YELLOW}General:${NC}
  -o DIR    Output directory          [default: ${OUTDIR}]
  -t INT    Threads                   [default: auto = ${THREADS}]
  -m STR    Memory limit              [default: ${MEMORY}]
  -s STR    Sample name               [default: ${SAMPLE}]
  -g STR    Expected genome size      [default: ${GENOME_SIZE}]
  -c INT    Min contig length (bp)    [default: ${MIN_CONTIG_LENGTH}]
  -P STR    PlasmidFinder DB          [default: ${PLASMID_DB}]
  -x        Skip plasmid detection
  -q        QC only (skip assembly)
  -r        Resume from last checkpoint
  -h|--help Show this help

${YELLOW}Examples:${NC}
  # Illumina PE
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

  # Nanopore
  $0 -l ont.fq.gz -s Salmonella -a flye -g 4.8m -t 16

  # Hybrid
  $0 -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz -s Klebsiella -a unicycler

  # PacBio HiFi
  $0 -l hifi.fastq.gz -s Pseudomonas -a flye --hifi -g 6.5m

  # Resume interrupted run
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
EOF
    exit 0
}

# ── Argument parsing ───────────────────────────────────────────────────────────
ARGS=()
for arg in "$@"; do
    case "${arg}" in
        --hifi)  HIFI_MODE=true ;;
        --help)  usage ;;
        *)       ARGS+=("${arg}") ;;
    esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

while getopts "1:2:l:o:t:s:g:a:m:c:P:xqrh" opt; do
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
        P) PLASMID_DB="${OPTARG}" ;;
        x) PLASMID_DETECTION=false ;;
        q) QC_ONLY=true ;;
        r) RESUME=true ;;
        h) usage ;;
        *) err "Invalid option: -${OPTARG}" ;;
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
            # spades/skesa never consume long reads: in hybrid mode they would
            # be QC'd (and required as a dependency) but silently never used
            # in assembly or polishing. Only unicycler performs true hybrid
            # assembly, so block the misleading combination outright.
            if [[ "${READ_TYPE}" == "hybrid" && "${ASSEMBLER}" != "unicycler" ]]; then
                err "${ASSEMBLER} does not use long reads — the -l input would be silently discarded in hybrid mode. Use -a unicycler for hybrid assembly, or drop -l to run ${ASSEMBLER} on Illumina reads only."
            fi
            ;;
        flye|canu|raven)
            [[ "${READ_TYPE}" == "illumina" ]] && \
                err "${ASSEMBLER} is a long-read assembler. For Illumina use spades, skesa, or unicycler."
            # flye/canu/raven assemble only the long reads; in hybrid mode the
            # Illumina reads are used later for Pilon polishing, not for an
            # integrated hybrid assembly graph. Not an error, but worth flagging.
            [[ "${READ_TYPE}" == "hybrid" ]] && \
                warn "${ASSEMBLER} assembles only the long reads; Illumina reads will be used for post-assembly Pilon polishing, not integrated hybrid assembly. Use -a unicycler for true hybrid assembly."
            ;;
        *) err "Unknown assembler: ${ASSEMBLER}" ;;
    esac

    [[ "${THREADS}" =~ ^[0-9]+$ ]] || err "Threads must be a positive integer."
    [[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]] || err "Min contig length must be a positive integer."
    [[ "${GENOME_SIZE,,}" =~ ^[0-9]+(\.[0-9]+)?[mg]$ ]] || \
        err "Genome size must be a number followed by m or g (e.g. 5m, 4.8m, 0.5g)."

    # ── Preflight tool checks (all at once) ───────────────────────────────────
    local required_tools=(fastqc fastp quast)
    [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]] && \
        required_tools+=(filtlong)
    case "${ASSEMBLER}" in
        spades)    required_tools+=(spades.py) ;;
        skesa)     required_tools+=(skesa) ;;
        unicycler) required_tools+=(unicycler) ;;
        flye)      required_tools+=(flye) ;;
        canu)      required_tools+=(canu) ;;
        raven)     required_tools+=(raven) ;;
    esac
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then
        required_tools+=(bowtie2 samtools pilon)
    fi
    check_tools "${required_tools[@]}"

    # ── Directory scaffold ────────────────────────────────────────────────────
    mkdir -p "${OUTDIR}"/{00_raw_data,01_qc/post_trim,02_assembly,\
03_polishing,05_reports/{quast,checkm,busco,plasmidfinder},logs}

    [[ -n "${R1:-}"   ]] && ln -sf "$(realpath "${R1}")"   "${OUTDIR}/00_raw_data/$(basename "${R1}")"
    [[ -n "${R2:-}"   ]] && ln -sf "$(realpath "${R2}")"   "${OUTDIR}/00_raw_data/$(basename "${R2}")"
    [[ -n "${LONG:-}" ]] && ln -sf "$(realpath "${LONG}")" "${OUTDIR}/00_raw_data/$(basename "${LONG}")"

    # Set derived paths now so they are available to all downstream stages,
    # including when those stages are resumed and run_qc is never called.
    _set_derived_paths

    ok "Input validation passed"
}

# ── QC ─────────────────────────────────────────────────────────────────────────
run_qc() {
    skip_if_done "qc" && return

    log "[1/4] Quality control..."

    # ── Illumina QC ──────────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then

        # Start FastQC on raw reads in background while fastp runs in parallel.
        log "FastQC on raw reads (background) + fastp trimming (foreground)..."
        fastqc "${R1}" "${R2}" \
            --outdir "${OUTDIR}/01_qc/" \
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
            --html "${OUTDIR}/01_qc/${SAMPLE}_fastp.html" \
            --json "${OUTDIR}/01_qc/${SAMPLE}_fastp.json" \
            --report_title "${SAMPLE}" \
            > "${OUTDIR}/logs/fastp.log" 2>&1

        # Wait for raw FastQC before running FastQC on trimmed reads
        wait "${fastqc_raw_pid}" || warn "FastQC (raw) finished with warnings; check fastqc_raw.log."

        log "FastQC on trimmed reads..."
        fastqc "${TRIMMED_R1}" "${TRIMMED_R2}" \
            --outdir "${OUTDIR}/01_qc/post_trim/" \
            --threads "${THREADS}" \
            --quiet \
            > "${OUTDIR}/logs/fastqc_post.log" 2>&1

        if command -v multiqc &>/dev/null; then
            multiqc "${OUTDIR}/01_qc/" \
                --outdir "${OUTDIR}/05_reports/" \
                --filename "${SAMPLE}_multiqc" \
                --no-megaqc-update \
                --quiet \
                > "${OUTDIR}/logs/multiqc.log" 2>&1
        else
            warn "multiqc not found; skipping MultiQC report."
        fi
    fi

    # ── Long-read QC ──────────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]]; then

        # Compute target bases for 20× coverage before launching background jobs.
        # Each branch captures independently to avoid BASH_REMATCH bleed-over.
        local _gs="${GENOME_SIZE,,}"
        local TARGET_BASES
        if [[ "${_gs}" =~ ^([0-9]+(\.[0-9]+)?)m$ ]]; then
            TARGET_BASES=$(awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]} * 1000000 * 20}")
        elif [[ "${_gs}" =~ ^([0-9]+(\.[0-9]+)?)g$ ]]; then
            TARGET_BASES=$(awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]} * 1000000000 * 20}")
        else
            err "Cannot parse genome size '${GENOME_SIZE}' for TARGET_BASES calculation."
        fi

        # Start NanoPlot in background while filtlong runs in foreground.
        local nanoplot_pid=""
        if command -v NanoPlot &>/dev/null; then
            log "NanoPlot QC (background) + Filtlong filtering (foreground)..."
            NanoPlot \
                --fastq "${LONG}" \
                --outdir "${OUTDIR}/01_qc/nanoplot/" \
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

# ── Assembly ───────────────────────────────────────────────────────────────────
run_assembly() {
    skip_if_done "assembly" && {
        _restore_assembly_path
        return
    }

    log "[2/4] Assembly with ${ASSEMBLER}..."

    local MEM_INT="${MEMORY//[^0-9]/}"

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
            # Auto-detect Flye version: --nano-hq requires Flye >= 2.9.
            # Flye <= 2.8 uses --nano-raw (or --nano-corr for corrected reads).
            local flye_read_flag
            if ${HIFI_MODE}; then
                flye_read_flag="--pacbio-hifi"
            else
                # Use POSIX `grep -E` under a forced C locale (portable across
                # GNU/BSD grep and locale configurations, unlike `grep -oP`,
                # which errors out on non-GNU grep or non-UTF8 locales and
                # previously crashed the whole run). Never let detection
                # failure abort the script; fall back to the safe default.
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
    local ctg_count
    ctg_count=$(grep -c '^>' "${ASSEMBLY}")
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
    # When resuming, the persisted .assembly_path already points to the Pilon
    # output for Illumina/hybrid runs (written after Pilon completes below).
    skip_if_done "polishing" && {
        _restore_assembly_path
        return
    }

    log "[3/4] Polishing assembly..."

    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then

        local avail_ram
        avail_ram=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo 2>/dev/null || echo 0)
        local MEM_INT="${MEMORY//[^0-9]/}"
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
    else
        # Long-read assembly: no external polishing needed.
        # ASSEMBLY is already set correctly from run_assembly(); persist it now
        # so that _restore_assembly_path works correctly on resume, and so that
        # run_quality_assessment / generate_report always receive a valid path.
        log "Long-read assembly — skipping external polishing (assembler handles consensus)."
    fi

    ok "Polishing complete: ${ASSEMBLY}"
    # Persist path in both branches (Illumina/hybrid: after Pilon; long-read: here).
    echo "${ASSEMBLY}" > "${OUTDIR}/logs/.assembly_path"
    stage_done "polishing"
}

# ── Quality assessment ─────────────────────────────────────────────────────────
run_quality_assessment() {
    skip_if_done "qa" && return
    log "[4/4] Quality assessment..."

    quast "${ASSEMBLY}" \
        -o "${OUTDIR}/05_reports/quast/" \
        -t "${THREADS}" \
        --min-contig "${MIN_CONTIG_LENGTH}" \
        --gene-finding \
        --rna-finding \
        --est-ref-size "${GENOME_SIZE}" \
        > "${OUTDIR}/logs/quast.log" 2>&1

    if command -v checkm2 &>/dev/null; then
        log "Running CheckM2..."
        mkdir -p "${OUTDIR}/05_reports/checkm2"
        checkm2 predict \
            --input "${ASSEMBLY}" \
            --output-directory "${OUTDIR}/05_reports/checkm2/" \
            --threads "${THREADS}" \
            --force \
            > "${OUTDIR}/logs/checkm2.log" 2>&1
    elif command -v checkm &>/dev/null; then
        log "Running CheckM (legacy)..."
        local checkm_input="${OUTDIR}/05_reports/checkm/input"
        mkdir -p "${checkm_input}"
        cp "${ASSEMBLY}" "${checkm_input}/${SAMPLE}.fasta"
        checkm lineage_wf \
            -x fasta -t "${THREADS}" \
            "${checkm_input}" \
            "${OUTDIR}/05_reports/checkm/" \
            > "${OUTDIR}/logs/checkm.log" 2>&1
        checkm qa \
            "${OUTDIR}/05_reports/checkm/lineage.ms" \
            "${OUTDIR}/05_reports/checkm/" \
            -o 2 --tab_table \
            > "${OUTDIR}/05_reports/checkm/checkm_results.tsv"
    else
        warn "Neither checkm2 nor checkm found; skipping completeness assessment."
    fi

    if command -v busco &>/dev/null; then
        log "Running BUSCO..."
        busco \
            -i "${ASSEMBLY}" \
            -o "${SAMPLE}_busco" \
            -l bacteria_odb10 \
            -m genome \
            -c "${THREADS}" \
            --out_path "${OUTDIR}/05_reports/busco/" \
            --quiet \
            > "${OUTDIR}/logs/busco.log" 2>&1
    else
        warn "BUSCO not found; skipping."
    fi

    if ${PLASMID_DETECTION}; then
        if command -v plasmidfinder.py &>/dev/null; then
            log "Running PlasmidFinder (db: ${PLASMID_DB})..."
            mkdir -p "${OUTDIR}/05_reports/plasmidfinder"
            plasmidfinder.py \
                -i "${ASSEMBLY}" \
                -o "${OUTDIR}/05_reports/plasmidfinder/" \
                -p "${PLASMID_DB}" \
                > "${OUTDIR}/logs/plasmidfinder.log" 2>&1
        else
            warn "plasmidfinder.py not found; skipping plasmid detection."
        fi
    fi

    ok "Quality assessment complete"
    stage_done "qa"
}

# ── Final report ───────────────────────────────────────────────────────────────
generate_report() {
    log "Generating final report..."

    local n_contigs total_len n50 largest
    n_contigs=$(grep -c '^>' "${ASSEMBLY}" 2>/dev/null || echo "N/A")

    if command -v assembly-stats &>/dev/null; then
        local _stats
        _stats=$(assembly-stats "${ASSEMBLY}" 2>/dev/null)
        # Field names vary by assembly-stats version; match the numeric value
        # after the keyword on the same line. Real assembly-stats output
        # terminates most numeric fields with a comma (e.g. "sum = 123, n = 1"),
        # so strip a trailing comma before validating each field as numeric —
        # otherwise the comma-suffixed value never matches and the loop picks
        # up the next unrelated number on the line (e.g. n50 became "1").
        total_len=$(awk '/^sum[ =]|sum[_ ]len|total[_ ]len/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}' <<< "${_stats}")
        n50=$(      awk '/^N50[ =]/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}'         <<< "${_stats}")
        largest=$(  awk '/largest/{for(i=1;i<=NF;i++){v=$i; gsub(/,$/,"",v); if(v~/^[0-9]+$/){print v; exit}}}'      <<< "${_stats}")
        total_len="${total_len:-N/A}"
        n50="${n50:-N/A}"
        largest="${largest:-N/A}"
    else
        total_len="N/A (install assembly-stats)"
        n50="N/A"; largest="N/A"
    fi

    cat > "${OUTDIR}/05_reports/SUMMARY.md" <<EOF
# Bacterial Genome Assembly Report

| Field           | Value                       |
|:----------------|:----------------------------|
| Sample          | ${SAMPLE}                   |
| Date            | $(date '+%Y-%m-%d %H:%M')   |
| Pipeline        | v${VERSION}                 |
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

## Output Files

| Step            | Path                                        |
|:----------------|:--------------------------------------------|
| Final assembly  | \`${ASSEMBLY}\`                             |
| QC reports      | \`${OUTDIR}/01_qc/\`                        |
| Quality reports | \`${OUTDIR}/05_reports/\`                   |
| Logs            | \`${OUTDIR}/logs/\`                         |

## Suggested Next Steps
1. **Annotation**     – \`Bakta\` (INSDC-ready) or \`PGAP\` for NCBI submission
2. **AMR screening**  – \`AMRFinderPlus\` or \`ResFinder\`
3. **MLST typing**    – \`mlst\` (Torsten Seemann)
4. **Phylogenetics**  – \`IQ-TREE2\` or \`FastTree\`
5. **Pan-genome**     – \`Panaroo\` (preferred) or \`Roary\`
EOF

    trap - ERR

    echo -e "\n${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}       PIPELINE COMPLETE ✓            ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "  Sample    : ${SAMPLE}"
    echo -e "  Assembly  : ${ASSEMBLY}"
    echo -e "  Contigs   : ${n_contigs}"
    echo -e "  N50       : ${n50}"
    echo -e "  Report    : ${OUTDIR}/05_reports/SUMMARY.md"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"
}

# ── Error trap ─────────────────────────────────────────────────────────────────
trap 'err "Pipeline failed at line ${LINENO}. Check logs in ${OUTDIR}/logs/"' ERR

# ── Main ───────────────────────────────────────────────────────────────────────
main() {
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}  Bacterial Genome Assembly Pipeline  ${NC}"
    echo -e "${GREEN}             v${VERSION}               ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"

    validate_inputs
    run_qc

    if ${QC_ONLY}; then
        trap - ERR
        ok "QC-only mode. Exiting."
        exit 0
    fi

    run_assembly
    run_polishing
    run_quality_assessment
    generate_report
}

main "$@"
