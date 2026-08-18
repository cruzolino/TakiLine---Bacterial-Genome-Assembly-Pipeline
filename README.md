# TakiLine - Bacterial Genome Assembly Pipeline

<img width="314" height="314" alt="Takiline" src="https://github.com/user-attachments/assets/0c909374-af2b-443a-ae72-6adaaeef48cc" />


A streamlined bash pipeline for *de novo* bacterial genome assembly from Illumina paired-end, Nanopore, PacBio HiFi, or hybrid reads. Produces a polished, QC-evaluated assembly FASTA ready for downstream annotation and comparative genomics.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
- [Output Structure](#output-structure)
- [Resume Mode](#resume-mode)
- [Suggested Next Steps](#suggested-next-steps)

---

## Features

- Four sequencing modes: **Illumina PE**, **Nanopore**, **PacBio HiFi**, **Hybrid**
- Seven assembler choices: SPAdes, SKESA, Unicycler, Flye, Canu, Raven, Trycycler
- Parallel QC: FastQC and NanoPlot run concurrently with trimming/filtering steps
- Automatic `pigz` detection for faster compression
- **Flye version auto-detection**: selects `--nano-hq` (Flye ≥ 2.9) or `--nano-raw` (Flye ≤ 2.8) at runtime
- Checkpoint-based **resume** (`-r`) — skip completed stages after interruption
- **Species identification**: inform it with `-S "Genus species"` (recommended), or Kraken2 auto-detects it from reads — see [Species identification](#species-identification)
- Optional typing/resistance analyses, off by default: PlasmidFinder (`--plasmid`), AMRFinderPlus (`--amr`), abricate (`--abricate`), mlst (`--mlst`) — see [Optional typing/resistance analyses](#optional-typingresistance-analyses-off-by-default)
- Genome completeness assessment via CheckM2/CheckM and BUSCO
- Summary report generated at completion, as both Markdown and standalone HTML
- All missing tools reported at startup in a single error — no fix-one-run-find-next cycle

---

## Requirements

### Mandatory (all modes)

| Tool | Purpose | Install |
|---|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Raw/trimmed read QC | `conda install -c bioconda fastqc` |
| [fastp](https://github.com/OpenGFP/fastp) | Illumina adapter trimming | `conda install -c bioconda fastp` |
| [QUAST](https://quast.sourceforge.net/) | Assembly quality metrics | `conda install -c bioconda quast` |

### Species identification (required unless `-S` is given)

| Tool | Purpose | Install |
|---|---|---|
| [Kraken2](https://github.com/DerrickWood/kraken2) | Auto-detects species from reads when `-S` isn't passed, and needs `-K <db_dir>` (or `$KRAKEN2_DB_PATH`) pointing at a Kraken2 database | `conda install -c bioconda kraken2` |

> Pass `-S "Escherichia coli"` if you already know the species — more accurate than auto-detection, and skips Kraken2 entirely. See [Species identification](#species-identification).

### Illumina / Hybrid

| Tool | Purpose | Install |
|---|---|---|
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | Read mapping for polishing | `conda install -c bioconda bowtie2` |
| [SAMtools](https://www.htslib.org/) | BAM sorting and indexing | `conda install -c bioconda samtools` |
| [Pilon](https://github.com/broadinstitute/pilon) | Illumina-based polishing | `conda install -c bioconda pilon` |

### Long-read / Hybrid

| Tool | Purpose | Install |
|---|---|---|
| [Filtlong](https://github.com/rrwick/Filtlong) | Long-read quality filtering | `conda install -c bioconda filtlong` |
| [Medaka](https://github.com/nanoporetech/medaka) | Nanopore-only consensus polishing (required for Nanopore mode; not used for PacBio HiFi) | `conda install -c bioconda medaka` |

### Assemblers (install only the one you need)

| Tool | Modes | Install |
|---|---|---|
| [SPAdes](https://github.com/ablab/spades) | Illumina, Hybrid | `conda install -c bioconda spades` |
| [SKESA](https://github.com/ncbi/SKESA) | Illumina | `conda install -c bioconda skesa` |
| [Unicycler](https://github.com/rrwick/Unicycler) | Illumina, Hybrid | `conda install -c bioconda unicycler` |
| [Flye](https://github.com/mikolmogorov/Flye) | Nanopore, HiFi | `conda install -c bioconda flye` |
| [Canu](https://github.com/marbl/canu) | Nanopore, HiFi | `conda install -c bioconda canu` |
| [Raven](https://github.com/lbcb-sci/raven) | Nanopore, HiFi | `conda install -c bioconda raven-assembler` |
| [Trycycler](https://github.com/rrwick/Trycycler) | Nanopore, HiFi | `conda install -c bioconda trycycler flye raven-assembler miniasm minipolish minimap2` — assembles each of its 12 subsamples with the latter four in rotation, so all are required even though only `-a trycycler` is passed |

> **Flye version note:** Flye ≥ 2.9 is recommended. The pipeline auto-detects the installed version and selects the appropriate Nanopore flag (`--nano-hq` for ≥ 2.9, `--nano-raw` for ≤ 2.8), but `--nano-hq` yields significantly better results with R10 chemistry reads.

### Optional (skipped gracefully if absent)

| Tool | Purpose | Install |
|---|---|---|
| [MultiQC](https://multiqc.info/) | Aggregated QC report | `conda install -c bioconda multiqc` |
| [NanoPlot](https://github.com/wdecoster/NanoPlot) | Long-read QC plots | `conda install -c bioconda nanoplot` |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | Fast contig length filtering | `conda install -c bioconda seqkit` |
| [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) | N50/total length stats | `conda install -c bioconda assembly-stats` |
| [CheckM2](https://github.com/chklovski/CheckM2) | Genome completeness (preferred) | `conda install -c bioconda checkm2` |
| [CheckM](https://ecogenomics.github.io/CheckM/) | Genome completeness (legacy) | `conda install -c bioconda checkm` |
| [BUSCO](https://busco.ezlab.org/) | Genome completeness (lineage, auto-selected via `--auto-lineage-prok`; override with `-L`) | `conda install -c bioconda busco` |
| [pigz](https://zlib.net/pigz/) | Parallel gzip (faster compression) | `conda install pigz` |

### Optional typing/resistance analyses (off by default)

Not best-effort like the table above: pass the flag and the tool becomes a hard requirement (fails fast at preflight, like a mandatory tool). Nothing here runs unless asked.

| Tool | Flag | Purpose | Install |
|---|---|---|---|
| [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/) | `--plasmid` | Plasmid replicon typing | `conda install -c bioconda plasmidfinder` |
| [AMRFinderPlus](https://github.com/ncbi/amr) | `--amr` | AMR genes + point mutations (`--plus` stress/virulence panel) | `conda install -c bioconda ncbi-amrfinderplus` |
| [abricate](https://github.com/tseemann/abricate) | `--abricate` | Secondary AMR screen, db `card` | `conda install -c bioconda abricate` |
| [mlst](https://github.com/tseemann/mlst) | `--mlst` | Sequence typing against PubMLST schemes | `conda install -c bioconda mlst` |

> **Why `card`, not abricate's own default `ncbi`:** `ncbi` mirrors AMRFinderPlus's own data — pairing it with `--amr` would just run the same database twice. CARD is independently curated, so the two flags give genuinely different results worth comparing. To use a different db instead, edit `ABRICATE_DB` near the top of `takiline.sh` (no CLI flag for this yet).

---

## Installation

```bash
# Clone or download the script
git clone https://github.com/cruzolino/TakiLine---Bacterial-Genome-Assembly-Pipeline
chmod +x takiline.sh

# Create the pinned environment (exact tested versions, not "whatever's newest today")
conda env create -f environment.yml
conda activate takiline
```

`environment.yml` covers everything except `-a trycycler` mode and Nanopore polishing — those need `environment-medaka.yml` and `environment-trycycler.yml` (see below). Versions are pinned to what's already proven to run together on this machine, not just "current bioconda release" — each file's header documents its own provenance.

> **Why pinned, not bare package names?** An unpinned `conda create` resolves to whatever's newest *the moment you run it* — the same command can silently give different tool versions months apart. A pinned file is a reviewable record of exactly which versions produced a given result; it does not give Docker-level OS isolation (see [Key design decisions](#key-design-decisions) for why that's skipped here).

> ### ⚑ `-a trycycler` needs a second environment, and it isn't proven to work yet
> TakiLine's `-a trycycler` path needs `trycycler`, `flye`, `minimap2`, `miniasm`, `minipolish`, and `raven` all on `$PATH` in one run. On this machine, `trycycler`/`minipolish` only ever existed in their own env (Python 3.13, numpy 2.5.x) — never alongside the rest (Python 3.12, numpy 1.26.x). That numpy 1.x-vs-2.x split is a classic bioconda solver-conflict signature. `environment-trycycler.yml` bundles both sides as a best-effort pin, but the combined solve hasn't actually been tested — treat it as a starting point, not a known-good environment.

Prefer managing versions yourself? The per-tool `conda install -c bioconda <name>` commands are still listed in [Requirements](#requirements) above.

---

## Quick Start

> Every run below needs a species — either `-S "Genus species"` if you already know it (recommended), or a Kraken2 database via `-K`/`$KRAKEN2_DB_PATH` so the pipeline can auto-detect it. See [Species identification](#species-identification).

```bash
# Illumina paired-end, species known
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -S "Escherichia coli" -t 16

# Illumina paired-end, species unknown — Kraken2 auto-detects it
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -K /path/to/kraken2_db

# Nanopore only
./takiline.sh -l ont.fq.gz -s Salmonella -S "Salmonella enterica" -a flye -g 4.8m -t 16

# Nanopore, high-confidence (needs ≥100x reads)
./takiline.sh -l ont.fq.gz -s Salmonella -S "Salmonella enterica" -a trycycler -g 4.8m -t 16

# Hybrid (Illumina + Nanopore)
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz -s Klebsiella -S "Klebsiella pneumoniae" -a unicycler

# PacBio HiFi
./takiline.sh -l hifi.fq.gz -s Pseudomonas -S "Pseudomonas aeruginosa" -a flye --hifi -g 6.5m

# Resume an interrupted run
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -S "Escherichia coli" -t 16 -r
```

---

## Usage

```
Usage:  takiline.sh [OPTIONS]

Input (at least one required):
  -1 FILE   Illumina forward reads (R1.fastq.gz)
  -2 FILE   Illumina reverse reads (R2.fastq.gz)
  -l FILE   Long reads – Nanopore or PacBio HiFi (.fastq[.gz])
  --hifi    Treat -l reads as PacBio HiFi (activates flye --pacbio-hifi)

Assembly:
  -a STR    Assembler: spades|skesa|unicycler (Illumina/hybrid)
                       flye|canu|raven|trycycler (long-read) [default: spades]

General:
  -o DIR    Output directory          [default: takiline]
  -t INT    Threads                   [default: auto-detect]
  -m STR    Memory limit              [default: 32G]
  -s STR    Sample name               [default: isolate]
  -g STR    Expected genome size      [default: 5m]
  -c INT    Min contig length (bp)    [default: 500]
  -S STR    Expected species, e.g. "Escherichia coli" [default: auto-detect via Kraken2]
  -K DIR    Kraken2 database path     [default: $KRAKEN2_DB_PATH env var]
  -P STR    PlasmidFinder DB          [default: enterobacteriaceae]
  -D DIR    PlasmidFinder DB path     [default: auto-detect]
  -M STR    Medaka model (Nanopore)   [default: auto-select via --bacteria]
  -L STR    BUSCO lineage odb10       [default: auto-select via --auto-lineage-prok]
  -q        QC only (skip assembly)
  -r        Resume from last checkpoint
  -k        Keep intermediates (skip automatic cleanup)
  -h        Show help

Optional typing/resistance analyses (off by default):
  --plasmid   PlasmidFinder replicon typing (uses -P/-D above)
  --amr       AMRFinderPlus (AMR genes + point mutations, --plus panel)
  --abricate  abricate secondary AMR screen (db: card)
  --mlst      mlst sequence typing (PubMLST schemes)
```

### Assembler selection guide

| Input | Recommended assembler | Notes |
|---|---|---|
| Illumina PE only | `spades` | `--isolate` mode; best for single isolates |
| Illumina PE only | `skesa` | Faster than SPAdes; fewer mis-assemblies on some datasets |
| Illumina + Nanopore | `unicycler` | Closes chromosomal and plasmid gaps using long reads |
| Nanopore (R10, standard) | `flye` | Robust for R9/R10 chemistry; auto-selects `--nano-hq` or `--nano-raw` |
| Nanopore (legacy R9.4) | `raven` | Lightweight alternative to Flye |
| Nanopore (deep coverage) | `canu` | Higher accuracy at the cost of longer runtime |
| Nanopore (≥100x, reliability-critical) | `trycycler` | Cross-validates 12 subsample assemblies to catch single-assembler misassembly. ~6-7x slower than `flye` for only marginal accuracy gain on well-behaved data — use it for the safety net, not for speed |
| PacBio HiFi | `flye --hifi` | `--pacbio-hifi` mode; produces near-perfect assemblies |

> **Hybrid mode (`-1`/`-2` + `-l`):** only `-a unicycler` does true integrated hybrid assembly. `-a flye`/`canu`/`raven` assemble the long reads only, then Pilon-polish with Illumina — valid, but a different strategy. `-a spades`/`skesa` are rejected in hybrid mode rather than silently dropping `-l`.

### Genome size format

Pass a number followed by `m` (megabases) or `g` (gigabases):

```
-g 5m      # 5 Mb  (typical E. coli)
-g 4.8m    # 4.8 Mb
-g 0.5g    # 500 Mb (unusual for bacteria; included for completeness)
```

---

## Pipeline Overview

```
Input reads
    │
    ▼
[1/5] Quality Control
    ├── Illumina: FastQC (raw, parallel) → fastp trimming → FastQC (trimmed) → MultiQC
    └── Long-read: NanoPlot (parallel) → Filtlong (20× coverage target)
    │         └── Quality threshold: Q7 (Nanopore) | Q20 (PacBio HiFi)
    │
    ▼
[2/5] Species Identification
    ├── -S "Genus species" given → recorded as-is, Kraken2 skipped (recommended)
    └── -S not given → Kraken2 on reads (required; error if missing/no DB)
    │                  top species-level hit becomes the run's species;
    │                  < 70% of classified reads on that hit → contamination warning
    │
    ▼
[3/5] Assembly
    ├── Illumina/Hybrid: SPAdes | SKESA | Unicycler
    ├── Nanopore/HiFi:   Flye (auto-versioned) | Canu | Raven
    │                    | Trycycler (12-subsample consensus: Flye/Miniasm+
    │                      Minipolish/Raven rotated → cluster → auto-drop
    │                      clusters that fail reconcile → MSA → consensus)
    └── Contig filtering (default: ≥ 500 bp)
    │
    ▼
[4/5] Polishing
    ├── Illumina/Hybrid: Bowtie2 mapping → Pilon SNP+indel correction
    ├── Nanopore:        Medaka consensus polishing (--bacteria model selection)
    └── PacBio HiFi:     skipped (Flye/Canu/Raven/Trycycler perform internal polishing)
    │                    assembly path persisted for downstream stages
    │
    ▼
[5/5] Quality Assessment
    ├── QUAST          — assembly statistics (N50, contig count, misassemblies)
    ├── CheckM2/CheckM — genome completeness and contamination
    ├── BUSCO          — lineage-specific gene completeness (auto-selected via
    │                    --auto-lineage-prok, or forced with -L)
    └── Optional (off by default, one flag each):
                         --plasmid  → PlasmidFinder (replicon typing)
                         --amr      → AMRFinderPlus (AMR genes + point mutations)
                         --abricate → abricate, db=card (2nd AMR opinion)
                         --mlst     → mlst (sequence typing)
    │
    ▼
Final Report (SUMMARY.md + SUMMARY.html)
```

### Species identification

TakiLine wants to know the isolate's species before assembling it, to catch an unidentified or mislabeled sample rather than silently assemble through it.

- **Recommended: inform it** with `-S "Genus species"` — more accurate than auto-detection, and skips Kraken2 entirely.
- **Otherwise Kraken2 runs automatically** on the read set and becomes the species of record (recorded in `SUMMARY.md` and `logs/.species`), taking the species-level hit with the highest percentage of classified reads.
- **Kraken2 is conditionally required, not skippable**: no `-S` and no working Kraken2 (missing binary, or `-K`/`$KRAKEN2_DB_PATH` not resolving to a database) is a hard error at preflight, same as a missing mandatory tool.
- **Low-confidence calls warn, not abort**: under 70% of classified reads on the top hit logs a possible-contamination warning and points at the Kraken2 report, but the run continues.
- **On resume**, the species is restored from `logs/.species` without re-running Kraken2; a different `-S` on resume is caught by the resume-fingerprint check.

### Key design decisions

**Pinned conda environments over Docker** — TakiLine is deliberately single-run and local-only; containerization mainly pays off across *different* environments (laptop, cluster, someone else's machine), none of which apply here, and conda already gives adequate isolation for one local user. nf-core-style pipelines get per-process containers nearly free because Nextflow dispatches each tool to its own container automatically — a bash script has no such dispatch, so bolting on one giant image would just inherit TakiLine's real cross-tool dependency conflicts (see the Trycycler callout in [Installation](#installation)) directly. Pinned `environment.yml` files are the right-sized fix for the actual risk (drift over time), revisited if multi-machine use, publishing, or a Nextflow move ever happens.

**Parallel QC** — FastQC on raw reads and NanoPlot are launched in the background while `fastp`/`filtlong` run in the foreground. This saves 30–120 s on typical datasets with no additional resource contention.

**Flye version auto-detection** — the pipeline queries `flye --version` at runtime and selects `--nano-hq` (Flye ≥ 2.9) or `--nano-raw` (Flye ≤ 2.8) automatically, with a warning prompting the user to upgrade if an older version is detected.

**Platform-aware long-read QC** — Filtlong applies `--min_mean_q 7` for Nanopore reads and `--min_mean_q 20` for PacBio HiFi reads, reflecting the distinct quality score distributions of each platform.

**Medaka polishing for Nanopore** — Flye/Canu/Raven/Trycycler's internal consensus still leaves ONT-specific indel error; a Medaka pass measurably improves gene completeness. `--bacteria` auto-model needs basecaller metadata in the read headers (absent from most SRA/ENA data) — use `-M` when auto-detection fails. Skipped for HiFi, whose consensus is already near-perfect.

**BUSCO auto-lineage** — defaults to `--auto-lineage-prok` instead of generic `bacteria_odb10`, placing the assembly against prokaryote marker trees for a tighter completeness score (e.g. `enterobacterales_odb10` for *E. coli*). `-L <lineage_odb10>` forces an exact lineage; a failed auto-placement retries once with `bacteria_odb10` and warns rather than aborting.

**Trycycler auto-pruning** — the pipeline runs non-interactively, so a cluster that fails reconcile has the offending contig dropped and reconcile retried automatically, instead of waiting for the manual review Trycycler's docs expect.

**Pilon memory guard** — the pipeline checks available RAM before running Pilon and warns if the `-m` limit exceeds free memory, rather than silently crashing mid-run.

**Typing/resistance analyses are opt-in, not opt-out** — PlasmidFinder, AMRFinderPlus, abricate, and mlst answer "what is this isolate," not "is this assembly any good" (QUAST/CheckM2/BUSCO's job, which stay on by default). All four require an explicit flag, and once passed, that tool joins the mandatory-tool preflight check — missing it is a hard error at startup, not a warning after a multi-hour assembly. CheckM2/BUSCO stay best-effort regardless, since they aren't flag-gated at all. See [Optional typing/resistance analyses](#optional-typingresistance-analyses-off-by-default) for the `--abricate` database choice.

**Annotation is out of scope** — keeping assembly and annotation separate allows each step to be run, updated, or repeated independently. See [Suggested Next Steps](#suggested-next-steps) for recommended tools.

---

## Output Structure

By default, once the full pipeline (not `-q`) completes, heavy intermediates are cleaned up automatically — `00_raw_data/`, the trimmed/filtered FASTQs, and the whole `02_assembly/`/`03_polishing/` working directories are removed, leaving only the reports, logs, and the final assembly:

```
takiline/
├── {input}.fasta                     # The final genome — named after the input
│                                      # file (e.g. cft073_55x.fastq.gz → cft073_55x.fasta)
├── qc/
│   ├── {sample}_fastp.html           # QC reports only; trimmed/filtered
│   ├── nanoplot/                     # FASTQs are removed after assembly
│   └── post_trim/
├── reports/
│   ├── SUMMARY.md                    # Human-readable pipeline report
│   ├── SUMMARY.html                  # Same report, as a standalone HTML page
│   ├── quast/
│   ├── checkm/ or checkm2/
│   ├── busco/
│   ├── kraken2/                      # only populated if -S was NOT passed
│   ├── plasmidfinder/                # only populated if --plasmid was passed
│   ├── amrfinderplus/                # only populated if --amr was passed
│   ├── abricate/                     # only populated if --abricate was passed
│   └── mlst/                         # only populated if --mlst was passed
└── logs/
    ├── .assembly_path                # Points at the final assembly FASTA
    ├── .species                      # Informed (-S) or Kraken2-detected species
    ├── .done_qa                      # qc/assembly/polishing sentinels are
    ├── fastp.log                     # removed along with their outputs, so
    ├── fastqc_raw.log                # a later -r resume redoes those stages
    ├── fastqc_post.log               # from the original input files
    ├── kraken2.log
    ├── {assembler}.log
    ├── bowtie2.log
    ├── pilon.log
    └── quast.log
```

Pass **`-k`** to keep intermediates instead (raw-data symlinks, trimmed/filtered FASTQs, `02_assembly/{assembler}/` working directories, `03_polishing/` BAMs) — useful for debugging an assembly. With `-k`, the final assembly stays at its normal per-assembler path (e.g. `02_assembly/{assembler}/contigs_filtered.fasta`, or `03_polishing/{sample}_pilon.fasta` for Illumina/hybrid) instead of being moved to `{input}.fasta`.

The **final assembly FASTA** path is printed at completion and recorded in `reports/SUMMARY.md`/`SUMMARY.html`.

---

## Resume Mode

If the pipeline is interrupted (node failure, timeout, manual kill), re-run the exact same command with `-r` appended:

```bash
# Original command
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

# Resume after interruption
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
```

Each stage writes a sentinel file (e.g. `logs/.done_assembly`) on successful completion. With `-r`, any stage whose sentinel exists is skipped entirely. To force a specific stage to re-run, delete its sentinel then resume:

```bash
# Force re-assembly only, keep QC results
rm takiline/logs/.done_assembly
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
```

---

## Suggested Next Steps

After the pipeline completes, the final polished assembly FASTA is ready for:

| Task | Tool |
|---|---|
| Gene annotation | [Bakta](https://github.com/oschwengers/bakta) (INSDC-ready GFF3/GBFF) |
| AMR gene detection | Built into TakiLine via `--amr`/`--abricate` — see [Optional typing/resistance analyses](#optional-typingresistance-analyses-off-by-default) |
| MLST typing | Built into TakiLine via `--mlst` ([mlst](https://github.com/tseemann/mlst), Torsten Seemann) |
| Plasmid replicon typing | Built into TakiLine via `--plasmid` ([PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/)) |
| NCBI submission | Bakta outputs + `table2asn` |
| Phylogenetics | [IQ-TREE2](http://www.iqtree.org/) or [FastTree](http://www.microbesonline.org/fasttree/) |
| Pan-genome | [Panaroo](https://github.com/gtonkinhill/panaroo) (preferred) or [Roary](https://sanger-pathogens.github.io/Roary/) |
| Secondary metabolites | [antiSMASH](https://antismash.secondarymetabolites.org/) |
