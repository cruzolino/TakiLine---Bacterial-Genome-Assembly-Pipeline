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
- Parallel QC: FastQC and NanoPlot run alongside trimming/filtering
- Automatic `pigz` detection for faster compression
- Flye version auto-detection: `--nano-hq` (Flye ≥ 2.9) or `--nano-raw` (≤ 2.8)
- Checkpoint-based resume (`-r`) — skip completed stages after interruption
- Species identification: `-S "Genus species"` (recommended), or Kraken2 auto-detects — see [Species identification](#species-identification)
- Optional typing/resistance analyses, off by default: PlasmidFinder (`--plasmid`), AMRFinderPlus (`--amr`), abricate (`--abricate`), mlst (`--mlst`)
- Genome completeness via CheckM2/CheckM and BUSCO
- Summary report at completion, as Markdown and standalone HTML
- All missing tools reported at startup in one error

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
| [Kraken2](https://github.com/DerrickWood/kraken2) | Auto-detects species when `-S` isn't passed; needs `-K <db_dir>` (or `$KRAKEN2_DB_PATH`) | `conda install -c bioconda kraken2` |

> Pass `-S "Escherichia coli"` if known — more accurate, and skips Kraken2 entirely.

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
| [Medaka](https://github.com/nanoporetech/medaka) | Nanopore consensus polishing (not used for HiFi) | `conda install -c bioconda medaka` |

### Assemblers (install only the one you need)

| Tool | Modes | Install |
|---|---|---|
| [SPAdes](https://github.com/ablab/spades) | Illumina, Hybrid | `conda install -c bioconda spades` |
| [SKESA](https://github.com/ncbi/SKESA) | Illumina | `conda install -c bioconda skesa` |
| [Unicycler](https://github.com/rrwick/Unicycler) | Illumina, Hybrid | `conda install -c bioconda unicycler` |
| [Flye](https://github.com/mikolmogorov/Flye) | Nanopore, HiFi | `conda install -c bioconda flye` |
| [Canu](https://github.com/marbl/canu) | Nanopore, HiFi | `conda install -c bioconda canu` |
| [Raven](https://github.com/lbcb-sci/raven) | Nanopore, HiFi | `conda install -c bioconda raven-assembler` |
| [Trycycler](https://github.com/rrwick/Trycycler) | Nanopore, HiFi | `conda install -c bioconda trycycler flye raven-assembler miniasm minipolish minimap2` — assembles 12 subsamples using the latter four in rotation, so all are required even though only `-a trycycler` is passed |

> **Flye:** ≥ 2.9 recommended. The pipeline auto-detects the version and picks `--nano-hq` or `--nano-raw`, but `--nano-hq` gives notably better results on R10 chemistry.

### Optional (skipped gracefully if absent)

| Tool | Purpose | Install |
|---|---|---|
| [MultiQC](https://multiqc.info/) | Aggregated QC report | `conda install -c bioconda multiqc` |
| [NanoPlot](https://github.com/wdecoster/NanoPlot) | Long-read QC plots | `conda install -c bioconda nanoplot` |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | Fast contig length filtering | `conda install -c bioconda seqkit` |
| [assembly-stats](https://github.com/sanger-pathogens/assembly-stats) | N50/total length stats | `conda install -c bioconda assembly-stats` |
| [CheckM2](https://github.com/chklovski/CheckM2) | Genome completeness (preferred) | Own env — see note below |
| [CheckM](https://ecogenomics.github.io/CheckM/) | Genome completeness (legacy) | `conda install -c bioconda checkm` |
| [BUSCO](https://busco.ezlab.org/) | Completeness (lineage, auto via `--auto-lineage-prok`; override with `-L`) | Own env — see note below |
| [pigz](https://zlib.net/pigz/) | Parallel gzip | `conda install pigz` |

### Optional typing/resistance analyses (off by default)

Passing the flag makes that tool a hard requirement at preflight, like a mandatory tool. Nothing here runs unless asked.

| Tool | Flag | Purpose | Install |
|---|---|---|---|
| [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/) | `--plasmid` | Plasmid replicon typing | `conda install -c bioconda plasmidfinder` |
| [AMRFinderPlus](https://github.com/ncbi/amr) | `--amr` | AMR genes + point mutations (`--plus` panel) | `conda install -c bioconda ncbi-amrfinderplus` |
| [abricate](https://github.com/tseemann/abricate) | `--abricate` | Secondary AMR screen, db `card` | `conda install -c bioconda abricate` |
| [mlst](https://github.com/tseemann/mlst) | `--mlst` | Sequence typing against PubMLST schemes | `conda install -c bioconda mlst` |

> **Why `card`, not abricate's default `ncbi`:** `ncbi` mirrors AMRFinderPlus's data, so pairing it with `--amr` would just duplicate results. CARD is independently curated. To use a different db, edit `ABRICATE_DB` near the top of `takiline.sh`.

> **kraken2 and skesa** install cleanly into `environment.yml`'s env if needed: `conda install -n takiline -c bioconda kraken2=2.17.1` (one package at a time).
>
> **CheckM2 and BUSCO need their own env each** — both conflict with `samtools=1.24`'s dependencies:
> ```bash
> conda create -n takiline-checkm2 -c bioconda checkm2=1.1.0
> conda run -n takiline-checkm2 checkm2 database --download   # ~1.7GB
> conda create -n takiline-busco -c bioconda busco=6.1.0
> ```
> They also can't both run automatically in the same `takiline.sh` invocation — both scripts use a `#!/usr/bin/env python3` shebang, so PATH-merging both envs makes whichever is listed first steal the other's interpreter. PATH-merge `takiline` with **one** of the two per run; run the other by hand afterward:
> ```bash
> conda activate takiline-busco
> busco -i <assembly.fasta> -o <name>_busco -l bacteria_odb10 -m genome \
>   -c <threads> --out_path <run_dir>/reports/busco/
> ```
>
> **abricate/mlst live in `environment-typing.yml`**, not `environment.yml` — both pull in perl-bioperl, which forces `samtools` down to a version too old for Pilon's `-@`/`-o` flags.

---

## Installation

```bash
git clone https://github.com/cruzolino/TakiLine---Bacterial-Genome-Assembly-Pipeline
chmod +x takiline.sh

conda env create -f environment.yml
conda activate takiline
```

`environment.yml` covers everything except `-a trycycler` mode, Nanopore polishing, and `--abricate`/`--mlst` — those need `environment-medaka.yml`, `environment-trycycler.yml`, and `environment-typing.yml` respectively. kraken2/skesa/checkm2/BUSCO aren't in any file — see the notes above.

> **Pinned, not bare package names:** an unpinned `conda create` resolves to whatever's newest at that moment, so the same command can give different versions months apart. A pinned file is a reviewable record of exactly what produced a given result.

> ### `-a trycycler` needs a second environment
> `-a trycycler` requires `trycycler`, `flye`, `minimap2`, `miniasm`, `minipolish`, and `raven` all on `$PATH` at once — `environment-trycycler.yml` covers this.

Prefer managing versions yourself? Per-tool `conda install -c bioconda <name>` commands are in [Requirements](#requirements) above.

---

## Quick Start

> Every run needs a species — `-S "Genus species"` if known (recommended), or a Kraken2 database via `-K`/`$KRAKEN2_DB_PATH` for auto-detection. See [Species identification](#species-identification).

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
| Nanopore (R10, standard) | `flye` | Robust for R9/R10; auto-selects `--nano-hq`/`--nano-raw` |
| Nanopore (legacy R9.4) | `raven` | Lightweight alternative to Flye |
| Nanopore (deep coverage) | `canu` | Higher accuracy, longer runtime |
| Nanopore (≥100x, reliability-critical) | `trycycler` | Cross-validates 12 subsample assemblies; ~6-7x slower than `flye` for marginal accuracy gain on well-behaved data |
| PacBio HiFi | `flye --hifi` | `--pacbio-hifi` mode; near-perfect assemblies |

> **Hybrid mode (`-1`/`-2` + `-l`):** only `-a unicycler` does true integrated hybrid assembly. `-a flye`/`canu`/`raven` assemble long reads only, then Pilon-polish with Illumina. `-a spades`/`skesa` are rejected in hybrid mode.

### Genome size format

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
    ├── -S "Genus species" given → recorded as-is, Kraken2 skipped
    └── -S not given → Kraken2 on reads (required; error if missing/no DB)
    │                  top species-level hit becomes the run's species;
    │                  < 70% classified reads on that hit → contamination warning
    │
    ▼
[3/5] Assembly
    ├── Illumina/Hybrid: SPAdes | SKESA | Unicycler
    ├── Nanopore/HiFi:   Flye (auto-versioned) | Canu | Raven
    │                    | Trycycler (12-subsample consensus → cluster →
    │                      auto-drop failed clusters → MSA → consensus)
    └── Contig filtering (default: ≥ 500 bp)
    │
    ▼
[4/5] Polishing
    ├── Illumina/Hybrid: Bowtie2 mapping → Pilon SNP+indel correction
    ├── Nanopore:        Medaka consensus polishing (--bacteria model selection)
    └── PacBio HiFi:     skipped (assemblers' internal polishing suffices)
    │
    ▼
[5/5] Quality Assessment
    ├── QUAST          — assembly statistics (N50, contig count, misassemblies)
    ├── CheckM2/CheckM — genome completeness and contamination
    ├── BUSCO          — lineage-specific gene completeness
    └── Optional (one flag each): --plasmid --amr --abricate --mlst
    │
    ▼
Final Report (SUMMARY.md + SUMMARY.html)
```

### Species identification

TakiLine confirms the isolate's species before assembling, to catch a mislabeled sample rather than silently assemble through it.

- **Recommended:** `-S "Genus species"` — more accurate, skips Kraken2.
- **Otherwise Kraken2 runs automatically**, and the species-level hit with the highest read share becomes the species of record (`SUMMARY.md`, `logs/.species`).
- **Kraken2 is conditionally required**: no `-S` and no working Kraken2 (missing binary or DB) is a hard error at preflight.
- **Low-confidence calls warn, not abort**: under 70% on the top hit logs a contamination warning but the run continues.
- **On resume**, species is restored from `logs/.species`; a different `-S` on resume is caught by the fingerprint check.

### Key design decisions

**Pinned conda over Docker** — TakiLine is single-run and local; conda gives adequate isolation for that. Containerization would mostly add overhead without solving the real risk (version drift over time), which the pinned `environment.yml` files already handle.

**Parallel QC** — FastQC/NanoPlot run in the background while `fastp`/`filtlong` run in the foreground, saving ~30–120 s with no added resource contention.

**Flye version auto-detection** — queries `flye --version` and selects `--nano-hq` (≥ 2.9) or `--nano-raw` (≤ 2.8), warning if outdated.

**Platform-aware long-read QC** — Filtlong uses `--min_mean_q 7` for Nanopore and `--min_mean_q 20` for HiFi, matching each platform's quality distribution.

**Medaka polishing for Nanopore** — assemblers' internal consensus still leaves ONT-specific indel error; Medaka measurably improves gene completeness. `--bacteria` auto-model needs basecaller metadata in read headers (often absent from SRA/ENA data) — use `-M` if auto-detection fails. Skipped for HiFi.

**BUSCO auto-lineage** — defaults to `--auto-lineage-prok` for a tighter completeness score than generic `bacteria_odb10`. `-L` forces an exact lineage; a failed auto-placement retries once with `bacteria_odb10`.

**Trycycler auto-pruning** — runs non-interactively: a cluster that fails reconcile has the offending contig dropped and reconcile retried automatically, instead of waiting for manual review.

**Pilon memory guard** — checks available RAM before running Pilon and warns if `-m` exceeds free memory, rather than crashing mid-run.

**Typing/resistance analyses are opt-in** — PlasmidFinder, AMRFinderPlus, abricate, mlst answer "what is this isolate," not assembly quality (QUAST/CheckM2/BUSCO's job, on by default). Once a flag is passed, that tool becomes a mandatory-tool preflight check.

**Annotation is out of scope** — keeping assembly and annotation separate lets each step be run, updated, or repeated independently. See [Suggested Next Steps](#suggested-next-steps).

---

## Output Structure

By default, once the full pipeline (not `-q`) completes, heavy intermediates are cleaned up automatically — `00_raw_data/`, trimmed/filtered FASTQs, and the `02_assembly/`/`03_polishing/` working directories are removed, leaving reports, logs, and the final assembly:

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

Pass **`-k`** to keep intermediates instead (raw-data symlinks, trimmed/filtered FASTQs, `02_assembly/{assembler}/`, `03_polishing/` BAMs) — useful for debugging. With `-k`, the final assembly stays at its per-assembler path (e.g. `02_assembly/{assembler}/contigs_filtered.fasta`, or `03_polishing/{sample}_pilon.fasta` for Illumina/hybrid) instead of being moved to `{input}.fasta`.

The final assembly FASTA path is printed at completion and recorded in `reports/SUMMARY.md`/`SUMMARY.html`.

---

## Resume Mode

If interrupted (node failure, timeout, manual kill), re-run the exact same command with `-r` appended:

```bash
# Original command
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

# Resume after interruption
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
```

Each stage writes a sentinel file (e.g. `logs/.done_assembly`) on success. With `-r`, any stage whose sentinel exists is skipped. To force a stage to re-run, delete its sentinel then resume:

```bash
# Force re-assembly only, keep QC results
rm takiline/logs/.done_assembly
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
```

---

## Suggested Next Steps

| Task | Tool |
|---|---|
| Gene annotation | [Bakta](https://github.com/oschwengers/bakta) (INSDC-ready GFF3/GBFF) |
| AMR gene detection | Built into TakiLine via `--amr`/`--abricate` |
| MLST typing | Built into TakiLine via `--mlst` ([mlst](https://github.com/tseemann/mlst), Torsten Seemann) |
| Plasmid replicon typing | Built into TakiLine via `--plasmid` ([PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/)) |
| Phylogenetics | [IQ-TREE2](http://www.iqtree.org/) or [FastTree](http://www.microbesonline.org/fasttree/) |
| Pan-genome | [Panaroo](https://github.com/gtonkinhill/panaroo) (preferred) or [Roary](https://sanger-pathogens.github.io/Roary/) |
| Secondary metabolites | [antiSMASH](https://antismash.secondarymetabolites.org/) |
