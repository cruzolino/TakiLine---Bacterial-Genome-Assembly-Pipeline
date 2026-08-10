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
- Six assembler choices: SPAdes, SKESA, Unicycler, Flye, Canu, Raven
- Parallel QC: FastQC and NanoPlot run concurrently with trimming/filtering steps
- Automatic `pigz` detection for faster compression
- **Flye version auto-detection**: selects `--nano-hq` (Flye ≥ 2.9) or `--nano-raw` (Flye ≤ 2.8) at runtime
- Checkpoint-based **resume** (`-r`) — skip completed stages after interruption
- Plasmid detection via PlasmidFinder (optional, on by default)
- Genome completeness assessment via CheckM2/CheckM and BUSCO
- Markdown summary report generated at completion
- All missing tools reported at startup in a single error — no fix-one-run-find-next cycle

---

## Requirements

### Mandatory (all modes)

| Tool | Purpose | Install |
|---|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Raw/trimmed read QC | `conda install -c bioconda fastqc` |
| [fastp](https://github.com/OpenGFP/fastp) | Illumina adapter trimming | `conda install -c bioconda fastp` |
| [QUAST](https://quast.sourceforge.net/) | Assembly quality metrics | `conda install -c bioconda quast` |

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

### Assemblers (install only the one you need)

| Tool | Modes | Install |
|---|---|---|
| [SPAdes](https://github.com/ablab/spades) | Illumina, Hybrid | `conda install -c bioconda spades` |
| [SKESA](https://github.com/ncbi/SKESA) | Illumina | `conda install -c bioconda skesa` |
| [Unicycler](https://github.com/rrwick/Unicycler) | Illumina, Hybrid | `conda install -c bioconda unicycler` |
| [Flye](https://github.com/mikolmogorov/Flye) | Nanopore, HiFi | `conda install -c bioconda flye` |
| [Canu](https://github.com/marbl/canu) | Nanopore, HiFi | `conda install -c bioconda canu` |
| [Raven](https://github.com/lbcb-sci/raven) | Nanopore, HiFi | `conda install -c bioconda raven-assembler` |

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
| [BUSCO](https://busco.ezlab.org/) | Genome completeness (lineage) | `conda install -c bioconda busco` |
| [PlasmidFinder](https://cge.food.dtu.dk/services/PlasmidFinder/) | Plasmid replicon detection | `conda install -c bioconda plasmidfinder` |
| [pigz](https://zlib.net/pigz/) | Parallel gzip (faster compression) | `conda install pigz` |

---

## Installation

```bash
# Clone or download the script
git clone https://github.com/cruzolino/TakiLine---Bacterial-Genome-Assembly-Pipeline
chmod +x takiline.sh

# Recommended: create a dedicated conda environment
conda create -n takiline \
    fastqc fastp quast bowtie2 samtools pilon \
    filtlong spades unicycler "flye>=2.9" \
    multiqc nanoplot seqkit assembly-stats \
    checkm2 busco plasmidfinder pigz
conda activate takiline
```

---

## Quick Start

```bash
# Illumina paired-end
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

# Nanopore only
./takiline.sh -l ont.fq.gz -s Salmonella -a flye -g 4.8m -t 16

# Hybrid (Illumina + Nanopore)
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz -s Klebsiella -a unicycler

# PacBio HiFi
./takiline.sh -l hifi.fq.gz -s Pseudomonas -a flye --hifi -g 6.5m

# Resume an interrupted run
./takiline.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
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
                       flye|canu|raven (long-read)        [default: spades]

General:
  -o DIR    Output directory          [default: takiline]
  -t INT    Threads                   [default: auto-detect]
  -m STR    Memory limit              [default: 32G]
  -s STR    Sample name               [default: isolate]
  -g STR    Expected genome size      [default: 5m]
  -c INT    Min contig length (bp)    [default: 500]
  -P STR    PlasmidFinder DB          [default: enterobacteriaceae]
  -D DIR    PlasmidFinder DB path     [default: auto-detect]
  -x        Skip plasmid detection
  -q        QC only (skip assembly)
  -r        Resume from last checkpoint
  -k        Keep intermediates (skip automatic cleanup)
  -h        Show help
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
| PacBio HiFi | `flye --hifi` | `--pacbio-hifi` mode; produces near-perfect assemblies |

> **Hybrid mode (`-1`/`-2` + `-l`) note:** only `-a unicycler` performs true integrated hybrid assembly (long reads close gaps in the assembly graph). `-a flye`/`canu`/`raven` in hybrid mode assemble the long reads only and use the Illumina reads for post-assembly Pilon polishing — a valid but different strategy. `-a spades`/`skesa` cannot use long reads at all and are rejected in hybrid mode to avoid silently discarding the `-l` input.

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
[1/4] Quality Control
    ├── Illumina: FastQC (raw, parallel) → fastp trimming → FastQC (trimmed) → MultiQC
    └── Long-read: NanoPlot (parallel) → Filtlong (20× coverage target)
    │         └── Quality threshold: Q7 (Nanopore) | Q20 (PacBio HiFi)
    │
    ▼
[2/4] Assembly
    ├── Illumina/Hybrid: SPAdes | SKESA | Unicycler
    ├── Nanopore/HiFi:   Flye (auto-versioned) | Canu | Raven
    └── Contig filtering (default: ≥ 500 bp)
    │
    ▼
[3/4] Polishing
    ├── Illumina/Hybrid: Bowtie2 mapping → Pilon SNP+indel correction
    └── Long-read only:  skipped (Flye/Canu/Raven perform internal polishing)
    │                    assembly path persisted for downstream stages
    │
    ▼
[4/4] Quality Assessment
    ├── QUAST          — assembly statistics (N50, contig count, misassemblies)
    ├── CheckM2/CheckM — genome completeness and contamination
    ├── BUSCO          — lineage-specific gene completeness (bacteria_odb10)
    └── PlasmidFinder  — plasmid replicon typing (disable with -x)
    │
    ▼
Final Report (SUMMARY.md)
```

### Key design decisions

**Parallel QC** — FastQC on raw reads and NanoPlot are launched in the background while `fastp`/`filtlong` run in the foreground. This saves 30–120 s on typical datasets with no additional resource contention.

**Flye version auto-detection** — the pipeline queries `flye --version` at runtime and selects `--nano-hq` (Flye ≥ 2.9) or `--nano-raw` (Flye ≤ 2.8) automatically, with a warning prompting the user to upgrade if an older version is detected.

**Platform-aware long-read QC** — Filtlong applies `--min_mean_q 7` for Nanopore reads and `--min_mean_q 20` for PacBio HiFi reads, reflecting the distinct quality score distributions of each platform.

**No external long-read polishing** — Flye, Canu, and Raven perform iterative internal polishing. A separate Medaka pass is redundant and a common source of pipeline crashes.

**Pilon memory guard** — the pipeline checks available RAM before running Pilon and warns if the `-m` limit exceeds free memory, rather than silently crashing mid-run.

**Annotation is out of scope** — keeping assembly and annotation separate allows each step to be run, updated, or repeated independently. See [Suggested Next Steps](#suggested-next-steps) for recommended tools.

---

## Output Structure

By default, once the full pipeline (not `-q`) completes, heavy intermediates are cleaned up automatically — `00_raw_data/`, the trimmed/filtered FASTQs, and the whole `02_assembly/`/`03_polishing/` working directories are removed, leaving only the reports, logs, and the final assembly:

```
takiline/
├── final_assembly.fasta              # The final genome — moved here by cleanup
├── 01_qc/
│   ├── {sample}_fastp.html           # QC reports only; trimmed/filtered
│   ├── nanoplot/                     # FASTQs are removed after assembly
│   └── post_trim/
├── 05_reports/
│   ├── SUMMARY.md                    # Human-readable pipeline report
│   ├── quast/
│   ├── checkm/ or checkm2/
│   ├── busco/
│   └── plasmidfinder/
└── logs/
    ├── .assembly_path                # Points at final_assembly.fasta
    ├── .done_qa                      # qc/assembly/polishing sentinels are
    ├── fastp.log                     # removed along with their outputs, so
    ├── fastqc_raw.log                # a later -r resume redoes those stages
    ├── fastqc_post.log               # from the original input files
    ├── {assembler}.log
    ├── bowtie2.log
    ├── pilon.log
    └── quast.log
```

Pass **`-k`** to keep intermediates instead (raw-data symlinks, trimmed/filtered FASTQs, `02_assembly/{assembler}/` working directories, `03_polishing/` BAMs) — useful for debugging an assembly. With `-k`, the final assembly stays at its normal per-assembler path (e.g. `02_assembly/{assembler}/contigs_filtered.fasta`, or `03_polishing/{sample}_pilon.fasta` for Illumina/hybrid) instead of being moved to `final_assembly.fasta`.

The **final assembly FASTA** path is printed at completion and recorded in `05_reports/SUMMARY.md`.

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
| AMR gene detection | [AMRFinderPlus](https://github.com/ncbi/amr) or [ResFinder](https://cge.food.dtu.dk/services/ResFinder/) |
| MLST typing | [mlst](https://github.com/tseemann/mlst) (Torsten Seemann) |
| NCBI submission | Bakta outputs + `table2asn` |
| Phylogenetics | [IQ-TREE2](http://www.iqtree.org/) or [FastTree](http://www.microbesonline.org/fasttree/) |
| Pan-genome | [Panaroo](https://github.com/gtonkinhill/panaroo) (preferred) or [Roary](https://sanger-pathogens.github.io/Roary/) |
| Secondary metabolites | [antiSMASH](https://antismash.secondarymetabolites.org/) |
