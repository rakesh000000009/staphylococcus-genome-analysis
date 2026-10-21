# Staphylococcus Genome Analysis Pipeline

[![Bioinformatics - WGS](https://img.shields.io/badge/Bioinformatics-WGS%20Assembly%20%26%20Annotation-blue.svg)](#)
[![Target Organism](https://img.shields.io/badge/Organism-Staphylococcus%20aureus-green.svg)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An end-to-end, production-grade bioinformatics pipeline tailored for Whole-Genome Sequencing (WGS) analysis of *Staphylococcus* species (including Methicillin-resistant *Staphylococcus aureus* / MRSA). 

The workflow performs automated paired-end quality control and trimming (`fastp`), *de novo* genome assembly (`SPAdes --isolate`), quality and completeness evaluation (`QUAST` & `BUSCO`), standardized genomic feature annotation (`Prokka`), and molecular epidemiological sequence typing (`mlst`), consolidating all metrics into a unified summary table.

---

## Directory Hierarchy

```
staphylococcus-genome-analysis/
├── README.md                  # Project overview and run guide
├── .gitignore                 # Exclusion rules for large genomic files & logs
├── metadata/
│   └── samples.tsv            # Tab-delimited sample sheet (sample_id, strain, fastq paths)
├── scripts/
│   ├── run_qc.sh              # Read trimming & quality control with fastp
│   ├── run_assembly.sh        # De novo genome assembly with SPAdes
│   ├── run_annotation.sh      # Genome annotation with Prokka & MLST typing
│   └── summarize_metrics.py   # Python script consolidating QUAST, BUSCO & Prokka stats
├── results/
│   ├── qc/                    # Cleaned FASTQ files and fastp HTML/JSON reports
│   ├── assemblies/            # SPAdes contigs, scaffolds, and QUAST reports
│   ├── annotations/           # Prokka GFF3, GBK, FFN, FAA, and TXT files
│   ├── mlst_report.tsv        # Multi-Locus Sequence Typing report
│   └── assembly_summary.tsv   # Consolidated multi-sample quality summary table
└── docs/
    └── workflow_methods.md    # Detailed scientific rationale and QC thresholds
```

---

## Tool Prerequisites

Ensure the following bioinformatics tools and libraries are installed. The recommended environment can be set up via **Conda** / **Mamba** from Bioconda and conda-forge:

| Tool | Version | Purpose |
| :--- | :--- | :--- |
| **fastp** | `≥ 0.23` | Ultra-fast adapter trimming, poly-G removal, and quality filtering |
| **SPAdes** | `≥ 3.15` | De novo isolate bacterial genome assembler |
| **QUAST** | `≥ 5.2` | Assembly continuity evaluation (N50, L50, length, GC%) |
| **BUSCO** | `≥ 5.4` | Genome completeness using `staphylococcus_odb10` / `odb12` |
| **Prokka** | `≥ 1.14` | Rapid prokaryotic genome feature annotation |
| **mlst** | `≥ 2.23` | Multi-Locus Sequence Typing from contigs (PubMLST) |
| **Python** | `≥ 3.8` | Metric aggregation with `pandas` |

### Setting Up the Conda Environment

You can install all dependencies into a dedicated conda environment:

```bash
conda create -n staph-wgs -c conda-forge -c bioconda \
    fastp spades quast busco prokka mlst pandas -y

conda activate staph-wgs
```

---

## Quickstart & Execution Guide

### 1. Configure Sample Metadata
Edit [`metadata/samples.tsv`](metadata/samples.tsv) to specify your paired-end FASTQ file locations:

```tsv
sample_id	strain	layout	fastq_1	fastq_2
SA_USA300	FPR3757	paired	data/raw_fastq/SA_USA300_R1.fastq.gz	data/raw_fastq/SA_USA300_R2.fastq.gz
SA_N315	N315	paired	data/raw_fastq/SA_N315_R1.fastq.gz	data/raw_fastq/SA_N315_R2.fastq.gz
SA_COL	COL	paired	data/raw_fastq/SA_COL_R1.fastq.gz	data/raw_fastq/SA_COL_R2.fastq.gz
```

### 2. Quality Control & Adapter Trimming
Run `fastp` to filter low-quality bases (Q < 20), remove poly-G tails, and clip adapter sequences:

```bash
./scripts/run_qc.sh metadata/samples.tsv 8
```
*Outputs: Cleaned FASTQs and HTML reports in `results/qc/`.*

### 3. De Novo Genome Assembly
Execute SPAdes in isolate mode to construct contigs and scaffolds:

```bash
./scripts/run_assembly.sh metadata/samples.tsv 8 16
```
*(Arguments: sample sheet, threads [8], memory in GB [16])*  
*Outputs: Contigs (`*_contigs.fasta`), scaffolds, and QUAST reports in `results/assemblies/`.*

### 4. Genome Annotation & MLST Typing
Annotate assemblies with Prokka (genus *Staphylococcus*) and call sequence types:

```bash
./scripts/run_annotation.sh metadata/samples.tsv 8
```
*Outputs: GFF3/GBK/FFN/FAA files in `results/annotations/` and `results/mlst_report.tsv`.*

### 5. Consolidate Metrics & Summary Report
Extract continuity, completeness (BUSCO), annotation counts, and STs into a single table:

```bash
python3 scripts/summarize_metrics.py
```
*Outputs: [`results/assembly_summary.tsv`](results/assembly_summary.tsv).*

---

## Expected Output Files

| Output Path | Description |
| :--- | :--- |
| `results/qc/{sample}_fastp.html` | Interactive QC report before and after trimming |
| `results/qc/{sample}_R1.trimmed.fastq.gz` | High-quality paired-end sequencing reads |
| `results/assemblies/{sample}_contigs.fasta` | Assembled genomic contigs |
| `results/assemblies/{sample}_quast/report.tsv` | Detailed QUAST assembly continuity statistics |
| `results/annotations/{sample}_prokka/{sample}.gff` | Standard GFF3 functional annotation |
| `results/annotations/{sample}_prokka/{sample}.gbk` | GenBank file for visualization in Artemis / IGV |
| `results/annotations/{sample}_prokka/{sample}.faa` | Translated protein sequences |
| `results/annotations/{sample}_prokka/{sample}.ffn` | Nucleotide sequences of all predicted coding genes |
| `results/mlst_report.tsv` | MLST sequence types (ST) and 7-gene allele calls |
| `results/assembly_summary.tsv` | Consolidated summary across all sequenced isolates |

---

## Documentation & Methodology

For detailed scientific descriptions of algorithms, parameter selection, and quality thresholds (e.g. expected *S. aureus* genome size ~2.8 Mb, GC content ~32.8%, and BUSCO benchmarks), refer to:
📖 [Workflow Methodology & Quality Standards](docs/workflow_methods.md)
