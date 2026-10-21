# Staphylococcus Genome Analysis

[![Bioinformatics - WGS](https://img.shields.io/badge/Bioinformatics-Galaxy%20WGS%20Analysis-blue.svg)](#)
[![Target Organism](https://img.shields.io/badge/Organism-Staphylococcus%20aureus-green.svg)](#)
[![Course](https://img.shields.io/badge/Course-AuroSeqTech%20Certificate-orange.svg)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Galaxy-based *de novo* whole-genome sequencing (WGS) analysis of a single *Staphylococcus aureus* isolate (`DRR610201`), completed as a learning project during the 3-day certificate course **"Whole De Novo Genome Analysis of Bacteria Using Galaxy Platform"**, conducted by AuroSeqTech, India (21–23 April 2026).

The analysis was performed interactively on [Galaxy Australia](https://usegalaxy.org.au). The complete published history is accessible here:  
🔗 **Published Galaxy History:** [bacterial-wgs (usegalaxy.org.au)](https://usegalaxy.org.au/u/rakeshsahoo/h/bacterial-wgs)

---

## My Contribution

* **Galaxy Analysis Execution:** Ran the Galaxy workflow tools interactively on Galaxy Australia, tracking intermediate data states and outputs for isolate `DRR610201`.
* **Output Interpretation:** Evaluated and interpreted assembly continuity statistics, biological gene completeness metrics (BUSCO), taxonomic classification, and functional annotation counts.
* **Documentation and Scripting:** Used an AI coding assistant (Antigravity) to draft the repository documentation and command-line helper scripts in `scripts/`.
* **Verification:** Performed local syntax and execution checks on all drafted scripts.

---

## Galaxy Analysis Workflow & Tools

The core analysis tools executed in Galaxy for this project comprise:

| Tool | Role in Galaxy Analysis |
| :--- | :--- |
| **SPAdes** | *De novo* assembly of paired-end Illumina short reads (`DRR610201`) into contigs and scaffolds. |
| **Bowtie2** | Alignment and mapping of paired-end reads back to contigs to assess read mapping consistency and depth. |
| **BUSCO** | Quantitative evaluation of biological completeness against the `staphylococcus_odb12` lineage dataset. |
| **Proksee** | Circular genome visualization, contig arrangement, and genomic feature mapping. |
| **TYGS** | Type (Strain) Genome Server analysis for digital DNA-DNA hybridization (dDDH) and phylogenomic species confirmation against type strains. |

*(Note: The published Galaxy history also includes datasets from FastQC, fastp, Kraken2, QUAST, Prokka, and MLST, which generated the preliminary metrics reported in `results/`.)*

---

## Summary Results for Isolate DRR610201

The table below summarizes measured assembly and annotation metrics for *Staphylococcus aureus* isolate **`DRR610201`** (NCBI SRA accession: `DRR610201`, paired-end Illumina dataset):

| Metric | Tool | Measured Value | Interpretation |
| :--- | :--- | :--- | :--- |
| **Taxonomic Identification** | Kraken2 | 89.30% *S. aureus* (Kraken2) | High proportion of reads assigned to *S. aureus*, supporting isolate identity; remainder represents minor unclassified or background taxa. |
| **Sequence Type (MLST)** | mlst | **ST8** (`3-3-1-1-4-4-3`) | 7-gene housekeeping allelic profile corresponds unambiguously to Sequence Type 8 (associated with the USA300 CA-MRSA lineage). |
| **Total Assembly Length** | SPAdes / QUAST | **3,063,262 bp** (~3.06 Mb) | Total assembly size is slightly higher than typical *S. aureus* genomes (~2.75–2.95 Mb); indicates possible small contig artifacts or minor contamination requiring screening. |
| **Contig Count** | SPAdes / QUAST | 895 contigs | Indicates a fragmented draft assembly with many short contigs, characteristic of short-read *de novo* assembly without long-read scaffolding. |
| **N50 Contig Size** | SPAdes / QUAST | **75,439 bp** (Scaffold N50: 106 kb) | Half of the assembled genome length is contained within contigs of 75,439 bp or larger (scaffolds ≥ 106 kb), showing reasonable contiguity for major chromosomal segments despite fragmentation. |
| **GC Content** | SPAdes / QUAST | **33.23%** | Consistent with the low-GC character of *Staphylococcus* genomes (general reference range ~32.7%–33.0%; QUAST reported 32.75% for contigs ≥ 500 bp). |
| **BUSCO Completeness** | BUSCO | **98.9%** (1,303 / 1,317 orthologs) | High biological completeness against `staphylococcus_odb12` (1,300 single-copy, 3 duplicated, 2 fragmented, 12 missing), showing core conserved genes are largely captured. |
| **Annotated Features** | Prokka | **2,686 CDS**, 62 tRNAs, 5 rRNAs | Total predicted coding sequences and RNA features align with expected gene counts for *S. aureus* (~2,500–2,800 CDS). |

---

## Limitations

* **Single Isolate:** The analysis examines only a single isolate dataset (`DRR610201`). Findings cannot be generalized to broader *Staphylococcus* populations without comparative multi-isolate cohorts.
* **Fragmented Draft Assembly:** With 895 contigs, the assembly remains fragmented. Repetitive regions (such as ribosomal RNA operons, insertion sequences, and prophages) cannot be fully resolved using short reads alone.
* **No Experimental Replicates or Reference-Guided Closure:** The dataset lacks biological/technical replicates. Long-read sequencing (e.g., Oxford Nanopore or PacBio) was not available to scaffold contigs into a finished, circular chromosome.
* **Elevated Assembly Length:** At ~3.06 Mb, the assembled length is slightly above the typical ~2.75–2.95 Mb range for *S. aureus*, suggesting the presence of small fragmented contigs, repetitive expansions, or minor assembly artifacts.
* **Galaxy Execution vs. Command-Line Scripts:** The primary results were generated through the Galaxy web interface; the scripts in `scripts/` are AI-assisted command-line drafts created for educational and reproducibility purposes, rather than the original execution environment.

---

## Directory Hierarchy

```
staphylococcus-genome-analysis/
├── README.md                  # Project overview, results, and limitations
├── LICENSE                    # MIT License (Rakesh Sahoo, 2026)
├── SUGGESTED_GITHUB_SETTINGS.md # Recommended GitHub repository description and topics
├── .gitignore                 # Exclusion rules for large genomic files, reads & logs
├── metadata/
│   └── samples.tsv            # Tab-delimited sample sheet (DRR610201)
├── scripts/
│   ├── run_qc.sh              # Read trimming with fastp (CLI draft)
│   ├── run_assembly.sh        # De novo genome assembly with SPAdes (CLI draft)
│   ├── run_annotation.sh      # Genome annotation with Prokka & MLST (CLI draft)
│   └── summarize_metrics.py   # Metric aggregator script (CLI draft)
├── results/
│   ├── qc/                    # Kraken2 taxonomic report for DRR610201
│   ├── assemblies/            # BUSCO completeness report for DRR610201
│   ├── annotations/           # Prokka summary report for DRR610201
│   ├── mlst_report.tsv        # Multi-Locus Sequence Typing report
│   └── assembly_summary.tsv   # Consolidated metrics table
└── docs/
    └── workflow_methods.md    # Detailed methodology and reference ranges
```

---

## Optional Command-Line Scripts (AI-Assisted Draft, Run Locally [date: [will update soon]])

The scripts in `scripts/` provide an alternative command-line implementation of read QC, assembly, annotation, and metric aggregation. They were drafted with AI assistance (Antigravity) and tested locally.

### Sample Metadata
The sample sheet [`metadata/samples.tsv`](metadata/samples.tsv) lists the active dataset:

```tsv
sample_id	strain	layout	fastq_1	fastq_2
DRR610201	USA300	paired	data/raw_fastq/DRR610201_1.fastq.gz	data/raw_fastq/DRR610201_2.fastq.gz
```

> **Note on Planned Samples:** Additional reference strains (`SA_USA300`, `SA_N315`, `SA_COL`) were initially considered for comparative benchmarking, but were not processed during the certificate course. They may be included in future comparative analyses.

### Conda Environment Setup (Optional)
To test or execute the command-line scripts locally:

```bash
conda create -n staph-wgs -c conda-forge -c bioconda \
    fastp spades quast busco prokka mlst pandas -y

conda activate staph-wgs
```

### Script Execution Overview

1. **Quality Control & Read Trimming:**
   ```bash
   ./scripts/run_qc.sh metadata/samples.tsv 4
   ```
2. **De Novo Assembly & Continuity Check:**
   ```bash
   ./scripts/run_assembly.sh metadata/samples.tsv 8 16
   ```
3. **Genome Annotation & MLST Typing:**
   ```bash
   ./scripts/run_annotation.sh metadata/samples.tsv 4
   ```
4. **Aggregate Summary Metrics:**
   ```bash
   python3 scripts/summarize_metrics.py
   ```

---

## Documentation & Methodology

For detailed scientific descriptions of algorithms, Galaxy tools, and general reference ranges, refer to:  
📖 [Scientific Methodology & Workflow Documentation](docs/workflow_methods.md)
