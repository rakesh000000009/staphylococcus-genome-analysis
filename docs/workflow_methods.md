# Scientific Methodology and Workflow Documentation

This document describes the scientific principles, computational tools, reference ranges, and analytical methodology for the **Staphylococcus Whole-Genome Sequencing (WGS) Analysis**.

The primary analysis was performed within the [Galaxy Australia](https://usegalaxy.org.au) web platform as part of the 3-day certificate course, *"Whole De Novo Genome Analysis of Bacteria Using Galaxy Platform"*, organized by AuroSeqTech, India (21–23 April 2026). Published Galaxy history: [bacterial-wgs](https://usegalaxy.org.au/u/rakeshsahoo/h/bacterial-wgs).

---

## 1. Biological and Genomic Background

*Staphylococcus aureus* is a clinically important Gram-positive bacterium associated with hospital-acquired and community-acquired infections. Whole-Genome Sequencing (WGS) facilitates:
- Identification of species identity and lineage classification.
- Detection of antimicrobial resistance (AMR) markers and mobile genetic elements (MGEs).
- Assessment of genomic completeness and sequence continuity in draft assemblies.

### General Genomic Reference Ranges for *Staphylococcus aureus*
* **Genome Size:** ~2.75 to 2.95 Mbp (megabases)
* **GC Content:** ~32.7% to 33.0% (low-GC organism)
* **Expected Coding Sequences (CDS):** ~2,500 – 2,800 genes
* **rRNA Operons:** 5–6 operons (5S, 16S, 23S)
* **tRNAs:** 55–62 tRNA genes

---

## 2. Galaxy Analysis Workflow

The primary analysis of paired-end sequencing reads from isolate **DRR610201** was carried out on Galaxy Australia using the following core tools:

```
                      DRR610201 Paired-End Reads
                                │
                                ▼
        ┌──────────────────────────────────────────────┐
        │ 1. De Novo Assembly (SPAdes)                 │
        │    • Multi-kmer de Bruijn graph assembly     │
        │    • Generates contigs and scaffolds         │
        └──────────────────────┬───────────────────────┘
                               │
                ┌──────────────┼──────────────┐
                ▼              ▼              ▼
   ┌──────────────────┐ ┌─────────────┐ ┌──────────────────┐
   │ 2. Read Mapping  │ │ 3. Quality  │ │ 4. Visualization │
   │    (Bowtie2)     │ │    (BUSCO)  │ │    (Proksee)     │
   │    • Read mapping│ │    • Core   │ │    • Circular    │
   │      to assembly │ │      genes  │ │      genome map  │
   └──────────────────┘ └─────────────┘ └──────────────────┘
                               │
                               ▼
        ┌──────────────────────────────────────────────┐
        │ 5. Taxonomic Confirmation (TYGS)             │
        │    • Digital DNA-DNA hybridization (dDDH)    │
        │    • Phylogenomic comparison to type strains │
        └──────────────────────────────────────────────┘
```

### Galaxy Tool Details

1. **SPAdes (v3.15+ on Galaxy)**:
   - **Role:** *De novo* assembly of paired-end Illumina reads (`DRR610201_1.fastq.gz`, `DRR610201_2.fastq.gz`).
   - **Parameters:** Default isolate mode / auto-selected k-mer sizes [exact Galaxy parameters: [will update soon]].
   - **Outputs:** Assembled contigs (`DRR610201_contigs.fasta`) and scaffolds.

2. **Bowtie2 (Galaxy)**:
   - **Role:** Re-aligns paired-end sequencing reads back to the assembled contigs to evaluate read representation, alignment consistency, and depth [exact Galaxy parameters: [will update soon]].

3. **BUSCO (Galaxy)**:
   - **Role:** Quantitative evaluation of biological assembly completeness.
   - **Lineage dataset:** `staphylococcus_odb12` (1,317 core ortholog groups).
   - **Outputs:** Summary of complete (single-copy and duplicated), fragmented, and missing orthologs.

4. **Proksee (Galaxy / Web)**:
   - **Role:** Circular genomic map generation and visual inspection of assembly contigs, GC content, and GC skew [exact Galaxy parameters: [will update soon]].

5. **TYGS (Type (Strain) Genome Server)**:
   - **Role:** Phylogenomic classification and digital DNA-DNA hybridization (dDDH) against curated bacterial type strains to verify taxonomic species assignment [exact Galaxy parameters: [will update soon]].

*(Note: The published Galaxy history also includes datasets from FastQC, fastp, Kraken2, QUAST, Prokka, and MLST, which provided additional preliminary metric outputs preserved in `results/`.)*

---

## 3. General Reference Ranges and Assembly Observations

The values below represent published general reference ranges for *Staphylococcus aureus* genomes, alongside the actual observed values for isolate DRR610201. These reference ranges provide biological context rather than rigid pass/fail thresholds.

| Parameter | General Reference Range (*S. aureus*) | Observed Value for DRR610201 | Context & Observations |
| :--- | :--- | :--- | :--- |
| **Total Assembly Size** | ~2.75 – 2.95 Mbp | 3,063,262 bp (~3.06 Mb) | Observed size is slightly above typical *S. aureus* genome size, indicating potential small contig artifacts, repeat expansion, or minor contig contamination. |
| **GC Content (%)** | ~32.7% – 33.0% | 33.23% (32.75% for contigs ≥ 500 bp) | Falls near expected low-GC range for *Staphylococcus aureus*. |
| **N50 Contig Size** | Typically > 50,000 bp (short-read draft) | 75,439 bp (Scaffold N50: 106 kb) | Indicates that half of the assembly is assembled into contigs of 75.4 kb or larger. |
| **Contig Count** | Typically < 100 contigs in contiguous assemblies | 895 contigs | Indicates a fragmented draft assembly; does not meet the continuity of finished or long-read assemblies due to short-read repeat breakages. |
| **BUSCO Completeness** | ≥ 98.0% (`staphylococcus_odb12`) | 98.9% (1,303 / 1,317 orthologs) | Core staphylococcal gene content is well-represented despite physical fragmentation. |
| **MLST Allele Matches** | 7/7 housekeeping loci resolved | ST8 (`3-3-1-1-4-4-3`) | Housekeeping alleles unambiguously match Sequence Type 8 (USA300 CA-MRSA lineage). |
| **Annotated CDS** | ~2,500 – 2,800 genes | 2,686 CDS (Prokka) | Coding sequence count falls within expected gene density for *S. aureus*. |

---

## 4. Optional Command-Line Scripts (AI-Assisted Draft, Run Locally [date: [will update soon]])

The repository includes command-line scripts in `scripts/` drafted with AI assistance (Antigravity) and run locally. These were not used as the primary execution pipeline for the certificate course results (which were run in Galaxy). They are provided for local testing and reproducibility:

1. **`scripts/run_qc.sh`**:
   - Executes `fastp` for quality filtering, adapter trimming, poly-G removal, and read length filtering (Q ≥ 20, min length 50 bp).
2. **`scripts/run_assembly.sh`**:
   - Executes `SPAdes` in `--isolate` mode on paired-end reads, followed by continuity assessment using `QUAST`.
3. **`scripts/run_annotation.sh`**:
   - Executes `Prokka` with the genus setting `--genus Staphylococcus` and calls sequence types using `mlst`.
4. **`scripts/summarize_metrics.py`**:
   - Python script utilizing `pandas` to aggregate QUAST metrics, BUSCO summaries, Prokka output counts, and MLST calls into `results/assembly_summary.tsv`.

---

## 5. References

1. **Bankevich, A., et al. (2012).** SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *Journal of Computational Biology*, 19(5), 455-477.
2. **Langmead, B., & Salzberg, S. L. (2012).** Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357-359.
3. **Manni, M., et al. (2021).** BUSCO Update: Novel features and streamlined workflows. *Molecular Biology and Evolution*, 38(10), 4647-4654.
4. **Grant, J. R., et al. (2023).** Proksee: in-depth characterization and visualization of bacterial genomes. *Nucleic Acids Research*, 51(W1), W484-W492.
5. **Meier-Kolthoff, J. P., & Göker, M. (2019).** TYGS: Type (Strain) Genome Server for genome-based taxonomy. *Nature Communications*, 10, 2182.
6. **Chen, S., et al. (2018).** fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890.
7. **Seemann, T. (2014).** Prokka: rapid prokaryotic genome annotation. *Bioinformatics*, 30(14), 2068-2069.
8. **Gurevich, A., et al. (2013).** QUAST: quality assessment tool for genome assemblies. *Bioinformatics*, 29(8), 1072-1075.
9. **Enright, M. C., et al. (2000).** Multilocus sequence typing of *Staphylococcus aureus* for characterization of methicillin-resistant clones. *Journal of Clinical Microbiology*, 38(3), 1008-1015.
