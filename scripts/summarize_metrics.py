#!/usr/bin/env python3
# Drafted with AI assistance; not used to generate the results in results/ unless stated in README.
"""
Summarize Assembly Continuity and Completeness Metrics
======================================================
Usage and requirements:
  Requirements:
    - Python 3.8+
    - pandas (pip install pandas or conda install pandas)
  Usage:
    python3 scripts/summarize_metrics.py [options]
  Options:
    --metadata PATH         Path to sample sheet (default: metadata/samples.tsv)
    --assemblies-dir PATH   Path to assembly directory (default: results/assemblies)
    --annotations-dir PATH  Path to annotations directory (default: results/annotations)
    --mlst-report PATH      Path to MLST report (default: results/mlst_report.tsv)
    --output PATH           Destination summary TSV (default: results/assembly_summary.tsv)
  Inputs:
    Assembly contigs, QUAST report directories, BUSCO short summaries,
    Prokka summary .txt files, and MLST tabular reports.
  Outputs:
    Consolidated summary table written to results/assembly_summary.tsv

This script parses assembly quality metrics (QUAST), gene completeness (BUSCO),
functional annotation counts (Prokka), and sequence types (MLST) across samples,
producing a consolidated tabular summary.
"""

import os
import sys
import glob
import json
import re
import argparse
from pathlib import Path
from typing import Dict, Any, Optional
import pandas as pd


def compute_fasta_stats(fasta_path: Path) -> Dict[str, Any]:
    """Calculate basic assembly metrics directly from a FASTA file as fallback."""
    lengths = []
    gc_count = 0
    total_bases = 0

    if not fasta_path.exists():
        return {}

    current_seq = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    seq_str = "".join(current_seq).upper()
                    lengths.append(len(seq_str))
                    gc_count += seq_str.count("G") + seq_str.count("C")
                    total_bases += len(seq_str)
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            seq_str = "".join(current_seq).upper()
            lengths.append(len(seq_str))
            gc_count += seq_str.count("G") + seq_str.count("C")
            total_bases += len(seq_str)

    if not lengths:
        return {}

    lengths.sort(reverse=True)
    half_len = total_bases / 2.0
    running_sum = 0
    n50 = 0
    l50 = 0
    for idx, l in enumerate(lengths, 1):
        running_sum += l
        if running_sum >= half_len:
            n50 = l
            l50 = idx
            break

    gc_percent = round((gc_count / total_bases) * 100, 2) if total_bases > 0 else 0.0

    return {
        "Contigs_Count": len(lengths),
        "Total_Length_bp": total_bases,
        "N50_bp": n50,
        "L50": l50,
        "GC_Percent": gc_percent,
    }


def parse_quast_report(quast_dir: Path) -> Dict[str, Any]:
    """Extract N50, total length, contig counts, and GC% from QUAST report."""
    metrics = {}
    report_tsv = quast_dir / "report.tsv"
    transposed_tsv = quast_dir / "transposed_report.tsv"

    if report_tsv.exists():
        with open(report_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    key = parts[0].strip()
                    val = parts[1].strip()
                    if key in ["# contigs", "# contigs (>= 0 bp)", "# contigs (>= 500 bp)"]:
                        metrics["Contigs_Count"] = int(val)
                    elif key in ["Total length", "Total length (>= 0 bp)", "Total length (>= 500 bp)"]:
                        metrics["Total_Length_bp"] = int(val)
                    elif key == "N50":
                        metrics["N50_bp"] = int(val)
                    elif key == "L50":
                        metrics["L50"] = int(val)
                    elif key == "GC (%)":
                        metrics["GC_Percent"] = float(val)
        return metrics

    if transposed_tsv.exists():
        df = pd.read_csv(transposed_tsv, sep="\t")
        if not df.empty:
            row = df.iloc[0]
            for col in df.columns:
                if "# contigs" in col:
                    metrics["Contigs_Count"] = int(row[col])
                elif "Total length" in col:
                    metrics["Total_Length_bp"] = int(row[col])
                elif col == "N50":
                    metrics["N50_bp"] = int(row[col])
                elif col == "L50":
                    metrics["L50"] = int(row[col])
                elif "GC" in col:
                    metrics["GC_Percent"] = float(row[col])
    return metrics


def parse_busco_metrics(busco_search_dir: Path, sample_id: str) -> Dict[str, Any]:
    """Parse BUSCO summary files (JSON or TXT) for gene completeness percentages."""
    metrics = {
        "BUSCO_Complete_%": "NA",
        "BUSCO_Single_%": "NA",
        "BUSCO_Duplicated_%": "NA",
        "BUSCO_Fragmented_%": "NA",
        "BUSCO_Missing_%": "NA",
    }

    # Search for json reports first
    json_files = list(busco_search_dir.glob(f"**/*{sample_id}*busco*/**/*.json")) + \
                 list(busco_search_dir.glob(f"**/*busco*{sample_id}*.json")) + \
                 list(busco_search_dir.glob(f"**/short_summary.*{sample_id}*.json"))

    for jpath in json_files:
        try:
            with open(jpath, "r", encoding="utf-8") as jf:
                data = json.load(jf)
                res = data.get("results", {})
                if "Complete percentage" in res:
                    metrics["BUSCO_Complete_%"] = res.get("Complete percentage")
                    metrics["BUSCO_Single_%"] = res.get("Single copy percentage")
                    metrics["BUSCO_Duplicated_%"] = res.get("Multi copy percentage")
                    metrics["BUSCO_Fragmented_%"] = res.get("Fragmented percentage")
                    metrics["BUSCO_Missing_%"] = res.get("Missing percentage")
                    return metrics
        except Exception:
            continue

    # Fallback to txt reports
    txt_files = list(busco_search_dir.glob(f"**/short_summary.*{sample_id}*.txt"))
    for tpath in txt_files:
        with open(tpath, "r", encoding="utf-8", errors="ignore") as tf:
            content = tf.read()
            # Match standard pattern: C:98.5%[S:98.0%,D:0.5%],F:0.5%,M:1.0%
            match = re.search(r"C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%", content)
            if match:
                metrics["BUSCO_Complete_%"] = float(match.group(1))
                metrics["BUSCO_Single_%"] = float(match.group(2))
                metrics["BUSCO_Duplicated_%"] = float(match.group(3))
                metrics["BUSCO_Fragmented_%"] = float(match.group(4))
                metrics["BUSCO_Missing_%"] = float(match.group(5))
                return metrics
    return metrics


def parse_prokka_summary(annot_dir: Path, sample_id: str) -> Dict[str, Any]:
    """Extract CDS, tRNA, and rRNA counts from Prokka txt report."""
    metrics = {"Prokka_CDS": "NA", "Prokka_tRNA": "NA", "Prokka_rRNA": "NA"}
    txt_path = annot_dir / f"{sample_id}_prokka" / f"{sample_id}.txt"
    if not txt_path.exists():
        txt_files = list(annot_dir.glob(f"**/{sample_id}.txt"))
        if txt_files:
            txt_path = txt_files[0]
        else:
            return metrics

    with open(txt_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = [p.strip() for p in line.split(":", 1)]
            if len(parts) == 2:
                key, val = parts[0], parts[1]
                if key == "CDS":
                    metrics["Prokka_CDS"] = int(val)
                elif key == "tRNA":
                    metrics["Prokka_tRNA"] = int(val)
                elif key == "rRNA":
                    metrics["Prokka_rRNA"] = int(val)
    return metrics


def parse_mlst_summary(mlst_tsv: Path) -> Dict[str, str]:
    """Parse Torsten Seemann's MLST output into a dict {sample_id: ST}."""
    st_dict = {}
    if not mlst_tsv.exists():
        return st_dict

    with open(mlst_tsv, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("FILE"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                filename = Path(parts[0]).name
                # Extract sample id (e.g. SA_USA300 from SA_USA300_contigs.fasta)
                sid = filename.replace("_contigs.fasta", "").replace(".fasta", "")
                st = parts[2]
                st_dict[sid] = st
    return st_dict


def main():
    parser = argparse.ArgumentParser(
        description="Consolidate Staphylococcus assembly, annotation, and completeness metrics."
    )
    parser.add_argument(
        "--metadata",
        default="metadata/samples.tsv",
        help="Path to sample sheet (default: metadata/samples.tsv)",
    )
    parser.add_argument(
        "--assemblies-dir",
        default="results/assemblies",
        help="Path to assembly directory (default: results/assemblies)",
    )
    parser.add_argument(
        "--annotations-dir",
        default="results/annotations",
        help="Path to annotations directory (default: results/annotations)",
    )
    parser.add_argument(
        "--mlst-report",
        default="results/mlst_report.tsv",
        help="Path to MLST report (default: results/mlst_report.tsv)",
    )
    parser.add_argument(
        "--output",
        default="results/assembly_summary.tsv",
        help="Destination summary TSV (default: results/assembly_summary.tsv)",
    )
    args = parser.parse_args()

    meta_path = Path(args.metadata)
    asm_dir = Path(args.assemblies_dir)
    annot_dir = Path(args.annotations_dir)
    mlst_path = Path(args.mlst_report)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect sample list from metadata if available, otherwise scan assembly dir
    samples = []
    sample_strains = {}
    if meta_path.exists():
        meta_df = pd.read_csv(meta_path, sep="\t")
        for _, row in meta_df.iterrows():
            sid = str(row["sample_id"]).strip()
            samples.append(sid)
            sample_strains[sid] = str(row.get("strain", "NA")).strip()
    else:
        for f in asm_dir.glob("*_contigs.fasta"):
            sid = f.name.replace("_contigs.fasta", "")
            samples.append(sid)
            sample_strains[sid] = "Unknown"

    mlst_data = parse_mlst_summary(mlst_path)

    rows = []
    for sid in samples:
        strain = sample_strains.get(sid, "NA")
        fasta_file = asm_dir / f"{sid}_contigs.fasta"
        quast_dir = asm_dir / f"{sid}_quast"

        # 1. Continuity Metrics: Prioritize QUAST; fallback to direct FASTA analysis
        quast_metrics = parse_quast_report(quast_dir)
        if not quast_metrics and fasta_file.exists():
            quast_metrics = compute_fasta_stats(fasta_file)

        # 2. Completeness Metrics (BUSCO)
        busco_metrics = parse_busco_metrics(asm_dir, sid)

        # 3. Annotation Metrics (Prokka)
        prokka_metrics = parse_prokka_summary(annot_dir, sid)

        # 4. MLST
        sequence_type = mlst_data.get(sid, "NA")

        row = {
            "Sample_ID": sid,
            "Strain": strain,
            "Sequence_Type (ST)": sequence_type,
            "Contigs_Count": quast_metrics.get("Contigs_Count", "NA"),
            "Total_Length_bp": quast_metrics.get("Total_Length_bp", "NA"),
            "N50_bp": quast_metrics.get("N50_bp", "NA"),
            "L50": quast_metrics.get("L50", "NA"),
            "GC_Percent": quast_metrics.get("GC_Percent", "NA"),
            "BUSCO_Complete_%": busco_metrics.get("BUSCO_Complete_%", "NA"),
            "BUSCO_Single_%": busco_metrics.get("BUSCO_Single_%", "NA"),
            "BUSCO_Duplicated_%": busco_metrics.get("BUSCO_Duplicated_%", "NA"),
            "BUSCO_Fragmented_%": busco_metrics.get("BUSCO_Fragmented_%", "NA"),
            "BUSCO_Missing_%": busco_metrics.get("BUSCO_Missing_%", "NA"),
            "Prokka_CDS": prokka_metrics.get("Prokka_CDS", "NA"),
            "Prokka_tRNA": prokka_metrics.get("Prokka_tRNA", "NA"),
            "Prokka_rRNA": prokka_metrics.get("Prokka_rRNA", "NA"),
        }
        rows.append(row)

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_path, sep="\t", index=False)

    print("=================================================================")
    print(f"[✓] Successfully generated summary report: {out_path}")
    print("=================================================================")
    print(df_out.to_string(index=False))


if __name__ == "__main__":
    main()
