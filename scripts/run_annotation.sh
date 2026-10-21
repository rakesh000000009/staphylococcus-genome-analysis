#!/usr/bin/env bash
# ==============================================================================
# Script: run_annotation.sh
# Purpose: Annotate Staphylococcus genomes with Prokka and perform MLST typing
# Usage: ./scripts/run_annotation.sh [sample_sheet] [threads]
# ==============================================================================

set -euo pipefail

# Configurations & Defaults
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
SAMPLE_SHEET="${1:-"${PROJECT_ROOT}/metadata/samples.tsv"}"
THREADS="${2:-4}"
ASSEMBLY_DIR="${PROJECT_ROOT}/results/assemblies"
ANNOT_OUTDIR="${PROJECT_ROOT}/results/annotations"
MLST_REPORT="${PROJECT_ROOT}/results/mlst_report.tsv"

mkdir -p "${ANNOT_OUTDIR}"

echo "=========================================================="
echo " Starting Staphylococcus Annotation & Typing Pipeline"
echo " Contigs dir   : ${ASSEMBLY_DIR}"
echo " Annotations   : ${ANNOT_OUTDIR}"
echo " MLST report   : ${MLST_REPORT}"
echo " Threads       : ${THREADS}"
echo "=========================================================="

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
    echo "[-] Error: Sample sheet not found: ${SAMPLE_SHEET}" >&2
    exit 1
fi

# Verify Prokka is installed
if ! command -v prokka &>/dev/null; then
    echo "[-] Error: 'prokka' command not found. Please activate your bioconda environment." >&2
    exit 1
fi

# Annotate each sample using Prokka
tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r sample_id strain layout fq1 fq2; do
    [[ -z "${sample_id}" ]] && continue

    contigs_file="${ASSEMBLY_DIR}/${sample_id}_contigs.fasta"
    sample_annot_dir="${ANNOT_OUTDIR}/${sample_id}_prokka"

    echo ""
    echo "[+] Annotating sample: ${sample_id} (Strain: ${strain})"

    if [[ ! -f "${contigs_file}" ]]; then
        echo "    [!] Warning: Contigs file not found: ${contigs_file}"
        echo "        Skipping Prokka annotation for ${sample_id}."
        continue
    fi

    echo "    Running Prokka with genus 'Staphylococcus'..."
    prokka \
        --outdir "${sample_annot_dir}" \
        --prefix "${sample_id}" \
        --genus "Staphylococcus" \
        --species "aureus" \
        --strain "${strain}" \
        --locustag "${sample_id}" \
        --cpus "${THREADS}" \
        --compliant \
        --force \
        "${contigs_file}"

    echo "    [✓] Prokka annotation complete for ${sample_id}:"
    echo "        GFF3: ${sample_annot_dir}/${sample_id}.gff"
    echo "        GBK:  ${sample_annot_dir}/${sample_id}.gbk"
    echo "        FFN:  ${sample_annot_dir}/${sample_id}.ffn"
done

# Perform MLST Typing across all assembled contigs
echo ""
echo "=========================================================="
echo " Running Multi-Locus Sequence Typing (MLST)"
echo "=========================================================="

if command -v mlst &>/dev/null; then
    # Find all contigs
    shopt -s nullglob
    contig_files=("${ASSEMBLY_DIR}"/*_contigs.fasta)
    shopt -u nullglob

    if [[ ${#contig_files[@]} -gt 0 ]]; then
        echo "[+] Found ${#contig_files[@]} assembly contig file(s). Calling MLST..."
        mlst --scheme saureus "${contig_files[@]}" > "${MLST_REPORT}" || mlst "${contig_files[@]}" > "${MLST_REPORT}"
        echo "[✓] MLST typing completed. Report written to: ${MLST_REPORT}"
        echo "--- MLST Preview ---"
        cat "${MLST_REPORT}"
    else
        echo "[!] No contig files found matching '${ASSEMBLY_DIR}/*_contigs.fasta'."
        echo "    Writing template MLST report to ${MLST_REPORT}."
        printf "FILE\tSCHEME\tST\tarcC\taroE\tglpF\tgmk\tpta\ttpi\tyqiL\n" > "${MLST_REPORT}"
    fi
else
    echo "[!] Warning: 'mlst' command not found. Writing template MLST report."
    printf "FILE\tSCHEME\tST\tarcC\taroE\tglpF\tgmk\tpta\ttpi\tyqiL\n" > "${MLST_REPORT}"
fi

echo ""
echo "=========================================================="
echo " Annotation & MLST stage completed successfully."
echo "=========================================================="
