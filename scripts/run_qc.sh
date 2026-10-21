#!/usr/bin/env bash
# ==============================================================================
# Script: run_qc.sh
# Purpose: Read quality control, adapter trimming, and poly-G filtering with fastp
# Usage: ./scripts/run_qc.sh [sample_sheet] [threads]
# ==============================================================================

set -euo pipefail

# Configurations & Defaults
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
SAMPLE_SHEET="${1:-"${PROJECT_ROOT}/metadata/samples.tsv"}"
THREADS="${2:-4}"
QC_OUTDIR="${PROJECT_ROOT}/results/qc"

mkdir -p "${QC_OUTDIR}"

echo "=========================================================="
echo " Starting Quality Control & Read Trimming Pipeline (fastp)"
echo " Sample sheet : ${SAMPLE_SHEET}"
echo " Output dir   : ${QC_OUTDIR}"
echo " Threads      : ${THREADS}"
echo "=========================================================="

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
    echo "[-] Error: Sample sheet not found: ${SAMPLE_SHEET}" >&2
    exit 1
fi

# Verify fastp is available
if ! command -v fastp &>/dev/null; then
    echo "[-] Error: 'fastp' command not found. Please install or activate bioconda environment." >&2
    exit 1
fi

# Iterate over sample sheet (skipping header)
tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r sample_id strain layout fq1 fq2; do
    # Skip empty lines
    [[ -z "${sample_id}" ]] && continue

    echo ""
    echo "[+] Processing sample: ${sample_id} (Strain: ${strain})"

    # Resolve relative paths against PROJECT_ROOT
    [[ "${fq1}" != /* ]] && fq1="${PROJECT_ROOT}/${fq1}"
    [[ "${fq2}" != /* ]] && fq2="${PROJECT_ROOT}/${fq2}"

    if [[ ! -f "${fq1}" || ! -f "${fq2}" ]]; then
        echo "    [!] Warning: Raw FASTQ file(s) missing for ${sample_id}:"
        echo "        R1: ${fq1}"
        echo "        R2: ${fq2}"
        echo "        Skipping sample. Please place raw FASTQs in the specified location."
        continue
    fi

    out_r1="${QC_OUTDIR}/${sample_id}_R1.trimmed.fastq.gz"
    out_r2="${QC_OUTDIR}/${sample_id}_R2.trimmed.fastq.gz"
    html_report="${QC_OUTDIR}/${sample_id}_fastp.html"
    json_report="${QC_OUTDIR}/${sample_id}_fastp.json"

    echo "    Executing fastp paired-end trimming..."
    fastp \
        --in1 "${fq1}" \
        --in2 "${fq2}" \
        --out1 "${out_r1}" \
        --out2 "${out_r2}" \
        --html "${html_report}" \
        --json "${json_report}" \
        --thread "${THREADS}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --length_required 50 \
        --correction \
        --trim_poly_g

    echo "    [✓] Completed ${sample_id}. Cleaned FASTQs and reports saved to ${QC_OUTDIR}"
done

echo ""
echo "=========================================================="
echo " Quality Control stage finished."
echo "=========================================================="
