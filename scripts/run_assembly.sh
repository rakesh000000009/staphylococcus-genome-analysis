#!/usr/bin/env bash
# ==============================================================================
# Script: run_assembly.sh
# Purpose: De novo bacterial genome assembly with SPAdes (--isolate mode)
# Usage: ./scripts/run_assembly.sh [sample_sheet] [threads] [memory_gb]
# ==============================================================================

set -euo pipefail

# Configurations & Defaults
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"
SAMPLE_SHEET="${1:-"${PROJECT_ROOT}/metadata/samples.tsv"}"
THREADS="${2:-8}"
MEMORY="${3:-16}"
QC_DIR="${PROJECT_ROOT}/results/qc"
ASSEMBLY_OUTDIR="${PROJECT_ROOT}/results/assemblies"

mkdir -p "${ASSEMBLY_OUTDIR}"

echo "=========================================================="
echo " Starting De Novo Genome Assembly Pipeline (SPAdes)"
echo " Sample sheet : ${SAMPLE_SHEET}"
echo " Assembly dir : ${ASSEMBLY_OUTDIR}"
echo " Threads      : ${THREADS} | Memory: ${MEMORY} GB"
echo "=========================================================="

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
    echo "[-] Error: Sample sheet not found: ${SAMPLE_SHEET}" >&2
    exit 1
fi

# Verify SPAdes is installed
if ! command -v spades.py &>/dev/null; then
    echo "[-] Error: 'spades.py' command not found. Please install SPAdes or activate bioconda environment." >&2
    exit 1
fi

tail -n +2 "${SAMPLE_SHEET}" | while IFS=$'\t' read -r sample_id strain layout fq1 fq2; do
    [[ -z "${sample_id}" ]] && continue

    echo ""
    echo "[+] Assembling sample: ${sample_id} (Strain: ${strain})"

    # Prioritize cleaned/trimmed reads from fastp; fallback to raw reads
    trimmed_r1="${QC_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    trimmed_r2="${QC_DIR}/${sample_id}_R2.trimmed.fastq.gz"

    if [[ -f "${trimmed_r1}" && -f "${trimmed_r2}" ]]; then
        in_r1="${trimmed_r1}"
        in_r2="${trimmed_r2}"
        echo "    Using trimmed FASTQs from results/qc/"
    else
        [[ "${fq1}" != /* ]] && fq1="${PROJECT_ROOT}/${fq1}"
        [[ "${fq2}" != /* ]] && fq2="${PROJECT_ROOT}/${fq2}"
        if [[ -f "${fq1}" && -f "${fq2}" ]]; then
            in_r1="${fq1}"
            in_r2="${fq2}"
            echo "    [!] Trimmed reads not found. Falling back to raw FASTQs."
        else
            echo "    [-] Error: Neither trimmed nor raw reads found for ${sample_id}. Skipping."
            continue
        fi
    fi

    sample_workdir="${ASSEMBLY_OUTDIR}/${sample_id}_spades"
    final_contigs="${ASSEMBLY_OUTDIR}/${sample_id}_contigs.fasta"
    final_scaffolds="${ASSEMBLY_OUTDIR}/${sample_id}_scaffolds.fasta"

    echo "    Running SPAdes in isolate mode..."
    spades.py \
        --isolate \
        -1 "${in_r1}" \
        -2 "${in_r2}" \
        -o "${sample_workdir}" \
        -t "${THREADS}" \
        -m "${MEMORY}"

    # Extract primary contigs and scaffolds
    if [[ -f "${sample_workdir}/contigs.fasta" ]]; then
        cp "${sample_workdir}/contigs.fasta" "${final_contigs}"
        cp "${sample_workdir}/scaffolds.fasta" "${final_scaffolds}"
        echo "    [✓] Assembly finished for ${sample_id}."
        echo "        Contigs:   ${final_contigs}"
        echo "        Scaffolds: ${final_scaffolds}"

        # Run QUAST evaluation if installed
        if command -v quast.py &>/dev/null || command -v quast &>/dev/null; then
            QUAST_CMD=$(command -v quast.py || command -v quast)
            quast_out="${ASSEMBLY_OUTDIR}/${sample_id}_quast"
            echo "    Evaluating assembly continuity with QUAST..."
            "${QUAST_CMD}" "${final_contigs}" \
                -o "${quast_out}" \
                --min-contig 500 \
                --threads "${THREADS}" \
                --labels "${sample_id}"
            echo "    [✓] QUAST report generated in: ${quast_out}"
        fi
    else
        echo "    [-] Error: SPAdes did not produce contigs.fasta for ${sample_id}."
    fi
done

echo ""
echo "=========================================================="
echo " De novo assembly stage finished."
echo "=========================================================="
