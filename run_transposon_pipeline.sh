#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
#  Transposon Sequencing Pipeline
#  Usage: run_transposon_pipeline.sh <input.fastq> <blast_db_name>
#
#  Example:
#    ./run_transposon_pipeline.sh sample1.fastq PA14
#
#  Expected directory structure:
#    pipeline_dir/
#    ├── run_transposon_pipeline.sh       ← this script
#    ├── FastQ-reads_to_barcode_FastA.pl
#    ├── blast_read_output_with_barcode_fixed.pl
#    ├── barcodes.txt
#    └── genomes/
#        └── PA14.*  (BLAST database files)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Argument handling ─────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input.fastq> <blast_db_name>"
    echo "  input.fastq    — raw FastQ reads"
    echo "  blast_db_name  — name of BLAST database in ./genomes/ (e.g. PA14)"
    exit 1
fi

FASTQ_INPUT="$1"
BLAST_DB_NAME="$2"

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_AND_BARCODES="$SCRIPT_DIR/scripts_and_barcodes"
FASTQ_TO_FASTA="$SCRIPTS_AND_BARCODES/FastQ-reads_to_barcode_FastA.pl"
BLAST_PARSER="$SCRIPTS_AND_BARCODES/blast_read_output_with_barcode_fixed.pl"
BARCODES="$SCRIPTS_AND_BARCODES/barcodes.txt"
BLAST_DB="$SCRIPT_DIR/genomes/$BLAST_DB_NAME"

# ── Derive output filenames from input fastq name ─────────────────────────────
BASENAME="$(basename "$FASTQ_INPUT")"          # e.g. sample1.fastq
OUTPUT_DIR="$(pwd)/output"                     # output/ in current working directory
mkdir -p "$OUTPUT_DIR"
FASTA_OUT="$OUTPUT_DIR/${BASENAME%.fastq}.fasta"   # e.g. output/sample1.fasta
BLAST_OUT="$OUTPUT_DIR/${BASENAME%.fastq}.blast"   # e.g. output/sample1.blast
FINAL_OUT="$OUTPUT_DIR/${BASENAME}.out"            # e.g. output/sample1.fastq.out

# ── Preflight checks ──────────────────────────────────────────────────────────
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Transposon Sequencing Pipeline"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

for f in "$FASTQ_TO_FASTA" "$BLAST_PARSER" "$BARCODES"; do
    if [[ ! -f "$f" ]]; then
        echo "❌ Missing required file: $f"
        exit 1
    fi
done

if [[ ! -f "$FASTQ_INPUT" ]]; then
    echo "❌ Input FastQ not found: $FASTQ_INPUT"
    exit 1
fi

# Check BLAST db exists (at least one file with the db prefix)
if ! ls "${BLAST_DB}".* &>/dev/null; then
    echo "❌ BLAST database not found: $BLAST_DB"
    echo "   Looking in: $SCRIPT_DIR/genomes/"
    echo "   Available files there:"
    ls "$SCRIPT_DIR/genomes/" 2>/dev/null || echo "   (genomes/ directory not found)"
    exit 1
fi

if ! command -v blastn &>/dev/null; then
    echo "❌ blastn not found in PATH — is BLAST+ installed?"
    exit 1
fi

echo ""
echo "  Input FastQ : $FASTQ_INPUT"
echo "  BLAST db    : $BLAST_DB"
echo "  Final output: $FINAL_OUT"
echo ""

# ── Step 1: FastQ → barcode FASTA ─────────────────────────────────────────────
echo "▶ Step 1/3 — Converting FastQ to barcode FASTA..."
perl "$FASTQ_TO_FASTA" "$FASTQ_INPUT" "$BARCODES" "$FASTA_OUT"

FASTA_READS=$(grep -c "^>" "$FASTA_OUT" 2>/dev/null || echo 0)
echo "  ✅ $FASTA_READS sequences written to $FASTA_OUT"
echo ""

if [[ "$FASTA_READS" -eq 0 ]]; then
    echo "⚠️  No sequences were output — check that your FastQ contains the expected barcode flanking sequences."
    exit 1
fi

# ── Step 2: BLAST ─────────────────────────────────────────────────────────────
echo "▶ Step 2/3 — Running BLAST against $BLAST_DB_NAME..."
blastn \
    -query "$FASTA_OUT" \
    -db "$BLAST_DB" \
    -out "$BLAST_OUT" \
    -outfmt 6 \
    -perc_identity 95 \
    -num_threads 4 \
    -max_target_seqs 1

BLAST_HITS=$(wc -l < "$BLAST_OUT" | tr -d ' ')
echo "  ✅ $BLAST_HITS BLAST hits written to $BLAST_OUT"
echo ""

if [[ "$BLAST_HITS" -eq 0 ]]; then
    echo "⚠️  No BLAST hits found — check your database and sequences."
    exit 1
fi

# ── Step 3: Parse BLAST output ────────────────────────────────────────────────
echo "▶ Step 3/3 — Parsing BLAST output..."
perl "$BLAST_PARSER" "$BLAST_OUT" "$FINAL_OUT"

FINAL_LINES=$(wc -l < "$FINAL_OUT" | tr -d ' ')
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  ✅ Pipeline complete!"
echo "  📄 Final output : $FINAL_OUT ($FINAL_LINES entries)"
echo "  🗂  Intermediates: output/$FASTA_OUT, output/$BLAST_OUT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
