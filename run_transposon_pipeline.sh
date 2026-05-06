#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
#  Transposon Sequencing Pipeline
#  Usage: run_transposon_pipeline.sh <input.fastq> <genome.gbff> [--pseudomonas]
#
#  Examples:
#    bash run_transposon_pipeline.sh sample1.fastq /path/to/genome.gbff
#    bash run_transposon_pipeline.sh sample1.fastq /path/to/28a24.gbff --pseudomonas
#
#  --pseudomonas  Also append a pseudomonas.com URL to the summary output.
#                 Only valid for P. stutzeri 28a24 (repliconid=352254).
#
#  The BLAST database is built automatically from the .gbff file on first run
#  and reused on subsequent runs (stored alongside the .gbff file).
#
#  Expected directory structure:
#    pipeline_dir/
#    ├── run_transposon_pipeline.sh
#    ├── plot_genomic_region.R
#    ├── barcodes.txt                        ← or in scripts_and_barcodes/
#    ├── FastQ-reads_to_barcode_FastA.pl     ← or in scripts_and_barcodes/
#    ├── blast_read_output_with_barcode_fixed.pl  ← or in scripts_and_barcodes/
#    └── genomes/
#        └── (BLAST db files auto-generated here)
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Argument handling ─────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input.fastq> <genome.gbff> [--pseudomonas]"
    echo "  input.fastq    — raw FastQ reads"
    echo "  genome.gbff    — path to GenBank flat file for the reference genome"
    echo "  --pseudomonas  — also output pseudomonas.com URLs (28a24 only)"
    exit 1
fi

FASTQ_INPUT="$1"
GBK_FILE="$2"
PSEUDOMONAS_URLS=false
if [[ "${3:-}" == "--pseudomonas" ]]; then
    PSEUDOMONAS_URLS=true
fi

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Look for scripts in the immediate directory first, then scripts_and_barcodes/
find_script() {
    local filename="$1"
    if [[ -f "$SCRIPT_DIR/$filename" ]]; then
        echo "$SCRIPT_DIR/$filename"
    elif [[ -f "$SCRIPT_DIR/scripts_and_barcodes/$filename" ]]; then
        echo "$SCRIPT_DIR/scripts_and_barcodes/$filename"
    else
        echo ""
    fi
}

FASTQ_TO_FASTA="$(find_script FastQ-reads_to_barcode_FastA.pl)"
BLAST_PARSER="$(find_script blast_read_output_with_barcode_fixed.pl)"
BARCODES="$(find_script barcodes.txt)"
PLOT_SCRIPT="$(find_script plot_genomic_region.R)"

# BLAST database lives next to the .gbff file, named <stem>.blastdb
GBK_DIR="$(cd "$(dirname "$GBK_FILE")" && pwd)"
GBK_STEM="$(basename "$GBK_FILE" | sed 's/\.[^.]*$//')"
BLAST_DB="$GBK_DIR/${GBK_STEM}.blastdb"
FASTA_FOR_DB="$GBK_DIR/${GBK_STEM}.fasta"

# ── Derive output filenames ───────────────────────────────────────────────────
BASENAME="$(basename "$FASTQ_INPUT")"
OUTPUT_DIR="$(pwd)/output"
mkdir -p "$OUTPUT_DIR"
FASTA_OUT="$OUTPUT_DIR/${BASENAME%.fastq}.fasta"
BLAST_OUT="$OUTPUT_DIR/${BASENAME%.fastq}.blast"
FINAL_OUT="$OUTPUT_DIR/${BASENAME}.out"
PLOT_TSV="$OUTPUT_DIR/${BASENAME}.plot_input.tsv"
PLOTS_DIR="$OUTPUT_DIR/plots"

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

if [[ ! -f "$GBK_FILE" ]]; then
    echo "❌ GenBank file not found: $GBK_FILE"
    exit 1
fi

if ! command -v blastn &>/dev/null; then
    echo "❌ blastn not found in PATH — is BLAST+ installed?"
    exit 1
fi

if ! command -v makeblastdb &>/dev/null; then
    echo "❌ makeblastdb not found in PATH — is BLAST+ installed?"
    exit 1
fi

echo ""
echo "  Input FastQ : $FASTQ_INPUT"
echo "  GenBank file: $GBK_FILE"
echo "  BLAST db    : $BLAST_DB"
echo "  Final output: $FINAL_OUT"
echo ""

# ── Step 0: Build BLAST database if needed ────────────────────────────────────
if ls "${BLAST_DB}".* &>/dev/null 2>&1; then
    echo "▶ Step 0/4 — BLAST database already exists, skipping build."
    echo "  ($BLAST_DB)"
else
    echo "▶ Step 0/4 — Building BLAST database from $GBK_FILE..."

    # Extract FASTA from the ORIGIN block using perl (reliable on all macOS versions)
    perl -e '
        my ($locus, $version, $in_seq, $seq) = ("", "", 0, "");
        while (<STDIN>) {
            chomp;
            if (/^LOCUS\s+(\S+)/)   { $locus   = $1 }
            if (/^VERSION\s+(\S+)/) { $version = $1 }
            if (/^ORIGIN/)          { $in_seq  = 1; next }
            if (/^\/\//) {
                if ($in_seq && $seq) {
                    my $id = $version || $locus;
                    print ">$id\n";
                    for (my $i = 0; $i < length($seq); $i += 60) {
                        print substr($seq, $i, 60) . "\n";
                    }
                }
                $in_seq = 0; $seq = ""; $version = "";
                next;
            }
            if ($in_seq) {
                s/^[\s0-9]+//;
                s/\s//g;
                $seq .= uc($_);
            }
        }
    ' < "$GBK_FILE" > "$FASTA_FOR_DB"

    FASTA_SEQS=$(grep -c "^>" "$FASTA_FOR_DB" 2>/dev/null || true)
    FASTA_SEQS=$(echo "$FASTA_SEQS" | tr -d '[:space:]')
    if [[ -z "$FASTA_SEQS" || "$FASTA_SEQS" -eq 0 ]]; then
        echo "❌ Failed to extract FASTA from $GBK_FILE"
        rm -f "$FASTA_FOR_DB"
        exit 1
    fi
    echo "  Extracted $FASTA_SEQS sequence(s) to $FASTA_FOR_DB"

    makeblastdb \
        -in     "$FASTA_FOR_DB" \
        -dbtype nucl \
        -out    "$BLAST_DB" \
        -title  "$GBK_STEM"
    echo "  ✅ BLAST database built: $BLAST_DB"
fi
echo ""

# ── Step 1: FastQ → barcode FASTA ─────────────────────────────────────────────
echo "▶ Step 1/4 — Converting FastQ to barcode FASTA..."
perl "$FASTQ_TO_FASTA" "$FASTQ_INPUT" "$BARCODES" "$FASTA_OUT"

FASTA_READS=$(grep -c "^>" "$FASTA_OUT" 2>/dev/null || echo 0)
FASTA_READS=$(echo "$FASTA_READS" | tr -d '[:space:]')
echo "  ✅ $FASTA_READS sequences written to $FASTA_OUT"
echo ""

if [[ "$FASTA_READS" -eq 0 ]]; then
    echo "⚠️  No sequences output — check that your FastQ contains expected barcode flanking sequences."
    exit 1
fi

# ── Step 2: BLAST ─────────────────────────────────────────────────────────────
echo "▶ Step 2/4 — Running BLAST against $GBK_STEM..."
blastn \
    -query          "$FASTA_OUT" \
    -db             "$BLAST_DB" \
    -out            "$BLAST_OUT" \
    -outfmt         6 \
    -perc_identity  95 \
    -num_threads    4 \
    -max_target_seqs 1

BLAST_HITS=$(wc -l < "$BLAST_OUT" | tr -d '[:space:]')
echo "  ✅ $BLAST_HITS BLAST hits written to $BLAST_OUT"
echo ""

if [[ "$BLAST_HITS" -eq 0 ]]; then
    echo "⚠️  No BLAST hits found — check your database and sequences."
    exit 1
fi

# ── Step 3: Parse BLAST output ────────────────────────────────────────────────
echo "▶ Step 3/4 — Parsing BLAST output..."
PSEUDO_FLAG=""
if [[ "$PSEUDOMONAS_URLS" == "true" ]]; then
    PSEUDO_FLAG="--pseudomonas"
fi
perl "$BLAST_PARSER" "$BLAST_OUT" "$FINAL_OUT" "$PLOT_TSV" $PSEUDO_FLAG

FINAL_LINES=$(wc -l < "$FINAL_OUT" | tr -d '[:space:]')
echo "  ✅ $FINAL_LINES location entries written to $FINAL_OUT"
echo ""

# ── Step 4: Generate genomic context plots ────────────────────────────────────
echo "▶ Step 4/4 — Generating genomic context plots..."

if [[ -z "$PLOT_SCRIPT" ]]; then
    echo "  ⚠️  plot_genomic_region.R not found — skipping plots."
elif [[ ! -f "$PLOT_TSV" ]]; then
    echo "  ⚠️  Plot input TSV not found ($PLOT_TSV) — skipping plots."
else
    mkdir -p "$PLOTS_DIR"

    tail -n +2 "$PLOT_TSV" | while IFS=$'\t' read -r label contig start stop count; do
        PLOT_OUT="$PLOTS_DIR/${label}.png"
        echo "  Plotting $label  ($contig : $start-$stop, n=$count reads)..."
        Rscript "$PLOT_SCRIPT" \
            --gbk    "$GBK_FILE" \
            --contig "$contig" \
            --start  "$start" \
            --end    "$stop" \
            --label  "Transposon Insertion Point" \
            --out    "$PLOT_OUT" \
            2>&1 | grep -v "^##" \
            || echo "    ⚠️  Plot failed for $label"
    done

    echo "  ✅ Plots written to $PLOTS_DIR/"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Pipeline complete!"
echo "  Final output : $FINAL_OUT ($FINAL_LINES entries)"
echo "  Plot input   : $PLOT_TSV"
echo "  Plots        : $PLOTS_DIR/"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
