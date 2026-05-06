#!/usr/bin/perl -w
use strict;

# Input file and output file
my $input_file  = $ARGV[0] or die "Usage: $0 <input_file> <output_file> [plot.tsv] [--pseudomonas]\n";
my $output_file = $ARGV[1] or die "Usage: $0 <input_file> <output_file> [plot.tsv] [--pseudomonas]\n";

# Optional arguments
my $plot_tsv         = $output_file . ".plot_input.tsv";
my $pseudomonas_urls = 0;

for my $i (2..$#ARGV) {
    if ($ARGV[$i] eq '--pseudomonas') {
        $pseudomonas_urls = 1;
    } elsif ($ARGV[$i] !~ /^-/) {
        $plot_tsv = $ARGV[$i];
    }
}

# Pseudomonas.com replicon ID for P. stutzeri 28a24
my $REPLICON_ID = "352254";

open(my $in,  "<", $input_file)  or die "Cannot open $input_file: $!";
open(my $out, ">", $output_file) or die "Cannot write to $output_file: $!";

# ── Data structures ───────────────────────────────────────────────────────────
# %line_counts : barcode -> { "contig\tstart\tstop" -> count }
# Also keep the full output line for the summary file.
my %line_counts;   # barcode -> { location_key -> count }
my %line_text;     # barcode -> { location_key -> full output line (for summary) }

# ── Parse BLAST outfmt 6 ──────────────────────────────────────────────────────
# Columns: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
while (my $line = <$in>) {
    chomp $line;
    $line =~ s/\r//;
    my @fields = split("\t", $line);
    next if @fields < 10;

    my $direction = $fields[0];   # e.g. "BC01_ATCGATCG_right"
    my $contig    = $fields[1];   # e.g. "CP007441.1"
    my $col9      = $fields[8];   # sstart
    my $col10     = $fields[9];   # send

    # Extract barcode ID and raw barcode sequence
    my ($barcode, $sequence) = $direction =~ /^([^_]+)_(.+)/;
    next unless ($barcode && $sequence);

    # Strip _left / _right suffix from sequence field
    $sequence =~ s/_(left|right)$//;

    my ($start, $stop);

    if ($direction =~ /right/) {
        if    ($col9 < $col10) { $start = $col9 - 50; $stop = $col9;       }
        elsif ($col9 > $col10) { $start = $col9;       $stop = $col9 + 50; }
    }
    elsif ($direction =~ /left/) {
        if    ($col9 > $col10) { $start = $col10 - 50; $stop = $col10;       }
        elsif ($col9 < $col10) { $start = $col10;       $stop = $col10 + 50; }
    }

    next unless (defined $start && defined $stop);

    # Key that uniquely identifies a location (contig + coordinates)
    my $loc_key  = "$contig\t$start\t$stop";

    # Full line for the human-readable summary output
    my $url_col = $pseudomonas_urls
        ? "\thttps://www.pseudomonas.com/feature/intergenic?repliconid=${REPLICON_ID}&start=${start}&stop=${stop}"
        : "";
    my $out_line = "$barcode\t$sequence\t$contig\t$start-$stop$url_col";

    $line_counts{$barcode}{$loc_key}++;
    $line_text{$barcode}{$loc_key} = $out_line;
}
close($in);

# ── Write full summary (all locations, all counts) ────────────────────────────
foreach my $barcode (sort keys %line_counts) {
    foreach my $loc_key (sort keys %{ $line_counts{$barcode} }) {
        my $count = $line_counts{$barcode}{$loc_key};
        print $out "$line_text{$barcode}{$loc_key}\t$count\n";
    }
}
close($out);
print "Output written to $output_file.\n";

# ── Build per-barcode best-hit TSV for R plot script ─────────────────────────
# For each barcode pick the location with the highest read count.
# "unknown" barcodes that share a location with a named barcode are merged;
# those with a unique location get sequential labels: unknown_1, unknown_2, ...

my %best;          # barcode -> { contig, start, stop, count }

foreach my $barcode (sort keys %line_counts) {
    my $locs   = $line_counts{$barcode};
    # Find location key with max count
    my $best_key = (sort { $locs->{$b} <=> $locs->{$a} } keys %$locs)[0];
    my $count    = $locs->{$best_key};
    my ($contig, $start, $stop) = split("\t", $best_key);
    $best{$barcode} = { contig => $contig,
                        start  => $start,
                        stop   => $stop,
                        count  => $count };
}

# Separate named barcodes from unknowns
my %named_locations;  # "contig:start:stop" -> barcode label (for dedup)
my @named_barcodes;
my @unknown_barcodes;

foreach my $bc (sort keys %best) {
    if ($bc =~ /^unknown$/i || $bc =~ /^unknown_/i) {
        push @unknown_barcodes, $bc;
    } else {
        my $loc_sig = "$best{$bc}{contig}:$best{$bc}{start}:$best{$bc}{stop}";
        $named_locations{$loc_sig} = $bc;
        push @named_barcodes, $bc;
    }
}

# Assign unknown labels: merge with named if same location, else number them
my $unknown_counter = 1;
my @plot_rows;

foreach my $bc (@named_barcodes) {
    push @plot_rows, { label  => $bc,
                       contig => $best{$bc}{contig},
                       start  => $best{$bc}{start},
                       stop   => $best{$bc}{stop},
                       count  => $best{$bc}{count} };
}

foreach my $bc (@unknown_barcodes) {
    my $loc_sig = "$best{$bc}{contig}:$best{$bc}{start}:$best{$bc}{stop}";
    my $label;
    if (exists $named_locations{$loc_sig}) {
        # Same location as a known barcode — use that barcode's label
        $label = $named_locations{$loc_sig} . "_unknown";
    } else {
        $label = "unknown_$unknown_counter";
        $unknown_counter++;
    }
    push @plot_rows, { label  => $label,
                       contig => $best{$bc}{contig},
                       start  => $best{$bc}{start},
                       stop   => $best{$bc}{stop},
                       count  => $best{$bc}{count} };
}

# Write the plot input TSV
# Columns: label  contig  start  stop  count
open(my $tsv, ">", $plot_tsv) or die "Cannot write to $plot_tsv: $!";
print $tsv "label\tcontig\tstart\tstop\tcount\n";
foreach my $row (sort { $a->{label} cmp $b->{label} } @plot_rows) {
    print $tsv "$row->{label}\t$row->{contig}\t$row->{start}\t$row->{stop}\t$row->{count}\n";
}
close($tsv);
print "Plot input TSV written to $plot_tsv.\n";
