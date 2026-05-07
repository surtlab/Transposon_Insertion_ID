#!/usr/bin/perl -w
use strict;

# Input and output files
my $input_file   = $ARGV[0] or die "Usage: $0 <input.fastq/.fasta> <barcode_file> <output.fasta>\n";
my $barcode_file = $ARGV[1] or die "Usage: $0 <input.fastq/.fasta> <barcode_file> <output.fasta>\n";
my $output_fasta = $ARGV[2] or die "Usage: $0 <input.fastq/.fasta> <barcode_file> <output.fasta>\n";

# Barcode reference file (barcode, forward_seq, reverse_seq)
my %barcode_lookup;

# Load barcode reference table (barcode ID, forward sequence, reverse sequence)
open(my $bfile, "<", $barcode_file) or die "Cannot open $barcode_file: $!";
while (<$bfile>) {
    chomp;
    my ($id, $fwd, $rev) = split(/\t/);
    $barcode_lookup{$fwd} = $id;  # Map forward barcode sequence to ID
    $barcode_lookup{$rev} = $id;  # Map reverse barcode sequence to ID
}
close($bfile);

# Temporary FASTA conversion if needed
my @fasta_lines;

if ($input_file =~ /\.fastq$/i) {
    open(my $fastq, "<", $input_file) or die "Cannot open $input_file: $!";
    while (my $header = <$fastq>) {
        my $seq  = <$fastq>;
        my $plus = <$fastq>;
        my $qual = <$fastq>;
        chomp($header, $seq);
        $header =~ s/^@/>/;  # Convert FASTQ to FASTA header
        push @fasta_lines, $header, $seq;
    }
    close($fastq);
} elsif ($input_file =~ /\.fasta$/i || $input_file =~ /\.fa$/i) {
    open(my $fasta, "<", $input_file) or die "Cannot open $input_file: $!";
    my $header = "";
    my $sequence = "";
    while (my $line = <$fasta>) {
        chomp($line);
        if ($line =~ /^>/) {
            if ($header) {
                push @fasta_lines, $header, $sequence;
            }
            $header = $line;
            $sequence = "";
        } else {
            $sequence .= $line;
        }
    }
    # Push last record
    push @fasta_lines, $header, $sequence if $header;
    close($fasta);
} else {
    die "Unsupported file format. Please provide a .fastq or .fasta file.\n";
}

# Output file (only the second output)
open(my $out, ">", $output_fasta) or die "Cannot write to $output_fasta: $!";

# Process each sequence
for (my $i = 0; $i < @fasta_lines; $i += 2) {
    my $header = $fasta_lines[$i];
    my $seq = $fasta_lines[$i + 1];
    my ($barcode, $barcode_start, $barcode_end);

    # Match forward orientation (look for the forward barcode in the sequence)
    if ($seq =~ /GCAGCGTACG(.*?)AGAGACCTCGTG/) {
        $barcode = $1;
        $barcode_start = index($seq, "GCAGCGTACG$barcode") + length("GCAGCGTACG");
        $barcode_end = $barcode_start + length($barcode);
    }
    # Match reverse orientation (look for the reverse barcode in the sequence)
    elsif ($seq =~ /GTCCACGAGGTCTCT(.*?)CGTACGCTGCAGGTCG/) {
        $barcode = $1;
        $barcode_start = index($seq, "GTCCACGAGGTCTCT$barcode") + length("GTCCACGAGGTCTCT");
        $barcode_end = $barcode_start + length($barcode);
    }

    next unless defined $barcode;

    # Check if barcode matches any in the reference file (forward or reverse)
    my $barcode_id = $barcode_lookup{$barcode} // "unknown";

    # Output - right side (based on a motif match)
    if ($seq =~ /CGTACGCTGCAGGTCG/) {
        my $after_side = substr($seq, $barcode_end + 50);
        print $out ">$barcode_id\_${barcode}_right\n$after_side\n";
    }

    # Output - left side (based on a motif match)
    if ($seq =~ /CCGGCCGTCGACCTGCAGCGTACG/) {
        my $start_pos_left = $barcode_start - 1000;
        $start_pos_left = 0 if $start_pos_left < 0;
        my $before_left = substr($seq, $start_pos_left, 950);
        print $out ">$barcode_id\_${barcode}_left\n$before_left\n";
    }
}

close($out);

print "✅ Output written to: $output_fasta\n";

