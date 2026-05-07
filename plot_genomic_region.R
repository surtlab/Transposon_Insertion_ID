#!/usr/bin/env Rscript
# ============================================================
# plot_genomic_region.R
#
# Visualize a genomic region from a GenBank/GBFF file, showing
# the highlighted region plus N upstream and downstream genes.
# Supports multi-record .gbff files via --contig selection.
#
# NO Bioconductor dependencies — pure base R parser.
#
# Dependencies (install once):
#   install.packages(c("optparse", "ggplot2", "gggenes", "dplyr"))
#
# Usage:
#   # List all contigs in a file:
#   Rscript plot_genomic_region.R --gbk genome.gbff --list-contigs
#
#   # Plot a region on a specific contig:
#   Rscript plot_genomic_region.R \
#     --gbk    genome.gbff \
#     --contig NZ_CP012345.1 \
#     --start  45000 \
#     --end    48000 \
#     --out    region_map.pdf \
#     --label  "BLAST hit"
#
#   # Single-record file — --contig can be omitted:
#   Rscript plot_genomic_region.R \
#     --gbk single.gbk --start 45000 --end 48000
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(gggenes)
  library(dplyr)
})

# ── CLI arguments ────────────────────────────────────────────
option_list <- list(
  make_option(c("-g", "--gbk"),
              type = "character",
              help = "Path to GenBank / GBFF file [required]"),

  make_option(c("-c", "--contig"),
              type = "character", default = NULL,
              help = paste("Contig / accession to use (e.g. NZ_CP012345.1).",
                           "Omit only if the file has a single record.",
                           "Use --list-contigs to see all available IDs.")),

  make_option("--list-contigs",
              action = "store_true", default = FALSE,
              help = "Print all contig IDs in the file and exit"),

  make_option(c("-s", "--start"),
              type = "integer", default = NULL,
              help = "Region start coordinate in bp [required for plotting]"),

  make_option(c("-e", "--end"),
              type = "integer", default = NULL,
              help = "Region end coordinate in bp   [required for plotting]"),

  make_option(c("-o", "--out"),
              type = "character", default = "genomic_region.pdf",
              help = "Output file (.pdf / .png / .svg) [default: genomic_region.pdf]"),

  make_option(c("-l", "--label"),
              type = "character", default = "Region of Interest",
              help = "Label printed above the highlighted region"),

  make_option(c("-n", "--flanking"),
              type = "integer", default = 2,
              help = "Number of flanking genes on each side [default: 2]"),

  make_option("--width",
              type = "double", default = 12,
              help = "Plot width  in inches [default: 12]"),

  make_option("--height",
              type = "double", default = 5,
              help = "Plot height in inches [default: 5]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$gbk))      stop("--gbk is required. Run with --help for usage.")
if (!file.exists(opt$gbk)) stop("File not found: ", opt$gbk)

# ╔══════════════════════════════════════════════════════════════╗
# ║  Pure base-R GenBank / GBFF parser — no Bioconductor needed ║
# ╚══════════════════════════════════════════════════════════════╝

# ── Split raw lines into per-record blocks on '//' ────────────
split_records <- function(lines) {
  ends   <- which(trimws(lines) == "//")
  if (length(ends) == 0) stop("No '//' record terminators found.")
  starts <- c(1L, head(ends + 1L, -1L))
  lapply(seq_along(ends), function(i) lines[starts[i]:ends[i]])
}

# ── Pull the second whitespace-delimited token from a tagged line
first_match_token <- function(lines, pattern, token = 2L) {
  hit <- lines[grepl(pattern, lines)]
  if (length(hit) == 0) return(NA_character_)
  toks <- strsplit(trimws(hit[1]), "\\s+")[[1]]
  if (length(toks) < token) return(NA_character_)
  toks[token]
}

# ── Build record metadata from raw lines ─────────────────────
record_meta <- function(lines, idx) {
  locus <- first_match_token(lines, "^LOCUS")
  acc   <- first_match_token(lines, "^ACCESSION")
  ver   <- first_match_token(lines, "^VERSION")
  id    <- if (!is.na(ver))   ver   else
           if (!is.na(acc))   acc   else
           if (!is.na(locus)) locus else paste0("record_", idx)
  list(locus = locus, accession = acc, version = ver, id = id, lines = lines)
}

# ── Parse the FEATURES table from one record's lines ─────────
#
# GenBank feature table format:
#   cols 1-5   : blank
#   cols 6-20  : feature key (e.g. "gene", "CDS")
#   cols 22+   : location or qualifier (/key="value")
#
# We collect gene, CDS, and misc_feature entries and extract
# /gene=, /locus_tag=, /product= qualifiers.
parse_features <- function(lines) {

  # Find FEATURES block boundaries
  feat_start <- which(grepl("^FEATURES", lines))
  if (length(feat_start) == 0) return(data.frame())
  # Block ends at ORIGIN or end of lines
  origin <- which(grepl("^ORIGIN", lines))
  feat_end <- if (length(origin) > 0) origin[1] - 1L else length(lines)
  feat_lines <- lines[feat_start:feat_end]

  # Feature key lines: non-blank in cols 6-20, blank in cols 1-5
  # In a standard GenBank flat file, feature keys start at column 6 (1-indexed)
  # and qualifiers start at column 22 with a leading '/'
  is_feature_key <- function(x) {
    grepl("^     \\S", x) & !grepl("^     /", x)
  }
  is_qualifier <- function(x) {
    grepl("^                     /", x) | grepl("^                      /", x)
  }
  is_continuation <- function(x) {
    grepl("^                     [^/]", x) | grepl("^                      [^/]", x)
  }

  # Target feature types
  target_keys <- c("gene", "CDS", "misc_feature", "mobile_element",
                   "repeat_region", "ncRNA", "rRNA", "tRNA")

  records <- list()
  current_key      <- NULL
  current_location <- NULL
  current_quals    <- list()

  flush_feature <- function() {
    if (is.null(current_key)) return()
    if (!current_key %in% target_keys) return()

    loc <- parse_location(current_location)
    if (is.null(loc)) return()

    # Capture locus_tag and gene symbol separately so we can show both
    locus_tag   <- current_quals[["locus_tag"]]
    gene_symbol <- current_quals[["gene"]]
    product     <- current_quals[["product"]]

    if (is.null(locus_tag)   || locus_tag   == "") locus_tag   <- NA_character_
    if (is.null(gene_symbol) || gene_symbol == "") gene_symbol <- NA_character_
    if (is.null(product)     || product     == "") product     <- NA_character_

    # Fallback display name if locus_tag is absent
    fallback <- if (!is.na(locus_tag)) locus_tag else
                if (!is.na(product))   product   else current_key

    records[[length(records) + 1L]] <<- data.frame(
      start       = loc$start,
      end         = loc$end,
      strand_dir  = loc$strand,
      locus_tag   = ifelse(is.na(locus_tag),   fallback, locus_tag),
      gene_symbol = ifelse(is.na(gene_symbol), NA_character_, gene_symbol),
      stringsAsFactors = FALSE
    )
  }

  # Parse location string → start, end, strand
  # Handles: 100..200  complement(100..200)  join(...)  <100..>200
  parse_location <- function(loc) {
    if (is.null(loc) || loc == "") return(NULL)
    strand <- 1L
    if (grepl("complement", loc)) strand <- -1L
    # Extract all numbers
    nums <- as.integer(regmatches(loc, gregexpr("[0-9]+", loc))[[1]])
    if (length(nums) < 2) {
      if (length(nums) == 1) return(list(start = nums[1], end = nums[1], strand = strand))
      return(NULL)
    }
    list(start = min(nums), end = max(nums), strand = strand)
  }

  active_qual_key  <- NULL
  active_qual_val  <- NULL

  flush_qualifier <- function() {
    if (!is.null(active_qual_key))
      current_quals[[active_qual_key]] <<- active_qual_val
    active_qual_key <<- NULL
    active_qual_val <<- NULL
  }

  for (ln in feat_lines) {

    if (is_feature_key(ln)) {
      # Save previous qualifier and feature
      flush_qualifier()
      flush_feature()
      current_key      <- trimws(substr(ln, 6, 20))
      current_location <- trimws(substr(ln, 22, nchar(ln)))
      current_quals    <- list()
      active_qual_key  <- NULL
      active_qual_val  <- NULL

    } else if (grepl("^                     /|^                      /", ln)) {
      flush_qualifier()
      qual_str <- trimws(sub("^\\s+/", "", ln))
      eq_pos   <- regexpr("=", qual_str, fixed = TRUE)
      if (eq_pos > 0) {
        active_qual_key <- substr(qual_str, 1, eq_pos - 1)
        val             <- substr(qual_str, eq_pos + 1, nchar(qual_str))
        # Strip surrounding quotes
        active_qual_val <- gsub('^"|"$', "", val)
      } else {
        # Flag qualifier with no value (e.g. /pseudo)
        active_qual_key <- qual_str
        active_qual_val <- "true"
      }

    } else if (!is.null(active_qual_key) &&
               grepl("^                     |^                      ", ln) &&
               !grepl("^ORIGIN|^//", ln)) {
      # Continuation line for a multi-line qualifier value
      continuation <- trimws(ln)
      # Remove closing quote if present and re-append
      active_qual_val <- paste0(gsub('"$', "", active_qual_val),
                                gsub('^"|"$', "", continuation))

    } else if (grepl("^                     ", ln) && !grepl("^\\s*/", ln)) {
      # Location continuation (wrapped location string)
      current_location <- paste0(current_location, trimws(ln))
    }
  }

  # Flush the last feature
  flush_qualifier()
  flush_feature()

  if (length(records) == 0) return(data.frame())
  do.call(rbind, records)
}

# ══════════════════════════════════════════════════════════════

# ── Read and split file ───────────────────────────────────────
cat("Reading:", opt$gbk, "\n")
raw_lines   <- readLines(opt$gbk, warn = FALSE)
record_list <- split_records(raw_lines)
all_records <- lapply(seq_along(record_list),
                      function(i) record_meta(record_list[[i]], i))
cat("Records found:", length(all_records), "\n")

# ── --list-contigs mode ───────────────────────────────────────
if (isTRUE(opt[["list-contigs"]])) {
  cat("\nAvailable contigs (pass the ID column to --contig):\n\n")
  cat(sprintf("  %-4s  %-26s  %-26s  %s\n", "#", "LOCUS", "ACCESSION", "VERSION / ID"))
  cat("  ", strrep("-", 72), "\n", sep = "")
  for (i in seq_along(all_records)) {
    r <- all_records[[i]]
    cat(sprintf("  %-4d  %-26s  %-26s  %s\n", i,
                ifelse(is.na(r$locus),     "-", r$locus),
                ifelse(is.na(r$accession), "-", r$accession),
                r$id))
  }
  cat("\n")
  quit(status = 0)
}

# ── Require coordinates ───────────────────────────────────────
if (is.null(opt$start) || is.null(opt$end))
  stop("--start and --end are required for plotting.\n",
       "Use --list-contigs to inspect the file first.")

region_start <- opt$start
region_end   <- opt$end
n_flank      <- opt$flanking

# ── Select target contig ──────────────────────────────────────
if (!is.null(opt$contig)) {
  query <- tolower(opt$contig)
  idx   <- which(sapply(all_records, function(r)
    tolower(r$locus)     == query |
    tolower(r$accession) == query |
    tolower(r$version)   == query |
    tolower(r$id)        == query
  ))
  if (length(idx) == 0) {
    cat("\nERROR: Contig '", opt$contig, "' not found.\n\n", sep = "")
    cat("Available IDs:\n")
    for (r in all_records)
      cat("  ", r$id, "\n", sep = "")
    stop("Use --list-contigs for the full table.")
  }
  chosen <- all_records[[idx[1]]]
} else {
  if (length(all_records) > 1)
    stop("File has ", length(all_records), " records. ",
         "Specify one with --contig <ID>. ",
         "Run --list-contigs to see all IDs.")
  chosen <- all_records[[1]]
}

cat("Using contig:", chosen$id,
    "(locus:", ifelse(is.na(chosen$locus), "?", chosen$locus), ")\n")

# ── Parse features ────────────────────────────────────────────
cat("Parsing annotations...\n")
feat_df <- parse_features(chosen$lines)

if (nrow(feat_df) == 0)
  stop("No gene/CDS features found in contig ", chosen$id,
       ". Check that the file has a populated FEATURES table.")

cat("Annotated features (raw):", nrow(feat_df), "\n")

# NCBI gbff files have both a 'gene' and 'CDS' feature for the same locus.
# Keep only one row per unique start/end/strand combination.
feat_df <- feat_df %>%
  arrange(start, end) %>%
  distinct(start, end, strand_dir, .keep_all = TRUE)

cat("Annotated features (deduplicated):", nrow(feat_df), "\n")

# ── Find flanking genes ───────────────────────────────────────
# First identify genes overlapping the insertion region
overlapping <- feat_df %>%
  filter(start < region_end & end > region_start)

# Determine the strand of the hit gene(s).
# If the insertion overlaps multiple genes on different strands,
# the majority strand wins; ties default to positive.
if (nrow(overlapping) > 0) {
  strand_votes <- sum(overlapping$strand_dir)
  hit_strand   <- if (strand_votes >= 0) 1L else -1L
} else {
  hit_strand <- 1L  # no overlapping gene — default to positive
}

cat("Hit gene strand:", if (hit_strand == 1) "positive (+)" else "negative (-)", "\n")

# Genes left of the region (lower coordinates)
left_genes <- feat_df %>%
  filter(end <= region_start) %>%
  arrange(desc(end))          # closest first

# Genes right of the region (higher coordinates)
right_genes <- feat_df %>%
  filter(start >= region_end) %>%
  arrange(start)              # closest first

# Assign upstream/downstream based on strand:
#   Positive strand: upstream = lower coords, downstream = higher coords
#   Negative strand: upstream = higher coords, downstream = lower coords
if (hit_strand == 1L) {
  upstream   <- slice_head(left_genes,  n = n_flank)
  downstream <- slice_head(right_genes, n = n_flank)
} else {
  upstream   <- slice_head(right_genes, n = n_flank)
  downstream <- slice_head(left_genes,  n = n_flank)
}

cat("Upstream   :", nrow(upstream),
    paste(upstream$locus_tag, collapse = ", "), "\n")
cat("In region  :", nrow(overlapping),
    paste(overlapping$locus_tag, collapse = ", "), "\n")
cat("Downstream :", nrow(downstream),
    paste(downstream$locus_tag, collapse = ", "), "\n")

# Pipeline summary block
cat("\n## PIPELINE_INFO_START\n")
cat("contig=",           chosen$id,                                    "\n", sep = "")
cat("region_start=",     region_start,                                 "\n", sep = "")
cat("region_end=",       region_end,                                   "\n", sep = "")
cat("upstream_genes=",   paste(upstream$locus_tag,    collapse = ","), "\n", sep = "")
cat("hit_genes=",        paste(overlapping$locus_tag, collapse = ","), "\n", sep = "")
cat("downstream_genes=", paste(downstream$locus_tag,  collapse = ","), "\n", sep = "")
cat("output_file=",      opt$out,                                      "\n", sep = "")
cat("## PIPELINE_INFO_END\n\n")

# ── Build plot data ───────────────────────────────────────────
plot_genes <- bind_rows(
  upstream    %>% mutate(category = "upstream"),
  overlapping %>% mutate(category = "highlight"),
  downstream  %>% mutate(category = "downstream")
) %>%
  distinct(start, end, .keep_all = TRUE) %>%
  arrange(start) %>%
  mutate(
    molecule   = factor("genome"),
    # Two-line label: gene symbol (if present) on top, locus_tag below
    label_text  = ifelse(
      !is.na(gene_symbol) & gene_symbol != "",
      paste0(gene_symbol, "\n", locus_tag),
      locus_tag
    ),
    has_symbol  = !is.na(gene_symbol) & gene_symbol != "",
    # nudge locus_tag down when a symbol is also shown, else centre it
    lt_nudge    = ifelse(!is.na(gene_symbol) & gene_symbol != "", -0.12, 0)
  )

if (nrow(plot_genes) == 0)
  stop("No genes to plot. Check --start / --end are valid for contig ", chosen$id)

plot_left  <- min(plot_genes$start) - 300
plot_right <- max(plot_genes$end)   + 300

# ── Palette ───────────────────────────────────────────────────
pal <- c(
  upstream   = "#5B9BD5",
  highlight  = "#E84855",
  downstream = "#70B77E"
)

# ── Plot ──────────────────────────────────────────────────────
region_mid  <- (region_start + region_end) / 2
bracket_gap <- (region_end - region_start) * 0.06

p <- ggplot(plot_genes,
            aes(xmin    = start,
                xmax    = end,
                y       = molecule,
                fill    = category,
                forward = strand_dir == 1,
                label   = label_text)) +

  annotate("rect",
           xmin = region_start, xmax = region_end,
           ymin = -Inf, ymax = Inf,
           fill = "#FFF3CD", alpha = 0.6) +

  annotate("segment",
           x = region_start, xend = region_start, y = 1.1, yend = 2.1,
           colour = "#C9A227", linewidth = 0.8, linetype = "dashed") +
  annotate("segment",
           x = region_end, xend = region_end, y = 1.1, yend = 2.1,
           colour = "#C9A227", linewidth = 0.8, linetype = "dashed") +

  annotate("segment",
           x = region_start, xend = region_mid - bracket_gap,
           y = 2.05, yend = 2.05,
           colour = "#C9A227", linewidth = 0.5) +
  annotate("segment",
           x = region_mid + bracket_gap, xend = region_end,
           y = 2.05, yend = 2.05,
           colour = "#C9A227", linewidth = 0.5) +

  annotate("text",
           x = region_mid, y = 2.18,
           label  = opt$label,
           colour = "#9B6A00", size = 3.8, fontface = "bold", hjust = 0.5) +

  geom_gene_arrow(
    arrowhead_height  = unit(14, "mm"),
    arrowhead_width   = unit(5,  "mm"),
    arrow_body_height = unit(12, "mm"),
    colour    = "white",
    linewidth = 0.4
  ) +

  # Two-layer labelling inside arrows:
  # gene symbol (bold, nudged up) when present + locus_tag (smaller, nudged down)
  geom_text(
    data        = function(d) d[d$has_symbol, ],
    aes(x = (start + end) / 2, y = molecule, label = gene_symbol),
    size        = 3,
    fontface    = "bold",
    colour      = "white",
    nudge_y     = 0.13,
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(x = (start + end) / 2, y = molecule, label = locus_tag),
    size        = 2.3,
    fontface    = "plain",
    colour      = "white",
    nudge_y     = plot_genes$lt_nudge,
    inherit.aes = FALSE
  ) +

  geom_hline(yintercept = 1.5, colour = "grey55", linewidth = 0.5) +

  scale_fill_manual(
    values = pal,
    labels = c(upstream   = "Upstream (5' of insertion)",
               highlight  = "Gene at insertion",
               downstream = "Downstream (3' of insertion)"),
    name   = "Gene context"
  ) +
  scale_x_continuous(
    labels = function(x) paste0(formatC(x / 1000, format = "f", digits = 1), " kb"),
    expand = expansion(mult = 0.03)
  ) +
  scale_y_discrete(expand = expansion(add = 1.4)) +

  theme_classic() +
  theme(
    axis.title.y       = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.line.y        = element_blank(),
    axis.title.x       = element_text(size = 10, colour = "grey35"),
    axis.text.x        = element_text(size = 9,  colour = "grey35"),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9, face = "bold"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "grey50"),
    plot.margin        = margin(12, 16, 10, 16)
  ) +
  labs(
    title    = paste0("Genomic Context: ",
                      formatC(region_start, big.mark = ",", format = "d"),
                      " - ",
                      formatC(region_end,   big.mark = ",", format = "d"), " bp"),
    subtitle = paste0("Contig: ", chosen$id,
                      " | ", n_flank, " flanking genes each side",
                      " | ", basename(opt$gbk)),
    x = "Genomic position"
  )

# ── Save ──────────────────────────────────────────────────────
ext <- tolower(tools::file_ext(opt$out))
cat("Saving ->", opt$out, "\n")

if (ext == "png") {
  ggsave(opt$out, p, width = opt$width, height = opt$height, dpi = 300, bg = "white")
} else if (ext == "svg") {
  ggsave(opt$out, p, width = opt$width, height = opt$height, device = "svg", bg = "white")
} else {
  ggsave(opt$out, p, width = opt$width, height = opt$height, device = "pdf")
}

cat("Done!\n")
