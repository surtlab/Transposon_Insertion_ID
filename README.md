# One Shot Transposon Insertion Identification Pipeline

All of the required files to set up the One Shot Transposon Insertion Identification Pipeline are available through this Github site, along with example data that can be used to test the pipeline.

---

## Setup

As a first step, please download all files listed below into a separate directory.

**Scripts to set up environment:**
- `environment.yml`
- `setup_r_packages.R`

**Scripts that are run by the pipeline:**
- `FastQ-reads_to_barcode_FastA.pl`
- `blast_read_output_with_barcode.pl`
- `plot_genomic_region.R`

**Wrapper script for the pipeline:**
- `run_transposon_pipeline.sh`

---

## Environment Setup

After downloading, use the `environment.yml` file to set up a conda environment to carry out all downstream analyses. Open terminal and run:

```bash
conda env create -f environment.yml
conda activate transposon_pipeline
```

Once this environment is activated, download and set up all required R packages:

```bash
Rscript setup_r_packages.R
```

---

## Input Files

Once all environments are set up and activated, download a GenBank flat file (`.gbk` or `.gbff`) for the genome you would like to query for transposon insertions. This can be either a draft genome or a complete genome, but make sure nucleotide sequences are included in the `.gbk`/`.gbff` file. Also download Nanopore read files in FastQ format to submit to the pipeline.

---

## Running the Pipeline

Once these pieces are in place, all that is needed to run the pipeline is to run the wrapper script in terminal:

```bash
bash run_transposon_pipeline.sh <input.fastq> <genome.gbff>
```

---

## What the Pipeline Does

**Step 1 — Barcode extraction**
The first step screens FastQ reads for hard-coded sequences that bracket the barcode locations in the Tn5 transposons. Reads that include both sequence brackets to the barcode will be pulled into an output `.reads` file where the entire read will be in FastA format. These reads will be named with the barcode, if barcode sequences are matched from the `barcodes.txt` file. If there is no barcode matched (either because the barcode is novel or because of a sequencing error), the read will be labelled `unknown`.

**Step 2 — BLAST database setup and search**
BlastN is then used to query reads containing the barcodes against the genome of interest. To set up BLAST databases, the `.gbk`/`.gbff` file is converted to FastA nucleotide sequence, and `makeblastdb` is used to create the BLAST database.

**Step 3 — Insertion site identification**
Results from BlastN are output to a `.blast` file in `-outfmt 6`. Regions identified by this BLAST search will be the insertion sites of the transposon in the strains/barcodes of interest.

**Step 4 — Genomic context visualization**
BlastN results highlighting the transposon insertion site are input into an R script that uses the `.gbk`/`.gbff` file to create a figure of the insertion region in the genome.
