#!/bin/sh
# 01_mkCOInr.sh
# 12/06/2025

# 00.- Prerequisites:
# We will create our COI reference database using Meglecz 2023 COInr database and mkCOInr software.
# https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13756

# We have downloaded version 2025-05-23 of COInr
# https://zenodo.org/records/15515860
# We have downloaded saved and uncompressed this datafile at 'PUF-PAS_AirDNA/00_Databases/COInr_2025_05_23/'

# Finally make sure to have installed mkCOInr third-party auxiliary programs using the instructions found here:
# https://mkcoinr.readthedocs.io/en/latest/content/installation.html
# Since we followed those instuctions be sure your conda environment to be set at:
# conda activate mkCOInr

# I have added mkCOInr to my PATH by concatinating to my .bashrc:
# export PATH="/home/luis/Software/mkCOInr/scripts:$PATH"

# I made all mkCOInr perl scripts executable with:
# chmod +x /home/luis/Software/mkCOInr/scripts/*.pl

# And prepended a perl shebang to all of them using:
# for pl in ./*.pl; do     sed -i '1i\#!/usr/bin/env perl' "$pl"; done

# 01.- insilico PCR for our Leray primers (ca. 313 bp)
# https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34
# -e_pcr: 1 for insilico PCR algorithm
# -trim_error: Mismatches maximum proportion between the primers and targeted sequences
# -min_overlap: minimum overlap between primers and sequences during PCR
select_region.pl \
  -tsv ../../00_Databases/COInr_2025_05_23/COInr.tsv \
  -outdir ../../00_Databases/COInrLeray \
  -e_pcr 1 \
  -fw GGWACWGGWTGAACWGTWTAYCCYCC \
  -rv TAIACYTCIGGRTGICCRAARAAYCA \
  -trim_error 0.3 \
  -min_amplicon_length 280 \
  -max_amplicon_length 345 \
  -min_overlap 20

# We recovered 5748708 unique sequences from the 5933866 originals in the 2025 COInr database
# 02.- Transforming our database into qiime
format_db.pl \
  -tsv ../../00_Databases/COInrLeray/trimmed.tsv \
  -taxonomy ../../00_Databases/COInr_2025_05_23/taxonomy.tsv \
  -outdir ../../00_Databases/COInrLerayQIIME/ \
  -out COInrLerayQIIME \
  -outfmt qiime
