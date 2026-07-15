#!/bin/sh
# 02_in_silico_summaries.sh
# 15/07/2026

# 01.- I will use mkCOInr to transform the whole raw database into QIIME format. The taxonomy within this format is easier to use than the COInr native format.
mkdir -p ../../00_Databases/COInrQIIME
format_db.pl \
  -tsv ../../00_Databases/COInr_2025_05_23/COInr.tsv  \
  -taxonomy ../../00_Databases/COInr_2025_05_23/taxonomy.tsv \
  -outdir ../../00_Databases/COInrQIIME/ \
  -out COInrQIIME \
  -outfmt qiime
  
# 02.- Removing the COI sequences (we will only keep the taxa file) to explore taxonomic representation before and after the insilico PCR
rm ../../00_Databases/COInrQIIME/COInrQIIME_trainseq.fasta
