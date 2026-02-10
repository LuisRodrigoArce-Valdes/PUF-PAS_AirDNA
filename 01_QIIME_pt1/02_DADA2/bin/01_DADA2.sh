#!/bin/sh
# 01_DADA2.sh
# 21/10/2025

# 00.- This second analyses was performed after Deblur's one. Thus, we need some files that were created during that script.
# We will also work this in the qiime2 environment.

# 01.- Symbolic-linking input files.
# DADA2 works with *unpaired*, *non-quality filtered*, raw sequences without primers. 
# These are the output files of step 04 of the DEBLUR pipeline (Dragonflies_Metabarcoding/02_Deblur/data/02_cutadapt.qza).
# We will simbolic link this file to use it as input without copying it.
ln -s -f ../../01_Deblur/data/02_cutadapt.qza ../data/01_cutadapt.qza

# 02.- Denoising with DADA2
# We won't trim nor trunc at a specific lenght (trims and truncs with 0s). Instead, we will trim when qualities become below a treshold as we did it on Deblur.
# Aditionally, we will test several tresholds to evaluate parameter sensibility.
# As with Deblur we will test several minimum overalps for forward and reverse reads to be paired.
# Pairing will only be done when the expected number of errors both in forward and the revearse read is less than 1.
# All other parameters are the defaults.
mkdir -p ./logs
for t in 0 20 26 30
do
	for l in 10 50 100 200
	do
		echo "Denoising using t = ${t}; l = ${l}"
		qiime dada2 denoise-paired \
		  --i-demultiplexed-seqs ../data/01_cutadapt.qza \
		  --p-trim-left-f 0 \
		  --p-trim-left-r 0 \
		  --p-trunc-len-f 0 \
		  --p-trunc-len-r 0 \
		  --p-trunc-q ${t} \
		  --p-min-overlap ${l} \
		  --p-max-ee-f 1 \
		  --p-max-ee-r 1 \
		  --p-pooling-method independent \
		  --p-chimera-method consensus \
		  --p-n-reads-learn 1000000 \
		  --o-representative-sequences ../data/02_denoised_seqs_t${t}_l${l}.qza \
		  --o-table ../data/02_denoised_table_t${t}_l${l}.qza \
		  --o-denoising-stats ../data/02_denoised_stats_t${t}_l${l}.qza \
		  --p-n-threads 4 \
		  --verbose 2>&1 | tee ./logs/01_DADA2_t${t}_l${l}.log
		echo "Finishing denoising t = ${t}; l = ${l}"
	done
done

### *NOTE THAT THE MOST ASTRINGENT FILTERING SCHEMES COULDN'T BE COMPLETED! MOST PROBABLY DADA2 COULDN'T RECOVER READS WHEN BEING TO STRINGENT!

# Visualizing DADA2's results
for t in 0 20 26 30
do
	for l in 10 50 100 200
	do
		echo "Visualizing t = ${t}; l = ${l}"
		qiime metadata tabulate \
		  --m-input-file ../data/02_denoised_stats_t${t}_l${l}.qza \
		  --o-visualization ../results/01_denoised_stats_t${t}_l${l}.qzv \
		  --verbose
	
		  qiime feature-table summarize \
		    --i-table ../data/02_denoised_table_t${t}_l${l}.qza \
		    --o-visualization ../results/01_denoised_table_t${t}_l${l}.qzv \
		    --verbose
	    
		  qiime feature-table tabulate-seqs \
		    --i-data ../data/02_denoised_seqs_t${t}_l${l}.qza \
		    --o-visualization ../results/01_denoised_seqs_t${t}_l${l}.qzv \
		    --verbose
		echo echo "Visualizing t = ${t}; l = ${l}"
	done
done

echo "To visualize any result use:"
ls ../results/01_* | sed 's/^/qiime tools view /'

# We will keep the same filtering schemes as with DEBLUR to do a comparison one-by-one between the two algorithms.

# *SEEMS THAT DADA2 IS WORKING BETTTER THAN DEBLUR IN THIS CASE, IT IS CLEANING THE CONTROLS AND RETAINING A LOT OF READS PER SAMPLE!* :)

# Deleting data files of filtering schemes we are not using. Visualisations will be kept in the results folder for future consultation.
rm -f ../data/02_denoised_*_t0_l100.qza \
      ../data/02_denoised_*_t0_l200.qza \
      ../data/02_denoised_*_t0_l50.qza \
      ../data/02_denoised_*_t20_l100.qza \
      ../data/02_denoised_*_t20_l200.qza \
      ../data/02_denoised_*_t26_l100.qza \
      ../data/02_denoised_*_t26_l200.qza \
      ../data/02_denoised_*_t26_l50.qza \
      ../data/02_denoised_*_t30_l100.qza \
      ../data/02_denoised_*_t30_l10.qza \
      ../data/02_denoised_*_t30_l200.qza \
      ../data/02_denoised_*_t30_l50.qza

