#!/bin/sh
# 01_Deblur.sh
# 21/10/2025

# 00.- Prerequisites:
# We will work our bioinformatic analyses in the Qiime2 environment. To follow this scripts install that software locally to your own PC. We installed ours using conda:
# conda env create -n qiime2 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
# conda activate qiime2

# Our first quality control analysis will be done using FastQC. Be sure to also install it!
# conda install bioconda::fastqc

# Raw data files are placed in the 'PUF-PAS_AirDNA/00_RawData' directory.

# We downloaded the raw files from this library directly from the Illumina server, thus, they don't need any editions in their names.

# Creating log files directory
mkdir -p ./logs

# 01.- Evaluating raw data files quality using fastQC
echo "Copying files"
mkdir -p ../data/tmp
mkdir -p ../results/00_FastQC
cp ../../../00_RawData/PUFs_Pt1/*.gz ../data/tmp

for i in $(ls ../data/tmp)
	do
	echo "Working on sample $i"
	echo ""
	fastqc -o ../results/00_FastQC/ ../data/tmp/$i
done

#Removing zip files; we will keep only the htmls
echo "Removing .zips"
rm ../results/00_FastQC/*.zip
echo ""

# 02.- Importing files into qiime2.
# When importing into qiime double check the format of your sequences: https://docs.qiime2.org/2024.10/tutorials/importing/
echo "Importing into qiime"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../data/tmp \
  --output-path ../data/01_raw.qza \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --validate-level max
rm -r ../data/tmp
echo "Finishing import"

# 03.- Visualizing raw data files; --p-n: number of random sequences for which to analyse quality
echo "Beginning visualisation"
qiime demux summarize \
  --i-data ../data/01_raw.qza \
  --p-n 500000 \
  --o-visualization ../results/01_raw.qzv
echo "Finishing visualisation. To view results use:"
echo "qiime tools view ../results/01_raw.qzv"

# 04.- Trimming PCR primers (since we used degenerate primers, we can't trust the reads at those positions due to mismatch amplification):
# https://docs.qiime2.org/2024.10/plugins/available/cutadapt/trim-paired/
# https://frontiersinzoology.biomedcentral.com/articles/10.1186/1742-9994-10-34
# --p-front-f: Forward primer, found at 5' of the forward read.
# --p-front-r: Reverse primer, found at 5' of the reverse read.
# --o-trimmed-sequences: Discard reads in which no adapter was found.
echo "Filtering with cutadapt"
qiime cutadapt trim-paired \
  --p-cores 4 \
  --i-demultiplexed-sequences ../data/01_raw.qza \
  --p-front-f ^GGWACWGGWTGAACWGTWTAYCCYCC \
  --p-front-r ^TAIACYTCIGGRTGICCRAARAAYCA \
  --p-discard-untrimmed \
  --o-trimmed-sequences ../data/02_cutadapt.qza \
  --verbose 2>&1 | tee ./logs/01_cutadapt.log
echo "Finished filtering"

echo "Visualizing filtered dataset"
qiime demux summarize \
  --i-data ../data/02_cutadapt.qza \
  --p-n 500000 \
  --o-visualization ../results/02_cutadapt.qzv
echo "Finishing visualisation. To view results use:"
echo "qiime tools view ../results/02_cutadapt.qzv"

# 05.- Pair merging & sequence filtering.
# We amplified using the Leray's primers emplyed by Morril et al. 2021 and designed by Leray et al. 2013.
# They amplify a 313 bp COI fragment. Then, we sequenced using Illumina's v3 2x300. 
# Thus, there should be plenty of overlap between both R1 and R2 reads to merge them together. Following this first scheme we will merge them to increase sequence quality. 
# However, instad of using PEAR, we will use qiime vsearch!
# https://forum.qiime2.org/t/plugin-like-pear/23299
# Sequencing structure stats:
# + Expected fragment length: 313 pb
# + Primer length: 26 pb (both forward and reverse)
# + Total amplicon lenght: 313 + 26 + 26 = 365 pb
# + Read on each direction after primer removal: 300 - 26 = 274 pb
# + Lenght not overlaped that each sequencing direction doesn't reach: 313 - 274 = 39 pb
# + Expected overlap lenght: 365 - 26 - 26 - 39 - 39 = 235 pb

# The numbers in qiimes reports vary a little bit from these ones because plenty of samples were sequenced with a total lenght of 301 pb instead of the 2x300.
# This last nucleotide is of very low quality. See:
# qiime tools view ../results/01_raw.qzv

# Running qiime vsearch inside a double for loop to evaluate parameter sensibility
# --p-truncqual : threshold to truncate sequences after quality drops to or below it (e.g. 26 = 99.8% base call accuracy) (0 = no filter).
# --p-minovlen : minimum overlap after truncation to merge the two reads (default = 10).
# --p-maxee 1 : With this we will keep only merged fragments with a number of expected errors (sum of error pob per nucleotide) lower than 1. Thus, it should minimize erroneous reads.

# *Note that R2 reads lose quality significantly faster than R1, i. e., at smaller cycle numbers. This is common. But in our case is particularly evident.
# *Ideally, we want to keep as much sequences with the highest quality as possible.

# Preliminary tests showed sequences with a length much higher than the maximum that we would expect: 313 nt.
# These longer fragments had lower qualities at the extra nucleotides. Thus, they are probably miss alignments. We can filter them out!
# Additonally deblurs algorithm works only when sequences have the same lenghts! We will keep only sequences of 313.

# --p-minmergelen 313 : To also remove sequences shorter than this lenght.
# --p-maxmergelen 313 : We must not retain sequences with lengths higher than this.

# Let's try some parameters combinations!
for t in 0 20 26 30
do
	for l in 10 50 100 200
	do
		echo "Merging pairs t = ${t}; l = ${l}"
		qiime vsearch merge-pairs \
		  --i-demultiplexed-seqs ../data/02_cutadapt.qza \
		  --o-merged-sequences ../data/03_merged_t${t}_l${l}.qza \
		  --o-unmerged-sequences ../data/03_unmerged_t${t}_l${l}.qza \
		  --p-truncqual ${t} \
		  --p-minovlen ${l} \
		  --p-maxee 1 \
		  --p-minmergelen 313 \
		  --p-maxmergelen 313 \
		  --p-threads 4 \
		  --verbose 2>&1 | tee ./logs/02_vsearch_t${t}_l${l}.log
		echo "Finished merging t = ${t}; l = ${l}"

		echo "Visualizing merged dataset t = ${t}; l = ${l}"
		qiime demux summarize \
		  --i-data ../data/03_merged_t${t}_l${l}.qza \
		  --p-n 500000 \
		  --o-visualization ../results/03_merged_t${t}_l${l}.qzv
		echo "Finishing visualisation t = ${t}; l = ${l}."
	done
done

echo "Finished 16 parameter combinations!"
echo "To visualize any result use:"
ls ../results/*.qzv | sed 's/^/qiime tools view /'

# Qualities improved just by merging the two directions. However, we can't be very stringent during our filtering. Otherwise we will loose a lot of samples.
# Some of them they have just hundreds or a few thousand reads. I selected some not so stringent filters, trying to keep at least more than 100 reads in the samples (excluding the controls). 

# We will keep: 
# i) The most relaxed filtering (t0; l10)
# ii) a moderate filtering (t20; l10) 
# iii) a stricter filtering based on sequence overlap (t20; l50).
# iv) a strincter filtering based on quality trimming (t26; l10).

# Note that we are trimming quality scores down to t20 (99% certainty of assignment)!
# We will analyze these four datasets to evaluate the consistency of our results.

# Removing all other filtering schemes files
rm -f ../data/03_merged_t0_l100.qza \
      ../data/03_merged_t0_l200.qza \
      ../data/03_merged_t0_l50.qza \
      ../data/03_merged_t20_l100.qza \
      ../data/03_merged_t20_l200.qza \
      ../data/03_merged_t26_l100.qza \
      ../data/03_merged_t26_l200.qza \
      ../data/03_merged_t26_l50.qza \
      ../data/03_merged_t30_l100.qza \
      ../data/03_merged_t30_l10.qza \
      ../data/03_merged_t30_l200.qza \
      ../data/03_merged_t30_l50.qza \
      ../data/03_unmerged_t0_l100.qza \
      ../data/03_unmerged_t0_l200.qza \
      ../data/03_unmerged_t0_l50.qza \
      ../data/03_unmerged_t20_l100.qza \
      ../data/03_unmerged_t20_l200.qza \
      ../data/03_unmerged_t26_l100.qza \
      ../data/03_unmerged_t26_l200.qza \
      ../data/03_unmerged_t26_l50.qza \
      ../data/03_unmerged_t30_l100.qza \
      ../data/03_unmerged_t30_l10.qza \
      ../data/03_unmerged_t30_l200.qza \
      ../data/03_unmerged_t30_l50.qza
      
# 06.- Denoising and deprelicating with deblur
# Now we will run our denoising algorithm using deblur
# -p-trim-length -1 : We have quality filtered our database and we are happy with our filtering. All seqs are 313 nt long. The -1 value disables trimming on the 3'.
# --p-left-trim-len 0 : Again we are disabling trimming due to certainty in our quality filtering. All seqs are 313 nt long. Now we are disabling trimming on the 5'.
# --p-indel-max 3 : Maximum number of allowed indels of 3. A visual inspection of Leray's database shows some insertion / deletion of three nucleotides (one aminoacid).
# --p-min-size 1: We will allow singleton ASVs considering the low frequence of air DNA.
# -p-min-reads 1: We are also not going to filter by number of reads. After taxonomic asignment we will apply standardized filters.
for i in t0_l10 t20_l10 t20_l50 t26_l10
do
	echo "Denoising $i"
	qiime deblur denoise-other \
	  --i-demultiplexed-seqs ../data/03_merged_$i.qza \
	  --i-reference-seqs ../../../00_Databases/COInrLerayQIIME/COInrLerayQIIME_trainseq.qza \
	  --o-table ../data/04_denoised_table_$i.qza \
	  --o-representative-sequences ../data/04_denoised_seqs_$i.qza \
	  --o-stats ../data/04_denoised_stats_$i.qza \
	  --p-trim-length -1 \
	  --p-left-trim-len 0 \
	  --p-indel-max 3 \
	  --p-min-size 1 \
	  --p-min-reads 1 \
	  --p-sample-stats \
	  --p-jobs-to-start 6 \
	  --verbose 2>&1 | tee ./logs/03_denoising_$i.log
	  echo "Finishing denoising of $i"
done

# Moving deblur's log
mv deblur.log ./logs/03_deblur.log

# 09.- Visualizing deblur results
for i in t0_l10 t20_l10 t20_l50 t26_l10
do
	echo "Visualizing $i"
	qiime deblur visualize-stats \
	  --i-deblur-stats ../data/04_denoised_stats_$i.qza \
	  --o-visualization ../results/04_denoised_stats_$i.qzv \
	  --verbose
	
	  qiime feature-table summarize \
	    --i-table ../data/04_denoised_table_$i.qza \
	    --o-visualization ../results/04_denoised_table_$i.qzv \
	    --verbose
	    
	  qiime feature-table tabulate-seqs \
	    --i-data ../data/04_denoised_seqs_$i.qza \
	    --o-visualization ../results/04_denoised_seqs_$i.qzv \
	    --verbose
	echo "Finishing visualizing $i"
done

echo "To visualize any result use:"
ls ../results/04_* | sed 's/^/qiime tools view /'

