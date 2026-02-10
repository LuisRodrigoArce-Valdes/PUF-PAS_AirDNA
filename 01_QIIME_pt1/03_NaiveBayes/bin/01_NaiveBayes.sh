#!/bin/sh
# 01_Taxonomy.sh
# 21/10/2025
# We will use this script to asign taxonomies to the sequences we have denoised using both Deblur and DADA2.
# https://cryptick-lab.github.io/NGS-Analysis/_site/QIIME2-Taxonomy.html
  
# 01.- Creating symbolic links to the denoised databases from Deblur and DADA2
echo "Beginning with Deblur!"
for i in t0_l10 t20_l10 t20_l50 t26_l10
do
	ln -s -f ../../01_Deblur/data/04_denoised_seqs_$i.qza ../data/01_Deblur_$i.qza
done

echo "Now with DADA2"
for i in t0_l10 t20_l10 t20_l50 t26_l10
do
	ln -s -f ../../02_DADA2/data/02_denoised_seqs_$i.qza ../data/01_DADA2_$i.qza
done

# 02.- Classifing using the trained classifier
# --p-confidence: confidence treshold for limiting taxonomic depth. Set 0 to estimate without using p-condifence.
for i in Deblur DADA2
do
	for n in t0_l10 t20_l10 t20_l50 t26_l10
	do
		echo "Classifing ${i}_${n}"
		qiime feature-classifier classify-sklearn \
		  --i-classifier ../../../00_Databases/COInrLerayQIIME/COInrLeray_bayes_classifier.qza \
		  --i-reads ../data/01_${i}_${n}.qza \
		  --p-n-jobs 6 \
		  --p-confidence 0 \
		  --o-classification ../data/02_Bayes_${i}_${n}.qza \
		  --verbose
	done
done

# 03.- Visualizing results
for i in Deblur DADA2
do
	for n in t0_l10 t20_l10 t20_l50 t26_l10
	do
		qiime metadata tabulate \
		  --m-input-file ../data/02_Bayes_${i}_${n}.qza \
		  --o-visualization ../results/01_Bayes_${i}_${n}_taxa.qzv
		echo "Visualize using: qiime tools view ../results/01_Bayes_${i}_${n}_taxa.qzv"
	done
done

# 04.- Exporting results
# Creating symbolic links to the denoised tables from Deblur and DADA2
echo "Beginning with Deblur!"
for i in t0_l10 t20_l10 t20_l50 t26_l10
do
	ln -s -f ../../01_Deblur/data/04_denoised_table_$i.qza ../data/01_Deblur_table_$i.qza
done

echo "Now with DADA2"
for i in t0_l10 t20_l10 t20_l50 t26_l10
do
	ln -s -f ../../02_DADA2/data/02_denoised_table_$i.qza ../data/01_DADA2_table_$i.qza
done

mkdir -p ../results/03_Matrices/

for i in Deblur DADA2
do
	for n in t0_l10 t20_l10 t20_l50 t26_l10
	do
		qiime tools export \
		  --input-path ../data/01_${i}_table_${n}.qza \
		  --output-path ../data/
		
		biom convert \
		--to-tsv \
		-i ../data/feature-table.biom \
		-o ../results/02_feature-table_${i}_${n}.tsv
		
		rm ../data/feature-table.biom
		
		qiime tools export \
		  --input-path ../results/01_Bayes_${i}_${n}_taxa.qzv \
		  --output-path ../results/03_Matrices/${i}_${n}
		  
	done
done

