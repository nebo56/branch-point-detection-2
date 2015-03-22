#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash


PYTHONPATH=/home/skgthab/programs/Python-2.7.5/bin
PERL5LIB=/home/skgthab/programs/ActivePerl-5.16.3.1603-x86_64-linux-glibc-2.3.5-296746/perl/bin
export PATH=$PYTHONPATH:$PATH
export PATH=$PERL5LIB:$PATH
export PATH=/home/skgthab/programs/tophat-2.0.9.Linux_x86_64:$PATH
export PATH=/home/skgthab/programs/bowtie2-2.1.0:$PATH
export PATH=/home/skgthab/programs/samtools-0.1.19:$PATH
export PATH=/home/skgthab/programs/fastx_toolkit0.0.13:$PATH
export PATH=/home/skgthab/programs/bedtools-2.17.0/bin:$PATH

data=$1
path=/cluster/project9/ule-group/BranchPoints/branch-point-detection-2/
introns=/home/skgthab/annotations/regions/mm9-introns.bed

# unzip
gunzip ${path}${data}.fq.gz

# clip the adapter and discard non-clipped sequences and discard the sequences that are shorter then 15 nt + 5 random barcode + 4 experimental barcode
fastx_clipper -Q 33 -a AGATCGGAAG -c -n -l 24 -i  ${path}${data}.fq -o  ${path}${data}-clipped.fq

# fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-clipped.fq -o ${path}${data}-clipped.fa
rm ${path}${data}-clipped.fq

# remove the barcode and reads that are shorter then 17 nts
python ${path}swap_barcodes_to_header.py ${path}${data}-clipped.fa ${path}${data}-noBarcodes.fa
rm ${path}${data}-clipped.fa

# map on genome
bowtie2-align -x ~/bowtie-indexes/mm9/mm9 -f ${path}${data}-noBarcodes.fa -S ${path}${data}.sam
rm ${path}${data}-noBarcodes.fa

# filter reads woth more then 2 mismatches
samtools view -Sh ${path}${data}.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-2mis.sam
rm ${path}${data}.sam

# count of genomic transitions from all reads
python ${path}transition_ratio.py ${path}${data}-2mis.sam ${path}${data}-genomic_transitions-all.log

# keep only reads that ends with AG
python ${path}get_branch_point_candidates.py ${path}${data}-2mis.sam ${path}${data}-filtered.sam
rm ${path}${data}-2mis.sam

# count of genomic transitions from selected reads that ends with AG
python ${path}transition_ratio.py ${path}${data}-filtered.sam ${path}${data}-genomic_transitions-filteredAG.log

# trim SAM reads that starts with A mutation on genome
python ${path}trimSAM.py ${path}${data}-filtered.sam ${path}${data}-trimmed.sam
rm ${path}${data}-filtered.sam

# SAM to BAM
samtools view -hSb ${path}${data}-trimmed.sam > ${path}${data}-trimmed.bam
rm ${path}${data}-trimmed.sam

# remove duplicates (reads with the same random barcodes)
bedtools bamtobed -i ${path}${data}-trimmed.bam > ${path}${data}-trimmed.bed
sort -k1,1 -k2,2n -k6,6 ${path}${data}-trimmed.bed | uniq > ${path}${data}-trimmed-uniq.bed
rm ${path}${data}-trimmed.bed

# report reads which overlaps 100% with introns (intron .bed is from UCSC bed introns only)
bedtools intersect -s -f 1.00 -a ${path}${data}-trimmed-uniq.bed -b ${introns} | uniq > ${path}${data}-trimmed-uniq-introns.bed
rm ${path}${data}-trimmed-uniq.bed

# flank intron border on the same and anti strand
python ${path}flankBEDpositions.py ${introns} ${introns}-flanked_0_-2.bed 0 -2

# remove all reads that are 100% in flanked introns which means they are not next to the exon position
bedtools intersect -s -v -f 1.00 -a ${path}${data}-trimmed-uniq-introns.bed -b ${introns}-flanked_0_-2.bed | uniq > ${path}${data}-trimmed-uniq-introns-selected.bed
#rm ${introns}-flanked_0_-2.bed
rm ${path}${data}-trimmed-uniq-introns.bed

# set end positions of the read
python ${path}set_branch_point_position.py ${path}${data}-trimmed-uniq-introns-selected.bed ${path}${data}-trimmed-uniq-introns-selected-bp.bed
#rm ${path}${data}-trimmed-uniq-introns-selected.bed

# sum together reads that ends on the same position
cat ${path}${data}-trimmed-uniq-introns-selected-bp.bed | awk '{print $1 "\t" $2 "\t" $3 "\t\t\t" $6}' | sort -k1,1 -k2,2n -k6,6 > ${path}${data}-trimmed-uniq-introns-selected-bp-sorted.bed
python ${path}BEDsum.py ${path}${data}-trimmed-uniq-introns-selected-bp-sorted.bed ${path}${data}-branch_points.bed
rm ${path}${data}-trimmed-uniq-introns-selected-bp.bed
rm ${path}${data}-trimmed-uniq-introns-selected-bp-sorted.bed

# compress and move results to a new folder
gzip ${path}${data}.fq
mkdir ${path}${data}
mv ${path}${data}* ${path}${data}
