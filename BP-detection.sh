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
introns=mm9-introns.bed

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

# count of genomic transitions from all reads
python ${path}transition_ratio.py ${path}${data}-2mis.sam ${path}${data}-genomic_transitions-all.log

# keep only reads that ends with AG
python ${path}get_branch_point_candidates.py ${path}${data}-2mis.sam ${path}${data}-filtered.sam
#rm ${path}${data}.sam

# count of genomic transitions from selected reads that ends with AG
python ${path}transition_ratio.py ${path}${data}-filtered.sam ${path}${data}-genomic_transitions-filteredAG.log

# trim SAM reads that starts with A mutation on genome
python ${path}trimSAM.py ${path}${data}-filtered.sam ${path}${data}-trimmed.sam
rm ${path}${data}-filtered.sam

# SAM to BAM
samtools view -hSb ${path}${data}-trimmed.sam > ${path}${data}-trimmed.bam

# remove duplicates (reads with the same random barcodes)
bedtools bamtobed -i ${path}${data}-trimmed.bam > ${path}${data}-trimmed.bed
sort -k1,1 -k2,2n -k6,6 ${path}${data}-trimmed.bed | uniq > ${path}${data}-trimmed-uniq.bed

# SAM to BED with collapsed read count by random barcodes
#python ${path}SAMtoCollapsedSAMandBED.py ${path}${data}-trimmed.sam ${path}${data}-Barcodes.fa ${path}${data}-collapsed.sam ${path}${data}.bed
#python ${path}SAMtoCollapsedSAMandBED.py ${path}${data}.sam ${path}${data}-Barcodes.fa ${path}${data}-collapsed.sam ${path}${data}.bed
#rm ${path}${data}-trimmed.sam
#rm ${path}${data}-collapsed.sam
#rm ${path}${data}-Barcodes.fa



# we need to split bed by strand because bedtools intersect -s is not working properly
#awk '{if($6=="+") print}' ${path}${introns} > ${path}${introns}-same.bed
#awk '{if($6=="-") print}' ${path}${introns} > ${path}${introns}-anti.bed
#awk '{if($5=="+") print}' ${path}${data}.bed > ${path}${data}-same.bed
#awk '{if($5=="-") print}' ${path}${data}.bed > ${path}${data}-anti.bed
#rm ${path}${data}.bed

# report reads which overlaps 100% with introns (intron .bed is from UCSC bed introns only)
bedtools intersect -f 1.00 -a ${path}${data}-trimmed-uniq.bed -b ${path}${introns} | uniq > ${path}${data}-trimmed-uniq-introns.bed

#bedtools intersect -f 1.00 -a ${path}${data}-same.bed -b ${path}${introns}-same.bed | uniq > ${path}${data}-same-introns.bed
#bedtools intersect -f 1.00 -a ${path}${data}-anti.bed -b ${path}${introns}-anti.bed | uniq > ${path}${data}-anti-introns.bed
#rm ${path}${data}-same.bed
#rm ${path}${data}-anti.bed
#rm ${path}${introns}-same.bed
#rm ${path}${introns}-anti.bed

# flank intron border on the same and anti strand
python ${path}flankBEDpositions.py ${path}${introns} ${path}${introns}-flanked_2_0.bed 2 0

#python ${path}flankBEDpositions.py ${path}${introns} ${path}${introns}-flanked-same.bed 0 -2
#python ${path}flankBEDpositions.py ${path}${introns} ${path}${introns}-flanked-anti.bed 2 0

# remove all reads that are 100% in flanked introns which means they are not next to the exon position
bedtools intersect -v -f 1.00 -a ${path}${data}-trimmed-uniq-introns.bed -b ${path}${introns}-flanked_2_0.bed | uniq > ${path}${data}-trimmed-uniq-introns-selected.bed

#bedtools intersect -v -f 1.00 -a ${path}${data}-same-introns.bed -b ${path}${introns}-flanked-same.bed | uniq > ${path}${data}-selected_reads-same.bed
#bedtools intersect -v -f 1.00 -a ${path}${data}-anti-introns.bed -b ${path}${introns}-flanked-anti.bed | uniq > ${path}${data}-selected_reads-anti.bed
#rm ${path}${data}-same-introns.bed
#rm ${path}${data}-anti-introns.bed
#rm ${path}${introns}-flanked-same.bed
#rm ${path}${introns}-flanked-anti.bed

#merge both strands together
#cat ${path}${data}-selected_reads-same.bed ${path}${data}-selected_reads-anti.bed > ${path}${data}-selected_reads.bed
#rm ${path}${data}-selected_reads-same.bed
#rm ${path}${data}-selected_reads-anti.bed

# set end positions of the read
python ${path}set_branch_point_position.py ${path}${data}-trimmed-uniq-introns-selected.bed ${path}${data}-branch_points.bed
#rm ${path}${data}-selected_reads.bed

# SAM to BAM and BAI
#samtools view -Sb ${path}${data}-collapsed.sam > ${path}${data}-collapsed.bam
#rm ${path}${data}-collapsed.sam
#samtools sort ${path}${data}-collapsed.bam ${path}${data}-collapsed-sorted
#rm ${path}${data}-collapsed.bam
#samtools index ${path}${data}-collapsed-sorted.bam ${path}${data}-collapsed-sorted.bam.bai



