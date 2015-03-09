#!/bin/bash
#script will count a genomic transitions on first nucleotide from each read. Input file must be in SAM format.

total=$(cat $1 | grep -v @ | wc -l)

#same strand
sameA=$(cat $1 | awk '{if($2=="0") print $18}' | grep :0A | wc -l)
sameC=$(cat $1 | awk '{if($2=="0") print $18}' | grep :0C | wc -l)
sameT=$(cat $1 | awk '{if($2=="0") print $18}' | grep :0T | wc -l)
sameG=$(cat $1 | awk '{if($2=="0") print $18}' | grep :0G | wc -l)

#anti strand
antiA=$(cat $1 | awk '{if($2=="16") print $18}' | grep :T0 | wc -l)
antiC=$(cat $1 | awk '{if($2=="16") print $18}' | grep :G0 | wc -l)
antiT=$(cat $1 | awk '{if($2=="16") print $18}' | grep :A0 | wc -l)
antiG=$(cat $1 | awk '{if($2=="16") print $18}' | grep :C0 | wc -l)

echo transitions of first genomic nts > $2
echo A:   $((($sameA+$antiA))) >> $2
echo C:   $((($sameC+$antiC))) >> $2
echo T:   $((($sameT+$antiT))) >> $2
echo G:   $((($sameG+$antiG))) >> $2
echo number of all candidates: $total >> $2
