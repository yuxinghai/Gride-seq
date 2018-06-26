#!/usr/bin
gmask=$1
gmids=$2
ghits=$3
pkgene=$4
pkcvg=$5
## GeneMids and GeneHits (gene)　
# 1.remove the genes including in another one（de duplilicate）,the smaller
## ghits is chromatin-enriched RNAs RNA RPK (reads per Kb) ≥100, DNA RPK >10
## gmask is 
bedtools intersect -sorted -s -wo -a $pkgene -b $pkgene | \
awk '$4!=$10&&$3-$2==$13&&$4!~/protein_coding/' | \
cut -f1-6 | sort -uk1,1 -k2,2n > $gmask

awk '$4~/protein_coding/' $pkgene > $gmids

awk -v fmsk=$gmask '
BEGIN {OFS="\t"; trpk=100; mrpk=10; while((getline < fmsk)>0) G[$4]++}
!G[$4] && ($7/($3-$2)*1000>=trpk || $9/($3-$2)*1000>=mrpk){print $1,$2,$3,$4,$5,$6}
' $pkcvg > $ghits