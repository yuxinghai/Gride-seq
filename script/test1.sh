#!/usr/bin
netg=$1
netb=$2
bedu=$3
gmask=$4
ghits=$5
pkdna=$6
gmids=$7
bgdna=$8
gsize=$9
## RNA-DNA network (gene)
# 1. bedu not overlap with gmask but overlap with ghits

zcat $bedu | awk 'BEGIN{OFS="\t"} {print $1,$2,$3, $5 "@" $6 "@" $7 "@" $8, 0,$4}' | \
sort -k1,1 -k2,2n | bedtools intersect -sorted -s -v -f 0.5 -wa -a - -b $gmask | \
bedtools intersect -sorted -s -wa -wb -a - -b $ghits | \
awk 'BEGIN{OFS="\t"} {split($4,X,"@"); split($10,Y,"|");
    print X[1],X[2],X[3], $7 "@" $8 "@" $9 "@" Y[1] "@" $11 "@" $12 ,$11,X[4]}' | \
sort -k1,1 -k2,2n | bedtools intersect -sorted -wa -wb -a - -b $pkdna | \
awk '{rid=$4; did=$7 "@" $8; D[rid "~" did]++}
END{ OFS="\t";
    for(k in D) {
        split(k, X, "~");
        gsub("@", "\t", X[1]); gsub("@", "\t", X[2]);
        print X[1], X[2], D[k]
    }
}' | sort -k1,1 -k2,2n -k7,7 -k8,8n | gzip -c > $netg


zcat $bedu | awk 'BEGIN{OFS="\t"} {print $1,$2,$3, $5 "@" $6 "@" $7 "@" $8, 0,$4}' | \
sort -k1,1 -k2,2n | bedtools intersect -sorted -s -v -f 0.5 -wa -a - -b $gmask | \
bedtools intersect -sorted -s -wa -wb -a - -b $gmids | \
awk 'BEGIN{OFS="\t"} {split($4,X,"@"); split($10,Y,"|");
    print X[1],X[2],X[3], $7 "@" $8 "@" $9 "@" Y[1] "@" $11 "@" $12 ,$11,X[4]}' | \
sort -k1,1 -k2,2n | bedtools intersect -sorted -wa -wb -a - -b $pkdna | \
awk '{rid=$4; did=$7 "@" $8; D[rid "~" did]++}
END{ OFS="\t";
    for(k in D) {
        split(k, X, "~");
        gsub("@", "\t", X[1]); gsub("@", "\t", X[2]);
        print X[1], X[2], D[k]
    }
}' | sort -k1,1 -k2,2n -k7,7 -k8,8n | gzip -c > $netb


## Chromosome Backgrounds by mid genes
## Smooth bin (dm 10k, hs 100k)
## v is foldchange
w=100

zcat $netb | awk -v pkdna=$pkdna -v g=$gsize -v w=$w '
BEGIN{ OFS="\t";
    while((getline < g)>0) { C[$1] = $2 }; close(g)
}
$1 != $7{
    did = $7 ":" $8; D[did] += $9; N[$7] += $9
}
END{
    while((getline < pkdna)>0) {
        x = 0
        for(i=-w/2; i<w/2; i++) {
            j = ($2+i*1000)
            did = $1 ":" j
            v = D[did]/N[$1] * C[$1]/1e3
            x += v/w
        }
        print $1, $2, $3, $4, x, $6
    }
}' > $bgdna
