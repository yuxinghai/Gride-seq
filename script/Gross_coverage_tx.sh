#!/usr/bin
#### not consider the enriched RNA not overlap with gene and distance to gene <1000bp

## gene is bed6 file with largest tx of a gene, $4=(gene name + gene type)
## ginfo is bed6 file without chrM but with largest tx of a gene, $4=(gene name + gene type) 
## gcvg is interacted RNA coverage in each gene
## pkrna is peak of interacted RNA
## pkgene is peak of interacted RNA mathched gene with coverage
## kt is rid num, ks is rid=rid num, km is num of rid eq did 

ginfo=$1 
gcvg=$2
gsize=$3
gene=$4  # bed6 format
pkgene=$5
pkcvg=$6
pkrna=$7
bedu=$8
tx_bed=$9
echo "####variable########"
echo "ginfo:"$ginfo
echo "gcvg:"$gcvg
echo "gsize:"$gsize
echo "gene:"$gene
echo "pkgene:"$pkgene
echo "pkcvg:"$pkcvg
echo "pkrna:"$pkrna
echo "bedu:"$bedu
echo "#####variable#######"

script_dir=/data1/zhoulab/yuxinghai/cluster_data2/Gride-seq/script
if [[ ! -e $ginfo ]]; then
	awk -v g=$gsize ' BEGIN{OFS="\t"; while((getline < g)>0) G[$1] = $2}
	$1!="chrM" && G[$1]{print}' $gene | sort -k1,1 -k2,2n -S 10G > $ginfo
fi

zcat $bedu | awk 'BEGIN{OFS="\t"} !/chrM/{print $1,$2,$3,$9,0,$4}' | \
sort -k1,1 -k2,2n -S 10G | bedtools coverage -sorted -g $gsize -counts -s -a $ginfo -b - | \
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$7,$6}' > $gcvg

# ## Active gene and Hit genes
# ##  cf achive (Ni − Ni+n) ≥ n
cf=$( Rscript ${script_dir}/gridRsub.getGeneHitCutoff.R $gcvg )
echo "##################"
## Filtering the genes with peaks
awk '$5>5 && !/chrM/' $pkrna | sort -k1,1 -k2,2n -S 10G | \
bedtools closest -s -d -a stdin -b $ginfo | \
awk 'BEGIN{OFS="\t"}
    $13!=0 && $3-$2>1000{x++; print $1,$2,$3,"Tx." x,$5,$6}
    $13==0{print $7,$8,$9,$10,$11,$12}' | \
sort -k1,1 -k2,2n -S 10G | uniq > $pkgene
echo "##################"
## Calculating the gross coverage of gene
echo $cf
zcat $bedu | awk 'BEGIN{OFS="\t"} !/chrM/{print $1,$2,$3,$9,0,$4}' | \
sort -k1,1 -k2,2n -S 10G | bedtools coverage -sorted -g $gsize -counts -s -a $pkgene -b - |
awk -v cf=$cf 'BEGIN{OFS="\t"} $7>cf {print $1,$2,$3,$4,$7,$6}' > $pkcvg
mv $pkcvg $pkgene
echo "##################"
## Calculating the point coverage of gene
zcat $bedu | awk 'BEGIN{OFS="\t"} !/chrM/{print $1,$2,$3, $5 "@" $6 "@" $7 "@" $8, 0,$4}' | \
sort -k1,1 -k2,2n -S 10G | bedtools closest -s -d -a - -b $pkgene | \
awk 'BEGIN{OFS="\t"}
$13==0{	split($4,X,"@");
    print X[1],X[2],X[3], $7 "@" $8 "@" $9 "@" $10 "@" $11 "@" $12 ,0,X[4]}' | \
sort -k1,1 -k2,2n -S 10G | bedtools closest -d -a - -b $pkgene | \
awk '$13==0{rid=$4; did=$7 "@" $8 "@" $9 "@" $10 "@" $11 "@" $12;
    X[rid "~" did]++; T[rid]++
    S[rid] += rid == did ? 1:0
}
END{ OFS="\t";
    for(i in X) {split(i, K, "~"); M[K[1]] = M[K[1]] > X[i] ? M[K[1]] : X[i] }
    for(k in T) {
        x = k; gsub("@","\t",x);
        kt = T[k]
        ks = S[k] > 0 ? S[k] : 0
        km = M[k] > 0 ? M[k] : 0
        print x, kt, ks, km
    }
}' | sort -k1,1 -k2,2n -S 10G > $pkcvg

grep "Tx." $pkgene > $tx_bed



