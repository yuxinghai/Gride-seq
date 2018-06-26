RNA_rep1=$1
RNA_rep2=$2
pkdna=$3
chrom_size=$4
RNA_rep1_cov=$5
RNA_rep2_cov=$6
if [[ ! -e $pkdna ]]; then
    awk -v bin=1000 'BEGIN { OFS="\t" }
        !/chrM/{
            for(i=0;i<$2/bin-1;i++) {
                print $1, i*bin, (i+1)*bin, $1 ":" (i*bin) "-" ((i+1)*bin), 0, "+"
                    }
                }' ${chrom_size} > ${pkdna}
fi

bedtools bamtobed -i ${RNA_rep1} |grep -v "chrM" | sort -k1,1 -k2,2n -S 10G | \
bedtools coverage -sorted -g $chrom_size -counts -s -a $pkdna -b - > $RNA_rep1_cov
bedtools bamtobed -i ${RNA_rep2} |grep -v "chrM" | sort -k1,1 -k2,2n -S 10G | \
bedtools coverage -sorted -g $chrom_size -counts -s -a $pkdna -b - > $RNA_rep2_cov