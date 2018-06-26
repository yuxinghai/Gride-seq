## 20180625
## gene ~ IR region network, source from gene ~ enhancer, promoter network
# input
gsize=$1
IR_reg=$2
hits=$3
pknet=$4
# ouput
dtgt=$5
netx=$6
netm=$7
maxbed=$8
avgbed=$9

sort -k1,1 -k2,2n $hits > tmp 
grep -v "Chr" $IR_reg | awk 'BEGIN{OFS="\t"} {
    k = $4 ":" $1 ":" $2 "-" $3
    print $1,$2,$3,k,$5,$6
}'| sort -k1,1 -k2,2n  >> tmp
sort -u -k1,1 -k2,2n -k3,3n tmp > $dtgt

zcat $pknet | awk 'BEGIN{OFS="\t"} {print $7,$8,$8+1000,$4,$9,"."}' | \
sort -k1,1 -k2,2n | bedtools intersect -sorted -wa -wb -a - -b $dtgt | \
awk '{
    k = $4 "\t" $10 "#" $7 ":" $8 "-" $9
    v = $5/int(1+($9-$8)/1000)
    S[k] = $5 "\t" $11
    D[k] += v
}
END{ OFS="\t"; 
    for(k in D) {
        split(k, X, "#"); gv = D[k]; gk = X[1]
        E[gk] = gv > E[gk]? gv:E[gk]
    }
    for(k in E) print k, E[k]
}' | sort -k1,1 -k2,2 | awk '$1!=$2 {print $0}'| gzip -c > $netx
zcat $netx |grep "/ENSG" > $maxbed

zcat $pknet | awk 'BEGIN{OFS="\t"} {print $7,$8,$8+1000,$4,$9,"."}' | \
sort -k1,1 -k2,2n | bedtools intersect -sorted -wa -wb -a - -b $dtgt | \
awk '{
    k = $4 "\t" $10 "#" $7 ":" $8 "-" $9
    v = $5/int(1+($9-$8)/1000)

    D[k] += v
}
END{ OFS="\t"; 
    for(k in D) {
        split(k, X, "#"); gv = D[k]; gk = X[1]
        E[gk] += gv; N[gk] ++
    }
    for(k in E) print k, E[k]/N[k]
}' | sort -k1,1 -k2,2 |  awk '$1!=$2 {print $0}'|  gzip -c > $netm
zcat $netm |grep "/ENSG" > $avgbed