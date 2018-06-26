#!/usr/bin

## RNA covert matched DNA, SO bedu is RNA-DNA interaction pairs for every cols
RNA=$1
DNA=$2
bedu=$3
id=$4
bedtools bamtobed -i ${RNA} | \
awk 'BEGIN{OFS="\t"} {$6 = $6=="+"? "-":"+"; print}' > tmp${id}
bedtools bamtobed -i ${DNA} >> tmp${id}
awk 'BEGIN { OFS = "\t" } {if(K[$4]) print K[$4], $0; else K[$4]=$0;}' tmp${id} | \
awk 'length($1)<6 && length($7)<6{
			x = $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $12;
			D[x] = D[x]? D[x] "," $4 : $4 
		}
		END { OFS = "\t"; for(x in D) print x, D[x]} ' | gzip -c > ${bedu}
	
rm tmp${id}
