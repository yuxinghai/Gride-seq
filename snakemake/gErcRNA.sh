zcat ../results/03_interaction/mtx.gz |head -n 1 >enrichRNA
sed -r -e 's/\.[0-9]+\+[a-z_]*//g' enrichRNA >enrichRNA_ID.txt
zcat ../results/03_interaction/mtx.gz |cut -f1,5 >sfpq.txt

