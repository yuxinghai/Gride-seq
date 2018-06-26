pknetg=$1
pkinfo=$2
gflts=$3

## merge each part to broad peak
zcat $pknetg | awk -v gap=10 '
gid != $4 || $7 != chr || $8-pe > gap*1000 {
    vp = chr ":" ps ":" pe
    if(NR > 1) D[gid] = D[gid] ? D[gid] "," vp : vp
    ps = $8
}
{ 	gid = $4; chr = $7; pe = $8+1000;
    G[gid] = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
}
END {	OFS = "\t"
    vp = chr ":" ps ":" pe
    D[gid] = D[gid] ? D[gid] "," vp : vp

    for(k in D) print G[k], D[k]
}
' | sort -k1,1 -k2,2n -k4,4 | \
tee $pkinfo | cut -f1-6 > $gflts
