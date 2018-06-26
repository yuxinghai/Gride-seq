pknetg=$1
gsize=$2
pkinfo=$3
mtx=$4
zcat $pknetg | sort -k7,7 | awk -v g=$gsize -v pkinfo=$pkinfo '
BEGIN{ OFS="\t";
    while((getline < g)>0) C[$1] = $2; close(g)
    while((getline < pkinfo)>0) { gids = gids ? gids "\t" $4 : $4 }; close(hit)
    print "Chrom", gids
    ylen = split(gids, G, "\t")
}
NR>1 && gchr != $7 {
    xlen = int(C[gchr]/1000)

    for(i=0; i<xlen; i++) {
        vstr = gchr ":" (i*1000)
        for(j=1; j<=ylen; j++) {
            k = G[j] "|" i
            x = ( V[k] ? int(V[k]*100)/100 : 0 )
            vstr = vstr "\t" x
        }
        print vstr
    }

    delete V
}
{ V[$4 "|" ($8/1000)] = $9; gchr = $7 }
END {
    xlen = int(C[gchr]/1000)

    for(i=0; i<xlen; i++) {
        vstr = gchr ":" (i*1000)
        for(j=1; j<=ylen; j++) {
            k = G[j] "|" i
            x = ( V[k] ? int(V[k]*100)/100 : 0 )
            vstr = vstr "\t" x
        }
        print vstr
    }
}
' | gzip -c > $mtx