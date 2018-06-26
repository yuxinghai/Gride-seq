#!/usr/bin
netg=$1
netb=$2
bedu=$3
gmask=$4
ghits=$5
pkdna=$6
gmids=$7
bgdna=$8
pknetg=$9
gsize=$10
pkinfo=$11
gflts=$12
mtx=$13
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


## RNA-DNA network normalize and peakfilter

zcat $netg | awk -v g=$gsize -v fb=$bgdna -v win=10 -v nf=2 -v mwin=3 '
BEGIN{ OFS="\t";
    while((getline < g)>0) {nc += $2; G[$1]=$2}; close(g)
    while((getline < fb)>0) B[$1 ":" $2] = $5; close(fb)
}
gid != $4 || $7 != chr || $8-pe>win*2000 {
    delete D
    if (NR>1) {
        for(p=ps; p<=pe; p+=1000) {
            nwl = 0; nwm = 0; nwr = 0
            for(i=-win; i<-win/2; i++) nwl += R[p+i*1000] > nf
            for(i=-win/2; i<win/2; i++) nwm += R[p+i*1000] > nf
            for(i=win/2; i<win; i++) nwr += R[p+i*1000] > nf

            R[p] = nwl > mwin || nwm > mwin || nwr > mwin ? R[p]:0
        }
        for(p=ps-win*500; p<pe+win*500; p+=1000) {
            vf = 0
            for(i=-win/2; i<win/2; i++)
                vf += R[p+i*1000]/win
            D[p] = vf
        }
        for(p=ps-win*1000; p<pe+win*1000; p+=1000) {
            vf = 0
            for(i=-win/2; i<win/2; i++)
                vf += D[p+i*1000]/win
            if(vf > 0 && p>=0 && p<G[chr])
                print ginfo, chr, p, vf
        }
    }
    ps = $8; delete R
}
{
    gid = $4; chr = $7; pe = $8+1000
    ginfo = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
    x = ($9/$5 * nc/1e3 + 1) / (B[$7 ":" $8] + 1)
    if(x > nf) R[$8] = x
}
END {
    delete D
    for(p=ps; p<=pe; p+=1000) {
        nwl = 0; nwm = 0; nwr = 0
        for(i=-win; i<-win/2; i++) nwl += R[p+i*1000] > nf
        for(i=-win/2; i<win/2; i++) nwm += R[p+i*1000] > nf
        for(i=win/2; i<win; i++) nwr += R[p+i*1000] > nf

        R[p] = nwl > mwin || nwm > mwin || nwr > mwin ? R[p]:0
    }
    for(p=ps-win*500; p<pe+win*500; p+=1000) {
        vf = 0
        for(i=-win/2; i<win/2; i++)
            vf += R[p+i*1000]/win
        D[p] = vf
    }
    for(p=ps-win*1000; p<pe+win*1000; p+=1000) {
        vf = 0
        for(i=-win/2; i<win/2; i++)
            vf += D[p+i*1000]/win
        if(vf > 0 && p>=0 && p<G[chr])
            print ginfo, chr, p, vf
    }
}
' | gzip -c > $pknetg

## peak regions

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


## gene hits network --> matrix

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
