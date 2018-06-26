netg=$1
bgdna=$2
pknetg=$3
gsize=$4

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