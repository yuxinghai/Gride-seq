#!/usr/bin

## bedu is RNA-DNA interaction pairs for every cols
## bed1 is pos strand RNA of each interaction, bed0 is neg
## pkrna is peak of interacted RNA
## pkdna is bed file that divide all chrom into 1k window
bedu=$1
bed1=$2
bed0=$3
pkrna=$4
pkdna=$5
chrom_size=$6
workdir=$7
## RNA Peakcalling
zcat ${bedu} | awk 'BEGIN{OFS="\t"} $4=="+"{print $1,$2,$3,$9,0,$4}' | sort -S 10G -k1,1 -k2,2n > ${bed1}
zcat ${bedu} | awk 'BEGIN{OFS="\t"} $4=="-"{print $1,$2,$3,$9,0,$4}' | sort -S 10G -k1,1 -k2,2n > ${bed0}
macs2 callpeak --verbose 0 --broad --nolambda --nomodel --extsize 100 -g hs -n ${workdir}/RNA+ -t ${bed1}
awk 'BEGIN{OFS="\t"} $9>3{print $1,$2,$3,$1 ":" $2 "-" $3 ",+", $9, "+"}' ${workdir}/RNA+_peaks.broadPeak > RNA.pktmp
macs2 callpeak --verbose 0 --broad --nolambda --nomodel --extsize 100 -g hs -n ${workdir}/RNA- -t ${bed0}
awk 'BEGIN{OFS="\t"} $9>3{print $1,$2,$3,$1 ":" $2 "-" $3 ",-", $9, "-"}' ${workdir}/RNA-_peaks.broadPeak >> RNA.pktmp

sort -S 10G -k1,1 -k2,2n RNA.pktmp > ${pkrna}
rm RNA.pktmp
# DNA bin=1k
awk -v bin=1000 'BEGIN { OFS="\t" }
        !/chrM/{
            for(i=0;i<$2/bin-1;i++) {
                print $1, i*bin, (i+1)*bin, $1 ":" (i*bin) "-" ((i+1)*bin), 0, "+"
                    }
                }' ${chrom_size} > ${pkdna}

