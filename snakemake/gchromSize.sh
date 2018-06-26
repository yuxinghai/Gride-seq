#!/usr/bin/bash
 grep -v "_" /data2/zhoulab/yuxinghai/software/chrom_size/hg19.chrom.sizes_bcup|grep -v "chrM" |sort -k 1,1 -k2,2n >/data2/zhoulab/yuxinghai/software/chrom_size/hg19.chrom.sizes
