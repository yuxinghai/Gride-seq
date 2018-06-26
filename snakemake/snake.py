"""
gride-seq workflow
"""
import os
base="/data1/zhoulab/yuxinghai/cluster_data2"
proj = base + "/Gride-seq"
work_dir = proj + "/results"
index = "/data2/zhoulab/common/index/bowtie2"
bw_dir = proj + "/results/bw"
raw_dir = proj + "/data"
script_dir = proj + "/script"
chrom_size = base + "/software/chrom_size/hg19.chrom.sizes"
gene_bed = proj + "/anno/gencodev19_Ltx.bed6"
IR_bed = base + "/parasp/PSF_IR/results_genecodev19/further_analysis/psf_IR.bed6"
ginfo = "%s/02_peak_Calling/Gene_wo_chrM.bed" % work_dir,

id = ["1","2"]
type = ["cDNA","gDNA"]
cDNA = ["hela_rep1.cDNA","hela_rep2.cDNA"]
gDNA = ["hela_rep1.gDNA","hela_rep2.gDNA"]
samples = cDNA + gDNA




rule all:
    input:
        expand("%s/01_mapping/{samples}.bam" % work_dir, samples = samples),
        expand("%s/01_mapping/rna_dna_rep{id}.bedu.gz" % work_dir, id =id),
        "%s/02_peak_Calling/RNA+.bed" % work_dir,
        "%s/03_Gross_coverage/RNA.peakgene.cvg" % work_dir,
        "%s/04_GeneHits/GeneMask.bed" % work_dir,
        "%s/05_interaction/DNA.bin1k.bg.bed" % work_dir,
        "%s/05_interaction/ghits.pkbin.net.gz" % work_dir,
        "%s/05_interaction/flt.bed" % work_dir,
        "%s/05_interaction/mtx.gz" % work_dir,
        "%s/RNA_cor/RNA_cor.pdf" % work_dir,
        "%s/06_IR_interact/IR_ghits.gep.pkavg.net.gz" % work_dir,


rule bowtie2:
    """
    maping: RNA/DNA to genome
    filtering unimapped RNA mapq > 2
    """
    input:
        fa = "%s/{samples}.fq" % raw_dir,
        index = "%s" % index,
    output:
        bam = "%s/01_mapping/{samples}.bam" % work_dir,
    shell:
        """
        bowtie2 -p 32 --local -x {input.index}/hg19 -U {input.fa} | samtools view -Sbq 2 - > {output.bam}
        """

rule unique_readmate:
    """
    GRIDseq: unique readmates
    """
    input:
        RNA = "%s/01_mapping/hela_rep{id}.cDNA.bam" % work_dir,
        DNA = "%s/01_mapping/hela_rep{id}.gDNA.bam" % work_dir,
    output:
        bedu = "%s/01_mapping/rna_dna_rep{id}.bedu.gz" % work_dir,
    shell:
        """
        bash {script_dir}/unique_read.sh {input.RNA} {input.DNA} {output.bedu} {wildcards.id}
        """

rule peakcalling:
    input:
        bedu = expand("%s/01_mapping/rna_dna_rep{id}.bedu.gz" % work_dir, id=id),
    output:
        bed1 = "%s/02_peak_Calling/RNA+.bed" % work_dir,
        bed0 = "%s/02_peak_Calling/RNA-.bed" % work_dir,
        pkrna = "%s/02_peak_Calling/RNA.peaks.bed" % work_dir,
        pkdna = "%s/02_peak_Calling/DNA.bin1k.bed" % work_dir,
    run:
        bedu_str = "\"" + " ".join(input.bedu) + "\""
        command = " ".join(["bash", script_dir+"/peak_calling.sh", bedu_str, output.bed1,
        output.bed0, output.pkrna, output.pkdna, chrom_size, work_dir+"/02_peak_Calling"])
        shell(command)


rule Gross_coverage_tx:
    input:
        pkrna = "%s/02_peak_Calling/RNA.peaks.bed" % work_dir,
        bedu = expand("%s/01_mapping/rna_dna_rep{id}.bedu.gz" % work_dir, id=id),
    output:
        gcvg = "%s/03_Gross_coverage/genecvg.bed" % work_dir,
        pkgene = "%s/03_Gross_coverage/RNA.peakgene.bed" % work_dir,
        pkcvg = "%s/03_Gross_coverage/RNA.peakgene.cvg" % work_dir,
        tx_bed = "%s/anno/Tx.bed" % proj,
    run:
        bedu_str = "\"" + " ".join(input.bedu) + "\""
        command = " ".join(["bash", script_dir+"/Gross_coverage_tx.sh", ginfo, output.gcvg,
                            chrom_size, gene_bed, output.pkgene, output.pkcvg, input.pkrna, bedu_str,
                            output.tx_bed])
        shell(command)

rule GeneMids_GeneHits:
    input:
        pkgene = "%s/03_Gross_coverage/RNA.peakgene.bed" % work_dir,
        pkcvg = "%s/03_Gross_coverage/RNA.peakgene.cvg" % work_dir,
    output:
        gmask = "%s/04_GeneHits/GeneMask.bed" % work_dir,
        gmids = "%s/04_GeneHits/GeneMids.bed" % work_dir,
        ghits = "%s/04_GeneHits/GeneHits.bed" % work_dir,
    shell:
        """
        bash {script_dir}/GeneMids_GeneHits.sh {output.gmask} {output.gmids} {output.ghits} \
        {input.pkgene} {input.pkcvg}
	    """

rule interaction_network1:
    input:
        bedu = expand("%s/01_mapping/rna_dna_rep{id}.bedu.gz" % work_dir, id=id),
        gmask = "%s/04_GeneHits/GeneMask.bed" % work_dir,
        ghits = "%s/04_GeneHits/GeneHits.bed" % work_dir,
        pkdna = "%s/02_peak_Calling/DNA.bin1k.bed" % work_dir,
        gmids = "%s/04_GeneHits/GeneMids.bed" % work_dir,
    output:
        netg = "%s/05_interaction/ghits.net.gz" % work_dir,
        netb = "%s/05_interaction/gmids.net.gz" % work_dir,
        bgdna = "%s/05_interaction/DNA.bin1k.bg.bed" % work_dir,
    message: 'bash {script_dir}/test1.sh {output.netg} {output.netb}\
            "{input.bedu}" {input.gmask} {input.ghits} {input.pkdna} {input.gmids} {output.bgdna}\
             {chrom_size}'
    run:
        bedu_str = "\"" + " ".join(input.bedu) + "\""
        command = " ".join(["bash", script_dir+"/test1.sh", output.netg, output.netb,
                    bedu_str, input.gmask, input.ghits, input.pkdna, input.gmids, output.bgdna,
                    chrom_size])
        shell(command)

rule interaction_network2:
    input:
        bgdna = "%s/05_interaction/DNA.bin1k.bg.bed" % work_dir,
        netg = "%s/05_interaction/ghits.net.gz" % work_dir,
    output:
        pknetg = "%s/05_interaction/ghits.pkbin.net.gz" % work_dir,
    message: 'bash {script_dir}/test2.sh {input.netg} {input.bgdna} {output.pknetg} {chrom_size}'
    run:
        command = " ".join(["bash", script_dir+"/test2.sh", input.netg, input.bgdna,
                    output.pknetg,chrom_size])
        shell(command)

rule interaction_network3:
    input:
        pknetg = "%s/05_interaction/ghits.pkbin.net.gz" % work_dir,
    output:
        pkinfo = "%s/05_interaction/ghits.dnapeak.info" % work_dir,
        gflts = "%s/05_interaction/flt.bed" % work_dir,
    message: 'bash {script_dir}/test3.sh {input.pknetg} {output.pkinfo} {output.gflts}'
    run:
        command = " ".join(["bash", script_dir+"/test3.sh", input.pknetg, output.pkinfo,
                    output.gflts])
        shell(command)

rule interaction_network4:
    input:
        pknetg = "%s/05_interaction/ghits.pkbin.net.gz" % work_dir,
        pkinfo = "%s/05_interaction/ghits.dnapeak.info" % work_dir,
    output:
        mtx = "%s/05_interaction/mtx.gz" % work_dir,
    message: 'bash {script_dir}/test4.sh {input.pknetg} {chrom_size} {input.pkinfo} {output.mtx}'
    run:
        command = " ".join(["bash", script_dir+"/test4.sh", input.pknetg, chrom_size, input.pkinfo,
                    output.mtx])
        shell(command)

rule RNA_cor:
    input:
        RNA_rep1 = "%s/01_mapping/hela_rep1.cDNA.bam" % work_dir,
        RNA_rep2 = "%s/01_mapping/hela_rep2.cDNA.bam" % work_dir,
        pkdna = "%s/02_peak_Calling/DNA.bin1k.bed" % work_dir,
    output:
        RNA_rep1_cov = "%s/RNA_cor/RNA_rep1_cov.bed" % work_dir,
        RNA_rep2_cov = "%s/RNA_cor/RNA_rep2_cov.bed" % work_dir,
        pdf = "%s/RNA_cor/RNA_cor.pdf" % work_dir,
    run:
        command1 = " ".join(["bash", script_dir+"/gRNA_reads_rep.sh", input.RNA_rep1, input.RNA_rep2, 
        input.pkdna, chrom_size, output.RNA_rep1_cov, output.RNA_rep2_cov])
        command2 = " ".join(["Rscript", script_dir+"/RNA_cor.R", output.RNA_rep1_cov, output.RNA_rep2_cov, output.pdf])
        shell(command1)
        shell(command2)

rule IR_gene_interaction:
    input:
        gsize = chrom_size,
        IR_reg = IR_bed,
        ghits = "%s/04_GeneHits/GeneHits.bed" % work_dir,
        pknetg = "%s/05_interaction/ghits.pkbin.net.gz" % work_dir,
    output:
        dtgt="%s/06_IR_interact/IR_dtgt.bed" % work_dir,
        netx="%s/06_IR_interact/IR_ghits.gep.pkmax.net.gz" % work_dir,
	    netm="%s/06_IR_interact/IR_ghits.gep.pkavg.net.gz" % work_dir,
        maxbed="%s/06_IR_interact/IR_ghits.gep.pkmax.bed" % work_dir,
	    avgbed="%s/06_IR_interact/IR_ghits.gep.pkavg.bed" % work_dir,
        
    run:
        command1 = " ".join(["bash", script_dir+"/IR_Gross_coverage.sh", input.gsize, input.IR_reg, 
        input.ghits, input.pknetg, output.dtgt, output.netx, output.netm, output.maxbed, output.avgbed])
        shell(command1)

