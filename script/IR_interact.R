library(stringr)
base_dir <-"~/mount_dir/rna3_data1/cluster_data2/Gride-seq/results/06_IR_interact/"
setwd(base_dir)
# mk a datafram
IR_inter <- read.table("IR_ghits.gep.pkavg.bed", header = F, stringsAsFactors = F)
GID <-substr(IR_inter$V1,1,15)
inter_pkN <-IR_inter$V1[332:length(IR_inter$V1)]
IR_inter$V1[332:length(IR_inter$V1)] <-paste(paste(inter_pkN,"+",sep = ""), inter_pkN, sep = "")
G_type <-IR_inter$V1 %>% strsplit("\\+") %>%
  sapply( "[", 2 )
IRID <-IR_inter$V2 %>% strsplit("/") %>%
  sapply( "[", 2 ) %>% substr(start = 1, stop = 15)
IRN <-IR_inter$V2 %>% strsplit("/") %>%
  sapply( "[", 1) 
IRR <-IR_inter$V2 %>% strsplit("/") %>%
  sapply( "[", 3) 
# get transcript info
tx <- read.table("../../anno/gencodev19_Ltx.bed6", header = F, stringsAsFactors = F)
tx$tx_id <- tx$V4 %>% strsplit("\\+") %>%
  sapply( "[", 1 ) %>% substr(start = 1, stop = 15)
tx_manual <- read.table("../../anno/Tx.bed", header = F, stringsAsFactors = F)
new_tx <-tx[,c(1:3,6,7)]
colnames(new_tx) <-c("chr","start","end","strand","tx_id")
tx_m <-tx_manual[,c(1:3,6,4)]
colnames(tx_m) <-colnames(new_tx)
new_tx1 <- rbind(new_tx,tx_m)

IR_df <- data.frame(geneID=GID, genetype=G_type, IRID=IRID, IRN=IRN, IR_RPK=IR_inter$V3,
                    IR_region=IRR,stringsAsFactors = F)
IR_df<-IR_df[order(-IR_df$IR_RPK),]
gGeneInterLen <-function(gene1,gene2) {
  G1 <-new_tx1[new_tx1$tx_id==gene1,]
  G2 <-new_tx1[new_tx1$tx_id==gene2,]
  if (G1$chr == G2$chr) {
    len <-max(G1$start,G2$start) - min(G1$end,G2$end)
    if (len < 0) {
      len = 0 }
  } else {
    len <- 10e10}
  return(len)
} 
geneInterLen <-c()
for (i in 1:nrow(IR_df)){
  len <-gGeneInterLen(IR_df$geneID[i],IR_df$IRID[i])
  geneInterLen <-c(geneInterLen,len)
  
}

IR_df$geneInterLen <-geneInterLen
write.table(IR_df,"IRinter.tsv",quote = F,row.names = F,sep = "\t")
IR_df_dup<-IR_df[IR_df$geneID != IR_df$IRID,]
write.table(IR_df_dup,"IRinter_duplicate.tsv",quote = F,row.names = F,sep = "\t")

######## gene-geneLen  iteract Intensity ############
IR_plot <- IR_df_dup[IR_df_dup$geneInterLen < 10e8 & IR_df_dup$geneInterLen > 0,]
library(tidyverse)
library(ggplot2)
library(reshape2)
p <-ggplot(IR_plot,(aes(x=log10(geneInterLen),y=IR_RPK))) +
  geom_point(size=0.5)+
  geom_smooth(aes(x=log10(geneInterLen),y=IR_RPK, colour="#008B8B"))+
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major=element_line(colour=NA),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=10),
        axis.text.y = element_text(face="bold",  color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log10 (gene-gene length)", y= "IR region interact intensity")
pdf("gene-geneLenVS_iteractIntensity.pdf")
print(p)
dev.off()


######## expression level interact intensity ########
base_dir2 <- "~/mount_dir/rna3_data1/cluster_data2/"
setwd(base_dir2)
con1 <- read.delim("parasp/hg19_jlRNAseq/results/03_featurecount/G1_featureCounts.txt",
                   header = F,sep="\t", quote = "")
con2 <- read.delim("parasp/hg19_jlRNAseq/results/03_featurecount/G2_featureCounts.txt",
                   header = F,sep="\t", quote = "")
PSF1 <- read.delim("parasp/hg19_jlRNAseq/results/03_featurecount/F1_featureCounts.txt",
                   header = F,sep="\t", quote = "")
PSF2 <- read.delim("parasp/hg19_jlRNAseq/results/03_featurecount/F2_featureCounts.txt",
                   header = F,sep="\t", quote = "")
pre_exp<-list(con1,con2,PSF1,PSF2) %>% reduce(left_join,by="V1")
# filter 2 out of 4 samples > 5 and sum() >20
exp <- pre_exp[apply(pre_exp[,c(2:5)]>10,1,sum)>2 & apply(pre_exp[,c(2:5)],1,sum)>40,]
libsizes <- apply(exp[,c(2:5)],2,sum)
normfactor <- libsizes/libsizes[1]
exp_norm <- apply(exp[,c(2:5)],2,function(x) {x/normfactor})
exp$e_name <-substr(exp$V1,1,15)
# head(exp)
exp$log2fc<-log2((exp_norm[,3]+exp_norm[,4])/(exp_norm[,1]+exp_norm[,2]))
IR_pre <-IR_df_dup[IR_df_dup$IR_RPK > 0,]
colnames(exp) <-c("ID","con1","con2","kdPSF1","kdPSF2","e_name","log2fc")
IR_plot2 <- merge(IR_pre,exp,by.x = "geneID",by.y="e_name")
p <-ggplot(IR_plot2,(aes(x=IR_RPK,y=log2fc))) +
  geom_point(size=0.5)+
  geom_smooth(aes(x=IR_RPK,y=log2fc, colour="#008B8B"))+
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major=element_line(colour=NA),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=10),
        axis.text.y = element_text(face="bold",  color="black", size=10),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(y="Related genes log2Foldchange", x= "IR region interact intensity (RPM)")
pdf("Gride-seq/results/06_IR_interact/ExpresLevelVS_iteractIntensity.pdf")
print(p)
dev.off()
