args<-commandArgs(T)
library(ggplot2)
library(ggthemes)
r1 <- read.table(args[1],sep="\t",header =F ,comment.char = "#")$V7
r2 <- read.table(args[2],sep="\t",header =F ,comment.char = "#")$V7
dat <-data.frame(rep1=r1,rep2=r2)
dat_exp <-dat[dat$rep1 >0 & dat$rep2 >0,]
rep_lm <- lm(rep1 ~ rep2, data=dat_exp)
R2 <- summary(rep_lm)$r.squared 
# R2=0.9934815
pdf(args[3])
ggplot(dat, aes(x=log2(rep1), y=log2(rep2))) +
  geom_point() + labs(title = "",y="log2(rep2 num in every bin)",x="log2(rep1 num in every bin)")+
    theme_bw(base_size = 12, base_family = "Times")
dev.off()