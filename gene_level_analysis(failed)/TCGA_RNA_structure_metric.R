library(tidyverse)

cancer<-readRDS("RNAseqPanCancerTCGA.rds")
pheno<-readRDS("TCGABasicPhenotypeData.rds")

UTRs<-read.csv("3'_UTRs.csv", stringsAsFactors = F)[,c("Gene.Name", "Delta.G.per.Nucleotide")]
all_data<-merge(UTRs, cancer, by.x = "Gene.Name", by.y = 0)

### assign min dG/nt

HSU_gene_min<-all_data %>% group_by(Gene.Name) %>%
  filter(max(Delta.G.per.Nucleotide) < -.25) %>%
  slice(which.min(Delta.G.per.Nucleotide))
PSU_gene_min<-all_data %>% group_by(Gene.Name) %>%
  filter(min(Delta.G.per.Nucleotide) > -.25) %>%
  slice(which.min(Delta.G.per.Nucleotide))

HSU_gene_min[,3:ncol(all_data)]<-HSU_gene_min[,3:ncol(all_data)] * HSU_gene_min$Delta.G.per.Nucleotide
PSU_gene_min[,3:ncol(all_data)]<-PSU_gene_min[,3:ncol(all_data)] * PSU_gene_min$Delta.G.per.Nucleotide

RNA_structure_ratio_min<-as.data.frame(matrix(NA, nrow = ncol(all_data)-2, ncol = 4))
colnames(RNA_structure_ratio_min)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio_min$sample<-colnames(all_data[3:ncol(all_data)])
RNA_structure_ratio_min$sum_PSU<-apply(PSU_gene_min[,3:ncol(all_data)], 2, sum)
RNA_structure_ratio_min$sum_HSU<-apply(HSU_gene_min[,3:ncol(all_data)], 2, sum)
RNA_structure_ratio_min$ratio_HSU_PSU<-RNA_structure_ratio_min$sum_HSU/RNA_structure_ratio_min$sum_PSU
RNA_structure_ratio_min$G3BP1<-all_data[which(all_data$Gene.Name == "G3BP1")[1], 3:ncol(all_data)] %>% unlist()
RNA_structure_ratio_min$METTL3<-all_data[which(all_data$Gene.Name == "METTL3")[1], 3:ncol(all_data)] %>% unlist()

ggplot(RNA_structure_ratio_min, aes(x=log2(METTL3 + 1), y=ratio_HSU_PSU)) +
  geom_point(alpha = .2)+
  geom_smooth(method='lm', formula= y~x)

summary(lm(RNA_structure_ratio_min$ratio_HSU_PSU ~ log2(RNA_structure_ratio_min$METTL3)))$coefficients

ggplot(RNA_structure_ratio_min, aes(x=log2(G3BP1 + 1), y=ratio_HSU_PSU)) +
  geom_point(alpha = .2)+
  geom_smooth(method='lm', formula= y~x)

summary(lm(RNA_structure_ratio_min$ratio_HSU_PSU ~ log2(RNA_structure_ratio_min$G3BP1)))$coefficients

ggplot(RNA_structure_ratio_min, aes(x=log2(METTL3 + 1), y=log2(G3BP1+1))) +
   geom_point(alpha = .2)+
   geom_smooth(method='lm', formula= y~x)

summary(lm(log2(RNA_structure_ratio_min$G3BP1 +1) ~ log2(RNA_structure_ratio_min$METTL3 +1)))$coefficients

summary(lm(RNA_structure_ratio_min$ratio_HSU_PSU ~ log2(RNA_structure_ratio_min$G3BP1 +1) + log2(RNA_structure_ratio_min$METTL3 +1)))$coefficients

#### assign median dG/nt (this is what is ultimately used; not sig different from min dG/nt)

HSU_gene<-all_data %>% group_by(Gene.Name) %>%
  filter(max(Delta.G.per.Nucleotide) < -.25) %>%
  mutate(Delta.G.per.Nucleotide=median(Delta.G.per.Nucleotide)) %>%
  distinct()
PSU_gene<-all_data %>% group_by(Gene.Name) %>%
  filter(min(Delta.G.per.Nucleotide) > -.25) %>%
  mutate(Delta.G.per.Nucleotide=median(Delta.G.per.Nucleotide)) %>%
  distinct()

HSU_gene[,3:ncol(all_data)]<-HSU_gene[,3:ncol(all_data)] * HSU_gene$Delta.G.per.Nucleotide
PSU_gene[,3:ncol(all_data)]<-PSU_gene[,3:ncol(all_data)] * PSU_gene$Delta.G.per.Nucleotide

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = ncol(all_data)-2, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:ncol(all_data)])
RNA_structure_ratio$sum_PSU<-apply(PSU_gene[,3:ncol(all_data)], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU_gene[,3:ncol(all_data)], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU
RNA_structure_ratio$G3BP1<-all_data[which(all_data$Gene.Name == "G3BP1")[1], 3:ncol(all_data)] %>% unlist()

ggplot(RNA_structure_ratio, aes(x=log2(G3BP1 + 1), y=ratio_HSU_PSU)) +
  geom_point(alpha = .2)+
  geom_smooth(method='lm', formula= y~x)

summary(lm(RNA_structure_ratio$ratio_HSU_PSU ~ log2(RNA_structure_ratio$G3BP1)))$coefficients

write.csv(RNA_structure_ratio, "RNA_structure_ratio.csv", row.names = F)

## systematically look for proteins most associated with this trend
gene_ratio_corr<-as.data.frame(matrix(NA, ncol = 7, nrow = length(unique(all_data$Gene.Name))))
colnames(gene_ratio_corr)<-c("gene", "Estimate", "Std.Error", "t.value", "p.value", "corr", "var")
gene_ratio_corr$gene<-unique(all_data$Gene.Name)
all_HSU_PSU_genes<-c(HSU_gene$Gene.Name, PSU_gene$Gene.Name)

for(Gene.Name in gene_ratio_corr$gene){
  if(Gene.Name %in% all_HSU_PSU_genes){
    gene_ratio_corr[which(gene_ratio_corr$gene == Gene.Name), 7]<-var(unlist(all_data[which(all_data$Gene.Name == Gene.Name)[1], 3:ncol(all_data)]))
    if(gene_ratio_corr[which(gene_ratio_corr$gene == Gene.Name), 7] != 0){
      gene_ratio_corr[which(gene_ratio_corr$gene == Gene.Name), 2:6]<- c(summary(lm(RNA_structure_ratio$ratio_HSU_PSU ~ 
                                                                                    log2(unlist(all_data[which(all_data$Gene.Name == Gene.Name)[1],
                                                                                                         3:ncol(all_data)])+1)))$coefficients[2,1:4],
                                                                         cor(RNA_structure_ratio$ratio_HSU_PSU,
                                                                             log2(unlist(all_data[which(all_data$Gene.Name == Gene.Name)[1],
                                                                                                  3:ncol(all_data)])+1), method = "spearman"))
    }
  }
}
gene_ratio_corr<-gene_ratio_corr[!is.na(gene_ratio_corr$Estimate),]
gene_ratio_corr$abs_est<-abs(gene_ratio_corr$Estimate)

write.csv(gene_ratio_corr, "HSU_PSU_ratio.csv", row.names = F)

top_95_percent_var<-top_frac(gene_ratio_corr, .95, var)

top_5_percent_estimate<-top_frac(top_95_percent_var, .05, abs_est)
top_5_percent_estimate_neg<-top_frac(top_95_percent_var, -.05, Estimate)

write(top_5_percent_estimate$gene, "top_5_signal.txt")
write(top_5_percent_estimate_neg$gene, "top_5_signal_neg.txt")

