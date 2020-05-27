library(dplyr)
library(tidyr)
library(ggplot2)

UTRs<-read.csv("3'_UTRs.csv", stringsAsFactors = F)
gene_info<-read.delim("genes_tpm_all_samples_norm.txt", stringsAsFactors = F)
isoform_info<-read.delim("isoforms_tpm_all_samples_norm.txt", stringsAsFactors = F)
UTRs$Accession.Number<-gsub("\\..*","",UTRs$Accession.Number)

all_data<-merge(isoform_info, UTRs, by.x = "transcript_id", by.y = "Accession.Number")
# some gene names don't match up, but a few spot checks indicate this is due to
# changes in conventional naming rather than true mis-matches

#### all isoform data split along -.25 dG/nt
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.25)
PSU<-all_data %>% filter(Delta.G.per.Nucleotide >= -.25)

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU

ggplot(RNA_structure_ratio, aes(x=sample, y=ratio_HSU_PSU))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# no obvious trends

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)

#### re-do so isoforms are split < -.3 and > -.2 dG/nt
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.3)
PSU<-all_data %>% filter(Delta.G.per.Nucleotide > -.2)

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU

ggplot(RNA_structure_ratio, aes(x=sample, y=ratio_HSU_PSU))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)


### return to -.25 dG/nt cut-off and look at genes that have all HSU or PSU
HSU_gene<-all_data %>% group_by(Gene.Name) %>%
  filter(max(Delta.G.per.Nucleotide) < -.25)
PSU_gene<-all_data %>% group_by(Gene.Name) %>%
  filter(min(Delta.G.per.Nucleotide) > -.25)

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU_gene[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU_gene[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU
RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)


### try as above (all isoforms HSU or PSU) with -.2/-.3 split
HSU_gene<-all_data %>% group_by(Gene.Name) %>%
  filter(max(Delta.G.per.Nucleotide) < -.3)
PSU_gene<-all_data %>% group_by(Gene.Name) %>%
  filter(min(Delta.G.per.Nucleotide) > -.2)

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU_gene[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU_gene[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU
RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)


### look at this from gene-level rather than isoform-level
all_gene_data<-merge(UTRs, gene_info, by.x = "Gene.Name", by.y = "gene_id")
HSU_gene<-all_gene_data %>% group_by(Gene.Name) %>%
  filter(max(Delta.G.per.Nucleotide) < -.25)
PSU_gene<-all_gene_data %>% group_by(Gene.Name) %>%
  filter(min(Delta.G.per.Nucleotide) > -.25)

HSU_gene<-HSU_gene[,c(1,14:29)] %>% unique()
PSU_gene<-PSU_gene[,c(1,14:29)] %>% unique()

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(HSU_gene[2:17])
RNA_structure_ratio$sum_PSU<-apply(PSU_gene[,2:17], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU_gene[,2:17], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU
RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)


### mutliply UTR structure by TPM if struc < -.25 at isoform level
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.25)
PSU<-all_data %>% filter(Delta.G.per.Nucleotide >= -.25)

HSU[,3:18]<-HSU[,3:18] * HSU$Delta.G.per.Nucleotide
PSU[,3:18]<-PSU[,3:18] * PSU$Delta.G.per.Nucleotide

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)
### best separation so far ###

### mutliply UTR structure by TPM along -.2 and -.3 lines at isoform level
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.3)
PSU<-all_data %>% filter(Delta.G.per.Nucleotide > -.2)

HSU[,3:18]<-HSU[,3:18] * HSU$Delta.G.per.Nucleotide
PSU[,3:18]<-PSU[,3:18] * PSU$Delta.G.per.Nucleotide

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)

### mutliply UTR structure by TPM if struc < -.25 at isoform level
### but only multiply HSU by struc; PSUs all by -.25
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.25)
PSU<-all_data %>% filter(Delta.G.per.Nucleotide >= -.25)

HSU[,3:18]<-HSU[,3:18] * HSU$Delta.G.per.Nucleotide
PSU[,3:18]<-PSU[,3:18] * -.25

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 4))
colnames(RNA_structure_ratio)<-c("sample", "sum_HSU", "sum_PSU", "ratio_HSU_PSU")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_PSU<-apply(PSU[,3:18], 2, sum)
RNA_structure_ratio$sum_HSU<-apply(HSU[,3:18], 2, sum)
RNA_structure_ratio$ratio_HSU_PSU<-RNA_structure_ratio$sum_HSU/RNA_structure_ratio$sum_PSU

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)

### just plot sum of TPM*struc
plot_data<-all_data[,3:18]*all_data$Delta.G.per.Nucleotide

RNA_structure_ratio<-as.data.frame(matrix(NA, nrow = 16, ncol = 2))
colnames(RNA_structure_ratio)<-c("sample", "sum_structure")

RNA_structure_ratio$sample<-colnames(all_data[3:18])
RNA_structure_ratio$sum_structure<-apply(plot_data, 2, sum)

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=abs(sum_structure)))+
  geom_point(aes(color = cell), size = 4)+
  scale_y_log10()
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=sum_structure))+
  geom_point(size = 4)
