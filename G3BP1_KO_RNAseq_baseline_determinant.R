library(dplyr)
library(tidyr)
library(ggplot2)

UTRs<-read.csv("3'_UTRs.csv", stringsAsFactors = F)
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
# still no obvious trends

RNA_structure_ratio$cell<-substr(RNA_structure_ratio$sample, 13, 14)
RNA_structure_ratio$timepoint<-substr(RNA_structure_ratio$sample, 15,
                                      nchar(RNA_structure_ratio$sample)-3)
ggplot(RNA_structure_ratio[1:12,], aes(x=as.numeric(timepoint), y=ratio_HSU_PSU))+
  geom_point(aes(color = cell), size = 4)
ggplot(RNA_structure_ratio[13:16,], aes(x=timepoint, y=ratio_HSU_PSU))+
  geom_point(size = 4)



### return to -.25 dG/nt cut-off and look at genes that have all HSU or PSU
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.25)
temp<-HSU$gene_id
PSU<-all_data %>% filter(Delta.G.per.Nucleotide >= -.25)
HSU<-HSU[!(HSU$gene_id %in% PSU$gene_id),]
PSU<-PSU[!(PSU$gene_id %in% temp),]
rm(temp)

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


### try as above (all isoforms HSU or PSU) with -.2/-.3 split
HSU<-all_data %>% filter(Delta.G.per.Nucleotide < -.3)
temp<-HSU$gene_id
PSU<-all_data %>% filter(Delta.G.per.Nucleotide > -.2)
HSU<-HSU[!(HSU$gene_id %in% PSU$gene_id),]
PSU<-PSU[!(PSU$gene_id %in% temp),]
rm(temp)

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
