library(tidyverse)

# correlate G3BP1 with all other genes

cancer<-readRDS("RNAseqPanCancerTCGA.rds")

G3BP1_corr<-as.data.frame(matrix(NA, nrow = nrow(cancer), ncol = 4))
colnames(G3BP1_corr)<-c("gene", "correlation", "covariance", "variance")
G3BP1_corr$gene<-rownames(cancer)
for(gene in rownames(cancer)){
  cor<-cor(unlist(cancer["G3BP1",]), unlist(cancer[gene,]), method = "spearman")
  cov<-cov(unlist(cancer["G3BP1",]), unlist(cancer[gene,]), method = "spearman")
  var<-var(unlist(cancer[gene,]))
  G3BP1_corr[which(G3BP1_corr$gene == gene),2:4]<-c(cor, cov, var)
}

write.csv(G3BP1_corr, "G3BP1_correlation.csv", row.names = F)

## log2 transform
G3BP1_corr_log2<-as.data.frame(matrix(NA, nrow = nrow(cancer), ncol = 4))
colnames(G3BP1_corr_log2)<-c("gene", "correlation", "covariance", "variance")
G3BP1_corr_log2$gene<-rownames(cancer)
for(gene in rownames(cancer)){
  cor<-cor(log2(unlist(cancer["G3BP1",])+1), log2(unlist(cancer[gene,])+1), method = "spearman")
  cov<-cov(log2(unlist(cancer["G3BP1",])+1), log2(unlist(cancer[gene,])+1), method = "spearman")
  var<-var(unlist(cancer[gene,]))
  G3BP1_corr_log2[which(G3BP1_corr_log2$gene == gene),2:4]<-c(cor, cov, var)
}

write.csv(G3BP1_corr_log2, "G3BP1_correlation_log2.csv", row.names = F)