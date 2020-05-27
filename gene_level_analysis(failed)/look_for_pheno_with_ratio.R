library(ggplot2)

pheno<-readRDS("TCGABasicPhenotypeData.rds")
RNA_structure_ratio<-read.csv("RNA_structure_ratio.csv")

## look for obvious phenotype correlations

RNA_structure_pheno<-merge(pheno, RNA_structure_ratio,
                           by.x = "Sample ID", by.y = "sample")

# pair normal and tumor samples where available
sum(unlist(lapply(RNA_structure_ratio$sample, function(x){substr(x, 14, 15) == "11"})))
# only 30 normal samples; will ignore these

give.n <- function(x){
  return(c(y = max(x)*1.1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
ggplot(RNA_structure_pheno, aes(x=Sex, y=ratio_HSU_PSU)) +
  geom_violin(scale = "count")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)
ggplot(RNA_structure_pheno, aes(x=`Fraction Genome Altered`, y=ratio_HSU_PSU)) +
  geom_point()+
  geom_smooth(method='lm', formula= y~x)
ggplot(RNA_structure_pheno, aes(x=log(`Mutation Count`), y=ratio_HSU_PSU)) +
  geom_point(alpha = .2)+
  geom_smooth(method='lm', formula= y~x)
summary(lm(RNA_structure_pheno$ratio_HSU_PSU ~ log(RNA_structure_pheno$`Mutation Count`)))
ggplot(RNA_structure_pheno, aes(x=log(`Mutation Count`), y=G3BP1)) +
  geom_point(alpha = .2)+
  geom_smooth(method='lm', formula= y~x)
summary(lm(RNA_structure_pheno$G3BP1 ~ log(RNA_structure_pheno$`Mutation Count`)))
ggplot(RNA_structure_pheno, aes(x=`Diagnosis Age`, y=ratio_HSU_PSU)) +
  geom_point()+
  geom_smooth(method='lm', formula= y~x)
ggplot(RNA_structure_pheno, aes(x=reorder(`Cancer Type`, ratio_HSU_PSU, FUN = median), y=ratio_HSU_PSU)) +
  geom_boxplot()+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(RNA_structure_pheno, aes(x=`Sample Type`, y=ratio_HSU_PSU)) +
  geom_violin(scale = "count")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)
ggplot(RNA_structure_pheno, aes(x=`Race Category`, y=ratio_HSU_PSU)) +
  geom_violin(scale = "count")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



