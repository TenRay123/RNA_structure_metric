library(dplyr)

UTRs<-read.csv("3'_UTRs.csv", stringsAsFactors = F)
isoform_info<-read.delim("isoforms_tpm_all_samples_norm.txt", stringsAsFactors = F)
UTRs$Accession.Number<-gsub("\\..*","",UTRs$Accession.Number)

all_data<-merge(isoform_info, UTRs, by.x = "transcript_id", by.y = "Accession.Number")
# some gene names don't match up, but a few spot checks indicate this is due to
# changes in conventional naming rather than true mis-matches

