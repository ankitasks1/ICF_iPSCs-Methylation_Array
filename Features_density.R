#Features to consider
#Imprinted DMRs /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed
#PCDH PCDHcluster.txt
#TNXB TNXB_hypo2000.txt
setwd("/media/ankitv/Archivio2/ankit/Array21/Controls")
awk '{print $1"\t"$2"\t"$3"\t"$6}' /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.gtf > gencode.v29.annotation.gene_re.gtf
grep Start CGIs_hg38.txt -v | awk '{print $2,$3,$4,$5,$6}' OFS="\t" | sort -k1,1 -k2,2n  > CGIs_hg38_re.txt

#bedtools makewindows -b gencode.v29.annotation.gene_re.gtf -n 100 -i srcwinnum | sort -k1,1 -k2,2n > Gene_hg38_100.tx
#CGI
bedtools intersect -wa -wb -a Monk_human_ICR_hg38.bed -b CGIs_hg38_re.txt > Monk_human_ICR_CGIs_hg38.txt
bedtools intersect -wa -wb -a PCDHcluster.txt -b CGIs_hg38_re.txt > PCDHcluster_CGIs_hg38.txt
bedtools intersect -wa -wb -a TNXB_hypo2000.txt -b CGIs_hg38_re.txt > TNXB_hypo2000_CGIs_hg38.txt

#CG 
cp /media/ankitv/Archivio2/ankit/motif_analysis/marcos/motif_trovati/human/motif[CG].human_hg38.bed CG_hg38_re.txt
bedtools intersect -wa -wb -a Monk_human_ICR_hg38.bed -b CG_hg38_re.txt > Monk_human_ICR_CG_hg38.txt
bedtools intersect -wa -wb -a PCDHcluster.txt -b CG_hg38_re.txt > PCDHcluster_CG_hg38.txt
bedtools intersect -wa -wb -a TNXB_hypo2000.txt -b CG_hg38_re.txt > TNXB_hypo2000_CG_hg38.txt

Monk_human_ICR_CG_hg38  <- read.table("Monk_human_ICR_CG_hg38.txt")
head(Monk_human_ICR_CG_hg38)
detach("package:dplyr")
library(plyr)
Monk_human_ICR_CG_hg38count  <- count(Monk_human_ICR_CG_hg38, "V4")
head(Monk_human_ICR_CG_hg38count)
Monk_human_ICR_CG_hg38count <- Monk_human_ICR_CG_hg38count[order(-Monk_human_ICR_CG_hg38count$freq),]
Monk_human_ICR_CG_hg38count_merge <- merge(Monk_human_ICR_CG_hg38count,Monk_human_ICR_CG_hg38, "V4")
head(Monk_human_ICR_CG_hg38count_merge)
Monk_human_ICR_CG_hg38count_merge <- Monk_human_ICR_CG_hg38count_merge[order(-Monk_human_ICR_CG_hg38count_merge$freq),]
Monk_human_ICR_CG_hg38count_merge["GeneSpan"] <- Monk_human_ICR_CG_hg38count_merge$V3 -Monk_human_ICR_CG_hg38count_merge$V2
Monk_human_ICR_CG_hg38count_merge["Density"] <- Monk_human_ICR_CG_hg38count_merge$freq * 100/Monk_human_ICR_CG_hg38count_merge$GeneSpan
Monk_human_ICR_CG_hg38count_merge_uniq <- Monk_human_ICR_CG_hg38count_merge[,c(1,2,14)]
Monk_human_ICR_CG_hg38count_merge_uniq <- data.frame(unique(paste0(Monk_human_ICR_CG_hg38count_merge_uniq$V4,"%",Monk_human_ICR_CG_hg38count_merge_uniq$freq,"%",Monk_human_ICR_CG_hg38count_merge_uniq$Density)))
colnames(Monk_human_ICR_CG_hg38count_merge_uniq) <- "col"
head(Monk_human_ICR_CG_hg38count_merge_uniq)
Monk_human_ICR_CG_hg38count_merge_uniq <- cSplit(Monk_human_ICR_CG_hg38count_merge_uniq,"col","%")
Monk_human_ICR_CG_hg38count_merge_uniq <- data.frame(Monk_human_ICR_CG_hg38count_merge_uniq)
colnames(Monk_human_ICR_CG_hg38count_merge_uniq) <- c("Feature","Freq","Density")
Monk_human_ICR_CG_hg38count_merge_uniq <- Monk_human_ICR_CG_hg38count_merge_uniq[order(-Monk_human_ICR_CG_hg38count_merge_uniq$Density),]
Monk_human_ICR_CG_hg38count_merge_uniq["Type"] <- "Known_ICR"

PCDHcluster_CG_hg38  <- read.table("PCDHcluster_CG_hg38.txt")
head(PCDHcluster_CG_hg38)
library(plyr)
PCDHcluster_CG_hg38count  <- count(PCDHcluster_CG_hg38, "V4")
head(PCDHcluster_CG_hg38count)
PCDHcluster_CG_hg38count <- PCDHcluster_CG_hg38count[order(-PCDHcluster_CG_hg38count$freq),]
PCDHcluster_CG_hg38count_merge <- merge(PCDHcluster_CG_hg38count,PCDHcluster_CG_hg38, "V4")
head(PCDHcluster_CG_hg38count_merge)
PCDHcluster_CG_hg38count_merge <- PCDHcluster_CG_hg38count_merge[order(-PCDHcluster_CG_hg38count_merge$freq),]
PCDHcluster_CG_hg38count_merge["GeneSpan"] <- PCDHcluster_CG_hg38count_merge$V3 -PCDHcluster_CG_hg38count_merge$V2
PCDHcluster_CG_hg38count_merge["Density"] <- PCDHcluster_CG_hg38count_merge$freq * 100/PCDHcluster_CG_hg38count_merge$GeneSpan
PCDHcluster_CG_hg38count_merge_uniq <- PCDHcluster_CG_hg38count_merge[,c(1,4,5,2,12)]
PCDHcluster_CG_hg38count_merge_uniq <- PCDHcluster_CG_hg38count_merge_uniq[unique(PCDHcluster_CG_hg38count_merge_uniq$V4),]
PCDHcluster_CG_hg38count_merge_uniq <- PCDHcluster_CG_hg38count_merge_uniq[order(-PCDHcluster_CG_hg38count_merge_uniq$Density),]
PCDHcluster_CG_hg38count_merge_uniq["Feature"] <- paste0("PCDH_gene_cluster","_",
                                                         PCDHcluster_CG_hg38count_merge_uniq$V4,
                                                         "%",
                                                         PCDHcluster_CG_hg38count_merge_uniq$V2,
                                                         "%",
                                                         PCDHcluster_CG_hg38count_merge_uniq$V3)

PCDHcluster_CG_hg38count_merge_uniq <- PCDHcluster_CG_hg38count_merge_uniq[,c(6,4,5)]
PCDHcluster_CG_hg38count_merge_uniq["Type"] <- "PCDH_gene_cluster"
colnames(PCDHcluster_CG_hg38count_merge_uniq) <- c("Feature", "Freq","Density","Type")


TNXB_hypo2000_CG_hg38  <- read.table("TNXB_hypo2000_CG_hg38.txt")
head(TNXB_hypo2000_CG_hg38)
library(plyr)
TNXB_hypo2000_CG_hg38count  <- count(TNXB_hypo2000_CG_hg38, "V4")
head(TNXB_hypo2000_CG_hg38count)
TNXB_hypo2000_CG_hg38count <- TNXB_hypo2000_CG_hg38count[order(-TNXB_hypo2000_CG_hg38count$freq),]
TNXB_hypo2000_CG_hg38count_merge <- merge(TNXB_hypo2000_CG_hg38count,TNXB_hypo2000_CG_hg38, "V4")
head(TNXB_hypo2000_CG_hg38count_merge)
TNXB_hypo2000_CG_hg38count_merge <- TNXB_hypo2000_CG_hg38count_merge[order(-TNXB_hypo2000_CG_hg38count_merge$freq),]
TNXB_hypo2000_CG_hg38count_merge["GeneSpan"] <- TNXB_hypo2000_CG_hg38count_merge$V3 -TNXB_hypo2000_CG_hg38count_merge$V2
TNXB_hypo2000_CG_hg38count_merge["Density"] <- TNXB_hypo2000_CG_hg38count_merge$freq * 100/TNXB_hypo2000_CG_hg38count_merge$GeneSpan
TNXB_hypo2000_CG_hg38count_merge_uniq <- TNXB_hypo2000_CG_hg38count_merge[,c(1,4,5,2,12)]
TNXB_hypo2000_CG_hg38count_merge_uniq <- TNXB_hypo2000_CG_hg38count_merge_uniq[unique(TNXB_hypo2000_CG_hg38count_merge_uniq$V4),]
TNXB_hypo2000_CG_hg38count_merge_uniq <- TNXB_hypo2000_CG_hg38count_merge_uniq[order(-TNXB_hypo2000_CG_hg38count_merge_uniq$Density),]
TNXB_hypo2000_CG_hg38count_merge_uniq["Feature"] <- paste0("TNXB_hypo2000","_",
                                                           TNXB_hypo2000_CG_hg38count_merge_uniq$V4,
                                                           "%",
                                                           TNXB_hypo2000_CG_hg38count_merge_uniq$V2,
                                                           "%",
                                                           TNXB_hypo2000_CG_hg38count_merge_uniq$V3)

TNXB_hypo2000_CG_hg38count_merge_uniq <- TNXB_hypo2000_CG_hg38count_merge_uniq[,c(6,4,5)]
TNXB_hypo2000_CG_hg38count_merge_uniq["Type"] <- "TNXB_hpyo2000"
colnames(TNXB_hypo2000_CG_hg38count_merge_uniq) <- c("Feature", "Freq","Density","Type")

Monk_human_ICR_PCDH_TNXB_CG <- rbind.data.frame(Monk_human_ICR_CG_hg38count_merge_uniq, PCDHcluster_CG_hg38count_merge_uniq, TNXB_hypo2000_CG_hg38count_merge_uniq)
Monk_human_ICR_PCDH_TNXB_CG <- Monk_human_ICR_PCDH_TNXB_CG[order(-Monk_human_ICR_PCDH_TNXB_CG$Density),]
write.table(Monk_human_ICR_PCDH_TNXB_CG, "Monk_human_ICR_PCDH_TNXB_CG.txt", sep = "\t", append = F, quote = F, row.names = F)
ggplot(Monk_human_ICR_PCDH_TNXB_CG) + 
  geom_bar(aes(x=reorder(Feature, Density), y=Density, col=Feature, fill= Feature), stat = "identity", position = "dodge")+
  coord_flip()+theme_bw()+
  scale_fill_manual(values = c("#5573e8","#C7CEEA","#1d46ea","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#1d46ea","#1d46ea","#C7CEEA","#C7CEEA","#1d46ea","#C7CEEA","#C7CEEA","#C7CEEA","#5573e8","#1d46ea","#C7CEEA","#5573e8","#C7CEEA","#C7CEEA","#C7CEEA","#1d46ea","#C7CEEA","#C7CEEA","#5573e8","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#C7CEEA","#1d46ea","#1d46ea","#C7CEEA","darkgreen","orange"))+
  scale_color_manual(values = c(rep("white",50),"white","white"))
  
ggsave("Monk_human_ICR_PCDH_TNXB_CG.png", width=12, height=10, units="in", dpi=96)





#---------------------------  CGIs  ----------------------------#

Monk_human_ICR_CGIs_hg38  <- read.table("Monk_human_ICR_CGIs_hg38.txt")
head(Monk_human_ICR_CGIs_hg38)
library(plyr)
Monk_human_ICR_CGIs_hg38count  <- count(Monk_human_ICR_CGIs_hg38, "V4")
head(Monk_human_ICR_CGIs_hg38count)
Monk_human_ICR_CGIs_hg38count <- Monk_human_ICR_CGIs_hg38count[order(-Monk_human_ICR_CGIs_hg38count$freq),]
Monk_human_ICR_CGIs_hg38count_merge <- merge(Monk_human_ICR_CGIs_hg38count,Monk_human_ICR_CGIs_hg38, "V4")
head(Monk_human_ICR_CGIs_hg38count_merge)
Monk_human_ICR_CGIs_hg38count_merge <- Monk_human_ICR_CGIs_hg38count_merge[order(-Monk_human_ICR_CGIs_hg38count_merge$freq),]
Monk_human_ICR_CGIs_hg38count_merge["GeneSpan"] <- Monk_human_ICR_CGIs_hg38count_merge$V3 -Monk_human_ICR_CGIs_hg38count_merge$V2
Monk_human_ICR_CGIs_hg38count_merge["Density"] <- Monk_human_ICR_CGIs_hg38count_merge$freq * 100/Monk_human_ICR_CGIs_hg38count_merge$GeneSpan
Monk_human_ICR_CGIs_hg38count_merge_uniq <- Monk_human_ICR_CGIs_hg38count_merge[,c(1,2,13)]
Monk_human_ICR_CGIs_hg38count_merge_uniq <- data.frame(unique(paste0(Monk_human_ICR_CGIs_hg38count_merge_uniq$V4,"%",Monk_human_ICR_CGIs_hg38count_merge_uniq$freq,"%",Monk_human_ICR_CGIs_hg38count_merge_uniq$Density)))
colnames(Monk_human_ICR_CGIs_hg38count_merge_uniq) <- "col"
head(Monk_human_ICR_CGIs_hg38count_merge_uniq)
Monk_human_ICR_CGIs_hg38count_merge_uniq <- cSplit(Monk_human_ICR_CGIs_hg38count_merge_uniq,"col","%")
Monk_human_ICR_CGIs_hg38count_merge_uniq <- data.frame(Monk_human_ICR_CGIs_hg38count_merge_uniq)
colnames(Monk_human_ICR_CGIs_hg38count_merge_uniq) <- c("Feature","Freq","Density")
Monk_human_ICR_CGIs_hg38count_merge_uniq <- Monk_human_ICR_CGIs_hg38count_merge_uniq[order(-Monk_human_ICR_CGIs_hg38count_merge_uniq$Density),]
Monk_human_ICR_CGIs_hg38count_merge_uniq["Type"] <- "Known_ICR"

PCDHcluster_CGIs_hg38  <- read.table("PCDHcluster_CGIs_hg38.txt")
head(PCDHcluster_CGIs_hg38)
library(plyr)
PCDHcluster_CGIs_hg38count  <- count(PCDHcluster_CGIs_hg38, "V4")
head(PCDHcluster_CGIs_hg38count)
PCDHcluster_CGIs_hg38count <- PCDHcluster_CGIs_hg38count[order(-PCDHcluster_CGIs_hg38count$freq),]
PCDHcluster_CGIs_hg38count_merge <- merge(PCDHcluster_CGIs_hg38count,PCDHcluster_CGIs_hg38, "V4")
head(PCDHcluster_CGIs_hg38count_merge)
PCDHcluster_CGIs_hg38count_merge <- PCDHcluster_CGIs_hg38count_merge[order(-PCDHcluster_CGIs_hg38count_merge$freq),]
PCDHcluster_CGIs_hg38count_merge["GeneSpan"] <- PCDHcluster_CGIs_hg38count_merge$V3 -PCDHcluster_CGIs_hg38count_merge$V2
PCDHcluster_CGIs_hg38count_merge["Density"] <- PCDHcluster_CGIs_hg38count_merge$freq * 100/PCDHcluster_CGIs_hg38count_merge$GeneSpan
PCDHcluster_CGIs_hg38count_merge_uniq <- PCDHcluster_CGIs_hg38count_merge[,c(1,4,5,2,11)]
PCDHcluster_CGIs_hg38count_merge_uniq <- PCDHcluster_CGIs_hg38count_merge_uniq[unique(PCDHcluster_CGIs_hg38count_merge_uniq$V4),]
PCDHcluster_CGIs_hg38count_merge_uniq <- PCDHcluster_CGIs_hg38count_merge_uniq[order(-PCDHcluster_CGIs_hg38count_merge_uniq$Density),]
PCDHcluster_CGIs_hg38count_merge_uniq["Feature"] <- paste0("PCDH_gene_cluster","_",
                                                           PCDHcluster_CGIs_hg38count_merge_uniq$V4,
                                                           "%",
                                                           PCDHcluster_CGIs_hg38count_merge_uniq$V2,
                                                           "%",
                                                           PCDHcluster_CGIs_hg38count_merge_uniq$V3)

PCDHcluster_CGIs_hg38count_merge_uniq <- PCDHcluster_CGIs_hg38count_merge_uniq[,c(6,4,5)]
PCDHcluster_CGIs_hg38count_merge_uniq["Type"] <- "PCDH_gene_cluster"
colnames(PCDHcluster_CGIs_hg38count_merge_uniq) <- c("Feature", "Freq","Density","Type")


TNXB_hypo2000_CGIs_hg38  <- read.table("TNXB_hypo2000_CGIs_hg38.txt")
head(TNXB_hypo2000_CGIs_hg38)
library(plyr)
TNXB_hypo2000_CGIs_hg38count  <- count(TNXB_hypo2000_CGIs_hg38, "V4")
head(TNXB_hypo2000_CGIs_hg38count)
TNXB_hypo2000_CGIs_hg38count <- TNXB_hypo2000_CGIs_hg38count[order(-TNXB_hypo2000_CGIs_hg38count$freq),]
TNXB_hypo2000_CGIs_hg38count_merge <- merge(TNXB_hypo2000_CGIs_hg38count,TNXB_hypo2000_CGIs_hg38, "V4")
head(TNXB_hypo2000_CGIs_hg38count_merge)
TNXB_hypo2000_CGIs_hg38count_merge <- TNXB_hypo2000_CGIs_hg38count_merge[order(-TNXB_hypo2000_CGIs_hg38count_merge$freq),]
TNXB_hypo2000_CGIs_hg38count_merge["GeneSpan"] <- TNXB_hypo2000_CGIs_hg38count_merge$V3 -TNXB_hypo2000_CGIs_hg38count_merge$V2
TNXB_hypo2000_CGIs_hg38count_merge["Density"] <- TNXB_hypo2000_CGIs_hg38count_merge$freq * 100/TNXB_hypo2000_CGIs_hg38count_merge$GeneSpan
TNXB_hypo2000_CGIs_hg38count_merge_uniq <- TNXB_hypo2000_CGIs_hg38count_merge[,c(1,4,5,2,11)]
TNXB_hypo2000_CGIs_hg38count_merge_uniq <- TNXB_hypo2000_CGIs_hg38count_merge_uniq[unique(TNXB_hypo2000_CGIs_hg38count_merge_uniq$V4),]
TNXB_hypo2000_CGIs_hg38count_merge_uniq <- TNXB_hypo2000_CGIs_hg38count_merge_uniq[order(-TNXB_hypo2000_CGIs_hg38count_merge_uniq$Density),]
TNXB_hypo2000_CGIs_hg38count_merge_uniq["Feature"] <- paste0("TNXB_hypo2000","_",
                                                             TNXB_hypo2000_CGIs_hg38count_merge_uniq$V4,
                                                             "%",
                                                             TNXB_hypo2000_CGIs_hg38count_merge_uniq$V2,
                                                             "%",
                                                             TNXB_hypo2000_CGIs_hg38count_merge_uniq$V3)

TNXB_hypo2000_CGIs_hg38count_merge_uniq <- TNXB_hypo2000_CGIs_hg38count_merge_uniq[,c(6,4,5)]
TNXB_hypo2000_CGIs_hg38count_merge_uniq["Type"] <- "TNXB_hpyo2000"
colnames(TNXB_hypo2000_CGIs_hg38count_merge_uniq) <- c("Feature", "Freq","Density","Type")

Monk_human_ICR_PCDH_TNXB_CGIs <- rbind.data.frame(Monk_human_ICR_CGIs_hg38count_merge_uniq, PCDHcluster_CGIs_hg38count_merge_uniq, TNXB_hypo2000_CGIs_hg38count_merge_uniq)
Monk_human_ICR_PCDH_TNXB_CGIs <- Monk_human_ICR_PCDH_TNXB_CGIs[order(-Monk_human_ICR_PCDH_TNXB_CGIs$Density),]
write.table(Monk_human_ICR_PCDH_TNXB_CGIs, "Monk_human_ICR_PCDH_TNXB_CGIs.txt", sep = "\t", append = F, quote = F, row.names = F)
ggplot(Monk_human_ICR_PCDH_TNXB_CGIs) + 
  geom_bar(aes(x=reorder(Feature, Density), y=Density, col=Feature, fill= Feature), stat = "identity", position = "dodge")+
  coord_flip()+theme_bw()+
  scale_fill_manual(values = c(rep("#C7CEEA",37),"darkgreen","orange"))+
  scale_color_manual(values = c(rep("white",50),"white","white"))

ggsave("Monk_human_ICR_PCDH_TNXB_CGIs.png", width=12, height=10, units="in", dpi=96)



#Calculate length of DMRs
head(MERGE_myiCombat3reavg_human_ICR,2)

MERGE_myiCombat3reavg_human_ICR_length <- MERGE_myiCombat3reavg_human_ICR[,c(1:4,23:27)]
MERGE_myiCombat3reavg_human_ICR_length["DMRLength"] <- MERGE_myiCombat3reavg_human_ICR_length$end - MERGE_myiCombat3reavg_human_ICR_length$start

#Assign Back ICRs
#PG
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep <- cSplit(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH, "row", "%")
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep[,c(4:8,1:3)]
colnames(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep) <- c("chr","start","end","TargetID","Additional","Sample","Value","Color")
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep) 

#Get pG specific ICRs length
diff_imp_loci_pG_Hypo_CpGs_density <- merge(diff_imp_loci_pG_Hypo_CpGs[,c(1:6)], MERGE_myiCombat3reavg_human_ICR_length, by="TargetID", all.x=T)
dim(diff_imp_loci_pG_Hypo_CpGs_density)
head(diff_imp_loci_pG_Hypo_CpGs_density,1)
diff_imp_loci_pG_Hypo_CpGs_density <- diff_imp_loci_pG_Hypo_CpGs_density[,c(4:6,1,2,15)]
count_diff_imp_loci_pG_Hypo_CpGs_density <- count(diff_imp_loci_pG_Hypo_CpGs_density,"DMR.x")
diff_imp_loci_pG_Hypo_CpGs_density_count <- merge(count_diff_imp_loci_pG_Hypo_CpGs_density, diff_imp_loci_pG_Hypo_CpGs_density, by="DMR.x", all.x=T)
diff_imp_loci_pG_Hypo_CpGs_density_count["Density"] <- (diff_imp_loci_pG_Hypo_CpGs_density_count$freq * 100)/diff_imp_loci_pG_Hypo_CpGs_density_count$DMRLength

miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr <- merge(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep, diff_imp_loci_pG_Hypo_CpGs_density_count, by="TargetID", all.y=T)
dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr)
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr$Color != "All"),]
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr[,c(6,7,8,9,15)]
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr)
colnames(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr) <- c("Sample","Value","Color","Feature","Density")
dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr)

#Make similar file PCDH for pG
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep$Color == "PCDH"),]
dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh)
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh)
length(unique(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh$TargetID))
PCDHcluster <- read.table("PCDHcluster.txt")
PCDHcluster["Color"] <- "PCDH"
PCDHcluster["Length"] <- PCDHcluster$V3 - PCDHcluster$V2
PCDHcluster["Density"] <- (length(unique(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh$TargetID)) *100) / PCDHcluster$Length
PCDHcluster
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh$Density <- PCDHcluster$Density
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh,1)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh[,c(6,7,8,8,9)]
colnames(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh) <- c("Sample","Value","Color","Feature","Density")
dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh)

miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh <- rbind.data.frame(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr, miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_pcdh)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Samples_Feature <- paste0(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Feature, "_",miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Sample)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Density_Feature <- paste0(round(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Density,2), "(",miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Feature,")")

head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh,1)


#Sort by density
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh[order(-miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Density),]
ggplot(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh, aes(x=reorder(Feature, -Density), y=Value, col=Sample, fill= Sample,group=Samples_Feature))+
  geom_point(position = position_dodge(width  = 0.3), size = 1,alpha = 1/3)+
  geom_hline(yintercept = c(0, 0), colour = "blue", linetype="dashed")+ 
  scale_color_manual(values=c("#EC7063","#48C9B0","#52BE80","#EC7063","#48C9B0","#52BE80"))+
  scale_fill_manual(values=c("white","white","white","white","white","white"))+ylim(c(-1,1))+
  theme_classic()+coord_flip()+scale_x_discrete(labels= c(unique(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh$Density_Feature)))


ggsave("miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH_sep_dmr_pcdh.svg", width=5, height=8, units="in", dpi=300)


#Assign Back ICRs
#PR
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep <- cSplit(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH, "row", "%")
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep[,c(4:8,1:3)]
colnames(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep) <- c("chr","start","end","TargetID","Additional","Sample","Value","Color")
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep) 

#Get pR specific ICRs length
diff_imp_loci_pR_Hypo_CpGs_density <- merge(diff_imp_loci_pR_Hypo_CpGs[,c(1:6)], MERGE_myiCombat3reavg_human_ICR_length, by="TargetID", all.x=T)
dim(diff_imp_loci_pR_Hypo_CpGs_density)
head(diff_imp_loci_pR_Hypo_CpGs_density,1)
diff_imp_loci_pR_Hypo_CpGs_density <- diff_imp_loci_pR_Hypo_CpGs_density[,c(4:6,1,2,15)]
count_diff_imp_loci_pR_Hypo_CpGs_density <- count(diff_imp_loci_pR_Hypo_CpGs_density,"DMR.x")
diff_imp_loci_pR_Hypo_CpGs_density_count <- merge(count_diff_imp_loci_pR_Hypo_CpGs_density, diff_imp_loci_pR_Hypo_CpGs_density, by="DMR.x", all.x=T)
diff_imp_loci_pR_Hypo_CpGs_density_count["Density"] <- (diff_imp_loci_pR_Hypo_CpGs_density_count$freq * 100)/diff_imp_loci_pR_Hypo_CpGs_density_count$DMRLength

miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr <- merge(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep, diff_imp_loci_pR_Hypo_CpGs_density_count, by="TargetID", all.y=T)
dim(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr)
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr[which(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr$Color != "All"),]
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr[,c(6,7,8,9,15)]
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr)
colnames(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr) <- c("Sample","Value","Color","Feature","Density")
dim(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr)

#Make similar file PCDH for pR
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep[which(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep$Color == "PCDH"),]
dim(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh)
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh)
length(unique(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh$TargetID))
PCDHcluster <- read.table("PCDHcluster.txt")
PCDHcluster["Color"] <- "PCDH"
PCDHcluster["Length"] <- PCDHcluster$V3 - PCDHcluster$V2
PCDHcluster["Density"] <- (length(unique(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh$TargetID)) *100) / PCDHcluster$Length
PCDHcluster
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh$Density <- PCDHcluster$Density
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh,1)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh[,c(6,7,8,8,9)]
colnames(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh) <- c("Sample","Value","Color","Feature","Density")
dim(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh)

miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh <- rbind.data.frame(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr, miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_pcdh)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Samples_Feature <- paste0(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Feature, "_",miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Sample)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Density_Feature <- paste0(round(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Density,2), "(",miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Feature,")")

head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh,1)
 

#Sort by density
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh[order(-miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Density),]
ggplot(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh, aes(x=reorder(Feature, -Density), y=Value, col=Sample, fill= Sample,group=Samples_Feature))+
  geom_point(position = position_dodge(width  = 0.3), size = 1,alpha = 1/3)+
  geom_hline(yintercept = c(0, 0), colour = "blue", linetype="dashed")+ 
  scale_color_manual(values=c("#922B21","#148F77","#1E8449","#922B21","#148F77","#1E8449"))+
  scale_fill_manual(values=c("white","white","white","white","white","white"))+ylim(c(-1,1))+
  theme_classic()+ coord_flip()+scale_x_discrete(labels= c(unique(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh$Density_Feature)))


ggsave("miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH_sep_dmr_pcdh.svg", width=5, height=8, units="in", dpi=300)

