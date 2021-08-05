### R code used to perform the differential expression analysis of DS2 ###
### The section numbers are the same as the ones in Materials and methods section ###

### 2.6 Differential expression analysis ###
### 2.6.1 FeatureCounts ###
#Input: Alignement files obtained with STAR
setwd("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021")

#FeatureCounts with countMultiMappingReads = TRUE, fraction= TRUE, strandSpecific=2 options
fc.mm2 <- featureCounts(files=c("map_55_Aligned.out.bam","map_56_Aligned.out.bam","map_57_Aligned.out.bam","map_62_Aligned.out.bam","map_63_Aligned.out.bam","map_64_Aligned.out.bam",
                                "map_52_Aligned.out.bam", "map_53_Aligned.out.bam", "map_54_Aligned.out.bam", "map_59_Aligned.out.bam", "map_60_Aligned.out.bam", "map_61_Aligned.out.bam",
                                "map_51_Aligned.out.bam", "map_65_Aligned.out.bam", "map_66_Aligned.out.bam", "map_58_Aligned.out.bam", "map_67_Aligned.out.bam", "map_68_Aligned.out.bam"),
                        isPairedEnd=TRUE,annot.ext="YPS128_Sbay_rebuilt_concatenated_TEs.gtf",isGTFAnnotationFile=TRUE,requireBothEndsMapped=TRUE,countMultiMappingReads = TRUE, fraction= TRUE, strandSpecific=2,
                        allowMultiOverlap = FALSE, countChimericFragments=FALSE, minFragLength=5, GTF.attrType="transcript_id")

# Extract raw count matrice  
raw_counts <- fc.mm2$counts
# Change column names 
colnames(raw_counts) <- c("hov_55_cer_30","hov_56_cer_30","hov_57_cer_30","hov_62_cer_12","hov_63_cer_12","hov_64_cer_12",
                          "hov_52_uv_30","hov_53_uv_30","hov_54_uv_30","hov_59_uv_12","hov_60_uv_12","hov_61_uv_12",
                          "hov_51_hyb_12", "hov_65_hyb_12", "hov_66_hyb_12", "hov_58_hyb_30","hov_67_hyb_30","hov_68_hyb_30")
#Export raw count matrice
write.table(raw_counts, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/raw_counts_tot_fev.mm.5.txt", sep="\t")

### 2.6.2 DESeq2 ###
library("DESeq2")
library("ggplot2")
library("Rsubread")
library("edgeR")
library("dplyr")
library("rlang")
library("stringi")
library("stringr") 
library("tidyr")
library("ggbeeswarm")
library("ggrepel")
library("rstatix")
library("ggpubr")
library("cowplot")
library("gridExtra")
library("tibble")

### 2.6.2.1 Multifactor analysis ###
raw_counts_tot <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/raw_counts_tot_fev.mm.5.txt", header = TRUE)
#Round down the raw counts  
for (i in (1:18)){
  raw_counts_tot[,i]<-round(raw_counts_tot[,i])
  raw_counts_tot[,i]<-as.integer(raw_counts_tot[,i])
}
raw_counts_tot[is.na(raw_counts_tot)] <- 0
#Select subgenomes 
raw_counts_cer_hyb <-raw_counts_tot[, c(1,2,3,4,5,6,13,14,15,16,17,18)]
raw_counts_cer_hyb <-raw_counts_cer_hyb %>% filter(str_detect(rownames(raw_counts_cer_hyb), 'YPS128')| str_detect(rownames(raw_counts_cer_hyb), 'TY'))
raw_counts_uv_hyb <-raw_counts_tot[, c(7,8,9,10,11,12,13,14,15,16,17,18)]
raw_counts_uv_hyb <-raw_counts_uv_hyb %>% filter(str_detect(rownames(raw_counts_uv_hyb), 'Sbay')|str_detect(rownames(raw_counts_uv_hyb), 'Tsu4'))
raw_counts_parents <-raw_counts_tot[, c(1,2,3,4,5,6,7,8,9,10,11,12)]

col_data_cer <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/table3_cer.txt", header = TRUE)
col_data_uv <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/table3_uv.txt", header = TRUE)
col_data_parents <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/table3_parents.txt", header = TRUE)
col_data <- col_data_cer

str(col_data)
col_data$group<-as.factor(col_data$group)
col_data$species<-as.factor(col_data$species)
col_data$temperature <-as.factor(col_data$temperature)
str(col_data)
col_data<- column_to_rownames(col_data, var = "id")

### 2.6.2.1.1 S. cerevisiae subgenome
raw_counts_dds <- raw_counts_cer_hyb
all(rownames(col_data)%in% colnames(raw_counts_dds))
dds <- DESeqDataSetFromMatrix(countData=raw_counts_dds, 
                              colData=col_data, 
                              design= ~ temperature + species + species:temperature)
#Set S. cerevisiae as the denominator, the reference species
dds$species <- relevel(dds$species, 'cer')
#Prefiltering 
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
#Run Deseq
dds<-DESeq(dds)
res <-results(dds) 
mcols(res, use.names=TRUE)
res
write.csv(as.data.frame(res), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/results_multifactor/res_tot.csv")
summary(res)
resultsNames(dds)

# Extract results of temperature effect 
res_temp <- results( dds,contrast = c("temperature", "12", "30") )
res_temp <-as.data.frame(res_temp)
res_temp <- tibble::rownames_to_column(res_temp, "gene")
temp_TEs <- res_temp %>% filter(str_detect(gene, 'TY'))
temp_TEs <-subset(temp_TEs, gene %in% list("TY1_I", "TY2_I", "TY3_I", "TY4_I" ))
temp_TEs <- with(temp_TEs,  temp_TEs[order(gene) , ])

# Extract results of species effect 
res_hyb_cer <- results( dds,contrast = c("species", "hyb", "cer") ) 
res_hyb_cer <-as.data.frame(res_hyb_cer)
res_hyb_cer <- tibble::rownames_to_column(res_hyb_cer, "gene")
hyb_cer_TEs <- res_hyb_cer %>% filter(str_detect(gene, 'TY'))
hyb_cer_TEs <-subset(hyb_cer_TEs, gene %in% list("TY1_I", "TY2_I", "TY3_I", "TY4_I" ))
hyb_cer_TEs <- with(hyb_cer_TEs,  hyb_cer_TEs[order(gene) , ])

# Extract results of the interaction between temperature and species 
res_intercept <- results( dds,name="temperature30.specieshyb") 
res_intercept <-as.data.frame(res_intercept)
res_intercept <- tibble::rownames_to_column(res_intercept, "gene")
intercept_TEs <- res_intercept %>% filter(str_detect(gene, 'TY'))
intercept_TEs <-subset(intercept_TEs, gene %in% list("TY1_I", "TY2_I", "TY3_I", "TY4_I" ))

# Export TEs normalized data
norm_counts_cer_30_12 <- counts(dds, normalized=T) 
norm_counts_cer_30_12 <- as.data.frame(norm_counts_cer_30_12)
norm_counts_cer_30_12 <- tibble::rownames_to_column(norm_counts_cer_30_12, "gene")
write.table(norm_counts_cer30_12, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/results/results_multifactor/norm_counts_cer30_12.txt", sep="\t")
norm_counts_cer_30_12_TEs <-norm_counts_cer_30_12 %>% filter(str_detect(gene, 'TY'))
write.table(norm_counts_cer_30_12_TEs, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/results_multifactor/norm_counts_cer_30_12_TEs.txt", sep="\t")

### 2.6.2.1.2 S. uvarum subgenome
col_data <- col_data_uv
str(col_data)
col_data$group<-as.factor(col_data$group)
col_data$species<-as.factor(col_data$species)
col_data$temperature <-as.factor(col_data$temperature)
str(col_data)
col_data<- column_to_rownames(col_data, var = "id")

raw_counts_dds <- raw_counts_uv_hyb
all(rownames(col_data)%in% colnames(raw_counts_dds))
dds <- DESeqDataSetFromMatrix(countData=raw_counts_dds, 
                              colData=col_data, 
                              design= ~ temperature + species + species:temperature)
# Set uvarum as denominator 
dds$species <- relevel(dds$species, 'uv')
#Prefiltering 
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
#Run Deseq
dds<-DESeq(dds)
res <-results(dds) 
mcols(res, use.names=TRUE)
res
summary(res)
resultsNames(dds)

# Extract results of temperature effect 
res_temp_uv <- results( dds,contrast = c("temperature", "12", "30") )
res_temp_uv <-as.data.frame(res_temp_uv)
res_temp_uv <- tibble::rownames_to_column(res_temp_uv, "gene")
temp_TEs_uv <- res_temp_uv %>% filter(str_detect(gene, 'Tsu4'))
temp_TEs_uv <-subset(temp_TEs_uv, gene %in% list("Tsu4-b-I" ))

#combine results of Scer and Suva subgenomes  
res_temp_Scer_Suva <- bind_rows(res_temp, res_temp_uv)
temp_TEs_Scer_Suva <- bind_rows(temp_TEs, temp_TEs_uv)

volcano_temp_hovhan2 <- ggplot(data=res_temp_Scer_Suva, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=1.0, color="light grey")+ theme_classic() +
  ggtitle("DS2 (Sc x Su) : \ntemperature effect")+
  theme(plot.title = element_text( face = "bold", size=10.5))+
  ggeasy::easy_center_title()+
  xlab("log2 fold change (12°C/30°C)")+
  ylab("-log10 (adjusted p-value)")+
  scale_x_continuous(breaks=seq(-4,4,1), limits=c(-2,2)) +
  scale_y_continuous(breaks=seq(0,8,2),limits=c(0,8))+
  geom_point(data=temp_TEs_Scer_Suva,aes(x=log2FoldChange, y=-log10(padj), shape=gene), color='#8117CE', size=3.5) +
  geom_vline(xintercept=c(0), col="black", size=0.7, linetype="dashed")+
  geom_hline(yintercept=c(-log10(0.01),-log10(0.05)), col="black", size=0.7, linetype="dashed")+
  theme(legend.position = "top")+
  guides(shape=guide_legend(override.aes = list(size=2.5),nrow=3, ncol=2,byrow=FALSE, title = ""))+
  scale_shape_manual(name = "Legend",
                     labels = c( "Tsu4_uva","Ty1_cer","Ty2_cer" ,"Ty3_cer","Ty4_cer"),
                     values =  c(8,15, 16, 17, 18))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(text=element_text(size=9))
volcano_temp_hovhan2

volcano_temp_hovhan2_subplot <- ggplot(data=res_temp_Scer_Suva, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=0.5, color="black")+ theme_classic() + 
  xlab("")+
  ylab("")+
  scale_x_continuous(breaks=seq(-8,8,8), limits=c(-8,8)) +
  scale_y_continuous(breaks=seq(0,80,80),limits=c(0,80))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(panel.background = element_rect(colour = NA, fill = NA), plot.background = element_rect(colour = NA, fill = NA))
volcano_temp_hovhan2_subplot

volcano_temp_hovhan2_complet<-
  ggdraw() +
  draw_plot(volcano_temp_hovhan2) +
  draw_plot(volcano_temp_hovhan2_subplot, x = 0.60, y = .33, width = .4, height = .3, scale = 1)
volcano_temp_hovhan2_complet

# Extract results of species effect 
res_hyb_uv <- results( dds,contrast = c("species", "hyb", "uv") ) 
res_hyb_uv <-as.data.frame(res_hyb_uv)
res_hyb_uv <- tibble::rownames_to_column(res_hyb_uv, "gene")
hyb_uv_TEs <- res_hyb_uv %>% filter(str_detect(gene, 'Tsu4'))
hyb_uv_TEs <-subset(hyb_uv_TEs, gene %in% list("Tsu4-b-I"))

res_hyb_Scer_Suv <- bind_rows(res_hyb_cer, res_hyb_uv)
hyb_Scer_Suv_TEs <- bind_rows(hyb_cer_TEs,hyb_uv_TEs)
write.csv(as.data.frame(hyb_Scer_Suv_TEs), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/hyb_Scer_Suv_TEs.csv")
 
volcano_hyb_Scer_hovhan <- ggplot(data=res_hyb_Scer_Suv, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=1.0, color="light grey")+ theme_classic() +
  ggtitle("DS2 (Sc x Su) : \nspecies effect")+
  theme(plot.title = element_text( face = "bold", size=10.5))+
  ggeasy::easy_center_title()+
  xlab("log2 fold change (hybrid/parent)")+
  ylab("-log10 (adjusted p-value)")+
  scale_x_continuous(breaks=seq(-4,4,1), limits=c(-2,2)) +
  scale_y_continuous(breaks=seq(0,8,2),limits=c(0,8))+
  geom_point(data=hyb_Scer_Suv_TEs,aes(x=log2FoldChange, y=-log10(padj), shape=gene), color='#8117CE', size=3.5) +
  geom_vline(xintercept=c(0), col="black", size=0.7, linetype="dashed")+
  geom_hline(yintercept=c(-log10(0.01),-log10(0.05)), col="black", size=0.7, linetype="dashed")+
  theme(legend.position = "top",legend.box = "vertical")+
  guides(shape=guide_legend(override.aes = list(size=2.5),nrow=3, ncol=2,byrow=FALSE, title = ""))+
  scale_shape_manual(name = "Legend",
                     labels = c("Tsu4_uva","Ty1_cer","Ty2_cer" ,"Ty3_cer","Ty4_cer" ),
                     values =  c(8,15, 16, 17, 18))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(text=element_text(size=9))
volcano_hyb_Scer_hovhan

volcano_hyb_Scer_hovhan_subplot <- ggplot(data=res_hyb_Scer_Suv, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=0.5, color="black")+ theme_classic() +
  xlab("")+
  ylab("")+
  scale_x_continuous(breaks=seq(-8,8,8), limits=c(-8,8)) +
  scale_y_continuous(breaks=seq(0,50,50),limits=c(0,50))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(panel.background = element_rect(colour = NA, fill = NA), plot.background = element_rect(colour = NA, fill = NA))
volcano_hyb_Scer_hovhan_subplot

volcano_hyb_Scer_hovhan_complet<-
  ggdraw() +
  draw_plot(volcano_hyb_Scer_hovhan) +
  draw_plot(volcano_hyb_Scer_hovhan_subplot, x = 0.60, y = .33, width = .4, height = .3, scale = 1)
volcano_hyb_Scer_hovhan_complet

#Extract results of interaction between temperature and species 
res_intercept2 <- results( dds,name="temperature30.specieshyb") 
res_intercept2 <-as.data.frame(res_intercept2)
res_intercept2 <- tibble::rownames_to_column(res_intercept2, "gene")
intercept_TEs2 <- res_intercept2 %>% filter(str_detect(gene, 'Tsu4'))
intercept_TEs2 <-subset(intercept_TEs2, gene %in% list("Tsu4-b-I" ))

res_intercept_tot <- bind_rows(res_intercept, res_intercept2)
intercept_TEs_tot <- bind_rows(intercept_TEs,intercept_TEs2)

volcano_intercept <- ggplot(data=res_intercept_tot, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=1.0, color="light grey")+ theme_classic() +
  ggtitle("DS2 (Sc x Su) : \n temperature-species interaction")+
  theme(plot.title = element_text( face = "bold", size=10.5))+
  ggeasy::easy_center_title()+
  xlab("log2 fold change (temperature/species)")+
  ylab("-log10 (adjusted p-value)")+
  scale_x_continuous(breaks=seq(-4,4,1), limits=c(-2,2)) +
  scale_y_continuous(breaks=seq(0,8,2),limits=c(0,8))+
  geom_point(data=intercept_TEs_tot,aes(x=log2FoldChange, y=-log10(padj), shape=gene), color='#8117CE', size=3.5) +
  geom_vline(xintercept=c(0), col="black", size=0.7, linetype="dashed")+
  geom_hline(yintercept=c(-log10(0.01),-log10(0.05)), col="black", size=0.7, linetype="dashed")+
  theme(legend.position = "top",legend.box = "vertical")+
  guides(shape=guide_legend(override.aes = list(size=2.5),nrow=3, ncol=2,byrow=FALSE, title = "",size=2.5))+
  scale_shape_manual(name = "Legend",
                     labels = c("Tsu4_uva","Ty1_cer","Ty2_cer" ,"Ty3_cer","Ty4_cer" ),
                     values =  c(8,15, 16, 17, 18))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(text = element_text(size=9))
volcano_intercept

volcano_intercept_subplot <- ggplot(data=res_intercept_tot, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=0.5, color="black")+ theme_classic() +
  xlab("")+
  ylab("")+
  scale_x_continuous(breaks=seq(-6,6,6), limits=c(-6,6)) +
  scale_y_continuous(breaks=seq(0,12,12),limits=c(0,12))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(panel.background = element_rect(colour = NA, fill = NA), plot.background = element_rect(colour = NA, fill = NA))
volcano_intercept_subplot

volcano_intercept_complet<-
  ggdraw() +
  draw_plot(volcano_intercept) +
  draw_plot(volcano_intercept_subplot, x = 0.60, y = .33, width = .4, height = .3, scale = 1)
volcano_intercept_complet

figure_volcanoplot_hovhan_CDE2 <- plot_grid(volcano_hyb_Scer_hovhan_complet, volcano_temp_hovhan2_complet, volcano_intercept_complet,align = "h", nrow = 1, rel_widths = c(1/3, 1/3,1/3),labels = c('C', 'D','E'))
figure_volcanoplot_hovhan_CDE2
ggsave(filename = "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/figure_volcanoplot_hovhan_CDE2.2.jpg", 
       plot = figure_volcanoplot_hovhan_CDE2,
       width = 20, 
       height = 10,
       units = "cm",
       dpi = 300)

### 2.6.2.2 Separated analysis for each temperature ###
raw_counts_tot <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/raw_counts_tot_fev.mm.5.txt", header = TRUE)
for (i in (1:18)){
  raw_counts_tot[,i]<-round(raw_counts_tot[,i])
  raw_counts_tot[,i]<-as.integer(raw_counts_tot[,i])
}
raw_counts_tot[is.na(raw_counts_tot)] <- 0
raw_counts_cer_hyb <-raw_counts_tot[, c(1,2,3,4,5,6,13,14,15,16,17,18)]
raw_counts_cer_hyb <-raw_counts_cer_hyb %>% filter(str_detect(rownames(raw_counts_cer_hyb), 'YPS128')| str_detect(rownames(raw_counts_cer_hyb), 'TY'))
raw_counts_uv_hyb <-raw_counts_tot[, c(7,8,9,10,11,12,13,14,15,16,17,18)]
raw_counts_uv_hyb <-raw_counts_uv_hyb %>% filter(str_detect(rownames(raw_counts_uv_hyb), 'Sbay')| str_detect(rownames(raw_counts_cer_hyb), 'Tsu4'))

col_data_cer <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/table3_cer.txt", header = TRUE)
col_data <-col_data_cer
str(col_data)
col_data$group<-as.factor(col_data$group)
col_data$species<-as.factor(col_data$species)
col_data$temperature<-as.factor(col_data$temperature)
str(col_data)
col_data <- tibble::column_to_rownames(col_data, var = "id")

# Analysis with only the species factor 
raw_counts_tot <- raw_counts_cer_hyb
all(rownames(col_data)%in% colnames(raw_counts_tot))
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = raw_counts_tot,
  colData = col_data,
  design = ~ species)
ddsFullCountTable

### Redo DESeq2 analysis with samples grown at 30 degrees and than with samples grown at 12 degrees 
# with the same steps as in DS1 differential expression analysis 
#Select samples at different temperature with: 
dds30 <- ddsFullCountTable[ , ddsFullCountTable$temperature == "30" ]
as.data.frame( colData(dds30) )
dds12 <- ddsFullCountTable[ , ddsFullCountTable$temperature == "12" ]
as.data.frame( colData(dds12) )

#Create figure 3 
figure_volcanoplot_ABCDE <- plot_grid(figure_volcanoplot3, figure_volcanoplot_hovhan_CDE2, align = "v", nrow = 2,rel_heights = c(1.3/3, 1.7/3))
figure_volcanoplot_ABCDE
ggsave(filename = "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/figure_volcanoplot_ABCDE.jpg", 
       plot = figure_volcanoplot_ABCDE,
       width = 20, 
       height = 17,
       units = "cm",
       dpi = 300)

# Create figure 4
figure_Schraiber_hovhan_ABCD <- plot_grid(figure_Schraiber, Figure1_Hovhan, align = "v", nrow = 2)
figure_Schraiber_hovhan_ABCD
ggsave(filename = "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/figure_Schraiber_hovhan_ABCD.jpg", 
       plot = figure_Schraiber_hovhan_ABCD,
       width = 20, 
       height = 18,
       units = "cm",
       dpi = 300)

### 2.7 TPM calculation ###
gene_lenght_hovhan <-cbind(fc.mm2$annotation[1],fc.mm2$annotation[6])
colnames(gene_lenght_hovhan) <- c("gene","length")
#add lenght column to the raw counts file
raw_counts_hovhan_tot <-merge(raw_counts_TPM, gene_lenght_hovhan, by="gene",all = T)
write.table(raw_counts_hovhan_tot, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/counts_hovhan.5.length.txt", sep="\t")

#TPM normalization
raw_counts_hovhan_5 <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/counts_hovhan.5.length.txt", header = TRUE)
raw_counts_hovhan_5_TPM <-as.data.frame(raw_counts_hovhan_5)

for(i in colnames(raw_counts_hovhan_5_TPM[2:19])){
  raw_counts_hovhan_5_TPM[,i]<-raw_counts_hovhan_5_TPM[,i]/(raw_counts_hovhan_5_TPM[,'length']/1000)
  raw_counts_hovhan_5_TPM[,i]<-raw_counts_hovhan_5_TPM[,i]/sum(raw_counts_hovhan_5_TPM[,i], na.rm=TRUE)*1e6
}
write.table(raw_counts_hovhan_5_TPM, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/counts_hovhan_cer_bay_TPM.txt", sep="\t")
raw_counts_hovhan_5_TPM_TEs <- raw_counts_hovhan_5_TPM %>% filter(str_detect(gene, 'Tsu4')| str_detect(gene, 'TY'))
write.table(raw_counts_hovhan_5_TPM_TEs, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/results_nov_2020/dec_2020/fev_2021/counts_hovhan_cer_bay_TPM_TEs.txt", sep="\t")
