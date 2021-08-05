### R code used to perform the differential expression analysis of DS1 ###
### The section numbers are the same as the ones in Materials and methods section ###

### 2.6 Differential expression analysis ###
### 2.6.1 FeatureCounts ###
#Input: Alignement files obtained with Bowtie2 --local
library(Rsubread)
library(edgeR)
library(dplyr)
# 2.6.1.1 S. uvarum 
fc.cer.bay4 <- featureCounts(files=c("index_bowtie2_polyT_cer_bay_conc/IGV2/SRR515228_gathered2.TY.sorted.bam", "index_bowtie2_polyT_cer_bay_conc/IGV2/SRR515229_gathered2.TY.sorted.bam", 
                                     "index_bowtie2_polyT_cer_bay_conc/IGV2/SRR515230_gathered2.TY.sorted.bam","index_bowtie2_polyT_cer_bay_conc/IGV2/SRR515231_gathered2.TY.sorted.bam",
                                     "index_bowtie2_polyT_cer_bay_conc/IGV2/SRR835217_gathered2.TY.sorted.bam", "index_bowtie2_polyT_cer_bay_conc/IGV2/SRR835218_gathered2.TY.sorted.bam"),
                             isPairedEnd=TRUE,annot.ext="genome/Scer.Sbay.rebuilt.TEs.conc.gtf",isGTFAnnotationFile=TRUE,requireBothEndsMapped=TRUE,countMultiMappingReads = FALSE, 
                             allowMultiOverlap = FALSE, countChimericFragments=FALSE, minFragLength=5, GTF.attrType="transcript_id")

#Extraction of the raw count matrice 
raw_counts <- fc.cer.bay4$counts
#Change column names
raw_counts <- as.data.frame(raw_counts)
raw_counts <- tibble::rownames_to_column(raw_counts, "gene")
colnames(raw_counts) <- c("gene", "SRR515228_cer", "SRR515229_cer","SRR515230_bay","SRR515231_bay","SRR835217_hyb", "SRR835218_hyb")
#Exportation of the raw counts matrice
write.table(raw_counts, "/home/madro269/schraiber2013/index_bowtie2_polyT_cer_bay_conc/results_FC2/polyT_trim_counts_cer_bay.mm.4.txt", sep="\t")
#Results extraction
write.table(fc.cer.bay4$stat, "/home/madro269/schraiber2013/index_bowtie2_polyT_cer_bay_conc/results_FC2/fc.cer.bay.stat.mm.4.txt", sep="\t")
write.csv(as.data.frame(fc.cer.bay4$stat), file="/home/madro269/schraiber2013/index_bowtie2_polyT_cer_bay_conc/results_FC2/fc.cer.bay.stat.mm.4.csv")

# 2.6.1.1 S. paradoxus : repeat the same steps 
fc.cer.par4 <- featureCounts(files=c("index_bowtie2_polyT_cer_par_conc/IGV2/SRR515220_gathered2.TY.sorted.bam", "index_bowtie2_polyT_cer_par_conc/IGV2/SRR515221_gathered2.TY.sorted.bam", 
                                     "index_bowtie2_polyT_cer_par_conc/IGV2/SRR515222_gathered2.TY.sorted.bam","index_bowtie2_polyT_cer_par_conc/IGV2/SRR515223_gathered2.TY.sorted.bam",
                                     "index_bowtie2_polyT_cer_par_conc/IGV2/SRR835213_gathered2.TY.sorted.bam", "index_bowtie2_polyT_cer_par_conc/IGV2/SRR835214_gathered2.TY.sorted.bam"),
                             isPairedEnd=TRUE,annot.ext="genome/Scer.Spar.rebuilt.TEs.conc.2.gtf",isGTFAnnotationFile=TRUE,requireBothEndsMapped=TRUE,countMultiMappingReads = FALSE, 
                             allowMultiOverlap = FALSE, countChimericFragments=FALSE, minFragLength=5, GTF.attrType="transcript_id")
raw_counts <- fc.cer.par4$counts
raw_counts <- as.data.frame(raw_counts)
raw_counts <- tibble::rownames_to_column(raw_counts, "gene")
colnames(raw_counts) <- c("gene", "SRR515220_cer", "SRR515221_cer","SRR515222_par","SRR515223_par","SRR835213_hyb", "SRR835214_hyb")
write.table(raw_counts, "/home/madro269/schraiber2013/index_bowtie2_polyT_cer_par_conc/results_FC2/polyT_trim_counts_cer_par.mm.4.txt", sep="\t")
write.table(fc.cer.par4$stat, "/home/madro269/schraiber2013/index_bowtie2_polyT_cer_par_conc/results_FC2/fc.cer.par.stat.mm.4.txt", sep="\t")
write.csv(as.data.frame(fc.cer.par4$stat), file="/home/madro269/schraiber2013/index_bowtie2_polyT_cer_par_conc/results_FC2/fc.cer.par.stat.mm.4.csv")

### 2.6.2 Deseq2 ###
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

###2.6.2.1 Scer-hyb Analysis ####
#Import raw counts
raw_counts_tot <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/polyT_trim_counts_cer_bay.mm.4.txt", header = TRUE)
#Replace Na by 0 because DESeq doesn't accept NA 
raw_counts_tot[is.na(raw_counts_tot)] <- 0
#Select appropriate columns for each parental-hybrid comparison
raw_counts_cer_hyb <-raw_counts_tot[, c(1,2,5,6)]
raw_counts_cer_hyb <-raw_counts_cer_hyb %>% filter(str_detect(rownames(raw_counts_cer_hyb), 'Scer')| str_detect(rownames(raw_counts_cer_hyb), 'TY'))
raw_counts_uv_hyb <-raw_counts_tot[, c(3,4,5,6)]
raw_counts_uv_hyb <-raw_counts_uv_hyb %>% filter(str_detect(rownames(raw_counts_uv_hyb), 'Sbay')|str_detect(rownames(raw_counts_uv_hyb), 'Tsu4'))
#Import col_data 
col_data_cer <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/table_Scer_Sbay_Scer.txt", header = TRUE)
col_data <-col_data_cer
#Make sure group, species and temperature are factors
str(col_data)
col_data$species<-as.factor(col_data$species)
str(col_data)
#Format the dataframe for DESeq2
col_data <- tibble::column_to_rownames(col_data, var = "id")

### Scer-hyb Analysis ####
raw_counts_tot <- raw_counts_cer_hyb
all(rownames(col_data)%in% colnames(raw_counts_tot))
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_tot,
  colData = col_data,
  design = ~ species)
dds
as.data.frame(colData(dds))

#Running the pipeline
dds$species <- relevel(dds$species, 'cer')
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds
dds <- DESeq(dds)
res <- results(dds)
res
mcols(res, use.names=TRUE)
summary(res)
sizeFactors(dds)

# Export normalized counts of TEs
norm_counts_cer <- counts(dds, normalized=T)
norm_counts_cer <- as.data.frame(norm_counts_cer)
norm_counts_cer_2 <- tibble::rownames_to_column(norm_counts_cer, "gene")
write.table(norm_counts_cer_2, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/norm_counts_Scer_hyb_Sc_Su.txt", sep="\t")
#Filter to keep all TEs that begin with TY and Tsu4
norm_counts_cer_TEs <-norm_counts_cer_2 %>% filter(str_detect(gene, 'TY')| str_detect(gene, 'Tsu4'))
write.table(norm_counts_cer_TEs, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/norm_counts_Scer_hyb_TEs_cb.txt", sep="\t")

#Optional analysis with resOrdered 
resOrdered <- res[order(res$padj),]
resOrdered.log2 <- res[order(res$log2FoldChange),]
resOrdered.log2 <-as.data.frame(resOrdered.log2)
resOrdered.log2.2 <- tibble::rownames_to_column(resOrdered.log2, "gene")
write.csv(as.data.frame(resOrdered.log2.2), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Scer_hyb_cb_resOrdered.csv")
hyb_cer_TEs <- resOrdered.log2.2 %>% filter(str_detect(gene, 'TY'))
write.csv(as.data.frame(hyb_cer_TEs), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Scer_hyb_TEs_cb_resOrdered.csv")

# Extract data to show in ggplot : Scer of Scer-Sbay
norm_counts_TEs <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/norm_counts_Scer_hyb_TEs_cb.txt", header = TRUE)
# Subdivise column info in different columns 
norm_counts_TEs_long<-norm_counts_TEs %>%
  gather(info, norm_counts_TEs, SRR515228_cer:SRR835218_hyb)%>%
  separate(info, into = c("replicate", "species"), sep = "_")
colnames(norm_counts_TEs_long)<-c("gene", "replicate", "species","norm.counts")
#Keep only TEs to show on graph 
norm_counts_TEs_long_I <-subset(norm_counts_TEs_long, gene %in% list("TY1_cer", "TY2_cer", "TY3_cer", "TY4_cer", "TY3_par", "TY1_par", "TY5_par" ))
norm.counts.cer.cb <- norm_counts_TEs_long_I

###2.6.2.2 Suva-hyb Analysis ####
#We redo the same thing with Suva of Scer-Suva

raw_counts_tot <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/polyT_trim_counts_cer_bay.mm.4.txt", header = TRUE)
raw_counts_tot[is.na(raw_counts_tot)] <- 0
raw_counts_cer_hyb <-raw_counts_tot[, c(1,2,5,6)]
raw_counts_uv_hyb <-raw_counts_tot[, c(3,4,5,6)]

col_data_uv <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/table_Scer_Sbay_Sbay.txt", header = TRUE)
col_data <-col_data_uv

str(col_data)
col_data$species<-as.factor(col_data$species)
str(col_data)
col_data <- tibble::column_to_rownames(col_data, var = "id")

### Suva-hyb Analysis ####
raw_counts_tot <- raw_counts_uv_hyb
all(rownames(col_data)%in% colnames(raw_counts_tot))
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts_tot,
  colData = col_data,
  design = ~ species)
dds
as.data.frame(colData(dds))

#Running the pipeline
dds$species <- relevel(dds$species, 'bay')
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds
dds <- DESeq(dds)
res <- results(dds)
res
mcols(res, use.names=TRUE)
summary(res)
dds$sizeFactor

# Export normalize data
norm_counts_cer <- counts(dds, normalized=T)
norm_counts_cer <- as.data.frame(norm_counts_cer)
norm_counts_cer_2 <- tibble::rownames_to_column(norm_counts_cer, "gene")
write.table(norm_counts_cer_2, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/norm_counts_Sbay_hyb_Sc_Su.txt", sep="\t")
norm_counts_cer_TEs <-norm_counts_cer_2 %>% filter(str_detect(gene, 'TY')| str_detect(gene, 'Tsu4'))
write.table(norm_counts_cer_TEs, "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/norm_counts_Sbay_hyb_TEs_cb.txt", sep="\t")

#Optional analysis
resOrdered <- res[order(res$padj),]
resOrdered.log2 <- res[order(res$log2FoldChange),]
resOrdered.log2 <-as.data.frame(resOrdered.log2)
resOrdered.log2.2 <- tibble::rownames_to_column(resOrdered.log2, "gene")
write.csv(as.data.frame(resOrdered.log2.2), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Sbay_hyb_cb_resOrdered.csv")
hyb_cer_TEs <- resOrdered.log2.2 %>% filter(str_detect(gene, 'Tsu4'))
write.csv(as.data.frame(hyb_cer_TEs), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Sbay_hyb_TEs_cb_resOrdered.csv")

# Extract data to show in ggplot : Sbay of Scer-Sbay 
norm_counts_TEs <- read.table("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/norm_counts_Sbay_hyb_TEs_cb.txt", header = TRUE)
norm_counts_TEs_long<-norm_counts_TEs %>%
  gather(info, norm_counts_TEs, SRR515230_bay:SRR835218_hyb)%>%
  separate(info, into = c("replicate", "species"), sep = "_")
colnames(norm_counts_TEs_long)<-c("gene", "replicate", "species","norm.counts")
norm_counts_TEs_long_I <-subset(norm_counts_TEs_long, gene %in% list("Tsu4_bay"))
norm.counts.bay.cb <- norm_counts_TEs_long_I

###2.6.2.3 Concatenate the 2 dataframes ####
norm.counts.tot <- bind_rows(norm.counts.cer.cb,norm.counts.bay.cb)
norm.counts.tot.test <- norm.counts.tot

norm.counts.tot.test$species <- as.character(norm.counts.tot.test$species)
norm.counts.tot.test[norm.counts.tot.test == "bay"] <- "uv"
norm.counts.tot.test[norm.counts.tot.test == "Tsu4_bay"] <- "Tsu4_uva"
norm.counts.tot.test[norm.counts.tot.test == "cer"] <- "Sc"
norm.counts.tot.test[norm.counts.tot.test == "hyb"] <- "Sc x Su"
norm.counts.tot.test[norm.counts.tot.test == "uv"] <- "Su"
norm.counts.tot.test

###2.6.2.4 ggplot Graph Scer-Sbay####

TEs_names <- list(
  "TY1_cer"="Ty1_cer",
  "TY2_cer"="Ty2_cer",
  "TY3_cer"="Ty3_cer",
  "TY3_cer"="Ty3_par",
  "TY4_cer"="Ty4_cer",
  "TY5_cer"="Ty5_cer",
  "Tsu4_uva"="Tsu4_uva")
TEs_labeller <- function(variable,value){
  return(TEs_names[value])
}
stat.test <- stat.test %>% add_xy_position(x = "species")
bxp + stat_pvalue_manual(stat.test)
annotation_df1 <- data.frame(
  gene = c("Tsu4_uva", "TY1_cer","TY2_cer","TY3_cer","TY4_cer"),
  start = c("Sc x Su", "Sc","Sc","Sc","Sc"),
  end = c("Su", "Sc x Su","Sc x Su","Sc x Su","Sc x Su"),
  y = c(14.5,14.5,14.5,14.5,14.5),
  label = c("0.0474", "0.000195", "0.789", "0.616", "0.786")
)
#Reorder sample in each facet
norm.counts.tot.test$species = factor(norm.counts.tot.test$species, levels=c("Sc","Su","Sc x Su"), labels=c("Sc","Su","Sc x Su")) 

Graph_Schraiber_Scer_Sbay <- ggplot(data = norm.counts.tot.test,
                                    aes(x = species, y = log2(norm.counts+1),
                                        color = species)) +
  scale_shape_manual(values=c(1, 19))+
  scale_colour_manual(name= '',labels = c("Sc", "Su", "Sc x Su"),values = c("#0925CF","red1", "darkviolet"))+
  geom_beeswarm(size = 2, cex=10)+
  facet_grid(~gene, labeller = TEs_labeller,scales = "free_x")+
  theme_test() +
  scale_y_continuous(breaks=seq(0,18,2), limits=c(0,16))+
  scale_fill_manual(values = c("white", "grey", "black")) +
  ylab("log2 (Normalized Ty counts)")+
  xlab("Genome")+
  labs(family="Arial", color='black')+
  theme(legend.position="top")+
  ggtitle ("DS1 : Sc x Su")+
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"))+
  theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=0.3, color='black'))+
  theme(axis.text.y = element_text(color='black'))+
  theme(text = element_text(size = 10))+
  geom_signif(
    inherit.aes = FALSE,
    data = annotation_df1,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 2.2, vjust = -0.2,
    manual = TRUE)
Graph_Schraiber_Scer_Sbay

#### 2.6.2.5 volcanoplot ###
Sbay <- read.csv("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Sbay_hyb_cb_resOrdered.csv", header = TRUE)
Scer <- read.csv("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Scer_hyb_cb_resOrdered.csv", header = TRUE)
Scer_Sbay <- bind_rows(Scer, Sbay)
Scer_Sbay$Ty <- Scer_Sbay$gene
Scer_Sbay$Ty <- read.fwf(textConnection(Scer_Sbay$Ty), 4)
str(Scer_Sbay$Ty)
Sbay_TEs <- read.csv("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Sbay_hyb_TEs_cb_resOrdered.csv", header = TRUE)
Scer_TEs <- read.csv("C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/trim_bowtie2_cer_bay/Scer_hyb_TEs_cb_resOrdered.csv", header = TRUE)
Scer_Sbay_TEs <-bind_rows(Scer_TEs, Sbay_TEs)
Scer_Sbay_TEs$name <- c("","Ty4_cer", "Ty1_cer", "Ty3_cer", "Tsu4_uva")
write.csv(as.data.frame(Scer_Sbay_TEs), 
          file="C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/Scer_Sbay_TEs.csv")

volcano_Scer_Sbay2 <- ggplot(data=Scer_Sbay, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=1.0, color="light grey")+ theme_classic() +
  ggtitle("DS1 : Sc x Su")+
  theme(plot.title = element_text( face = "bold", size=10.5))+
  ggeasy::easy_center_title()+
  xlab("log2 fold change (hybrid/parent)")+
  ylab("-log10 (adjusted p-value)")+
  scale_x_continuous(breaks=seq(-4,4,1), limits=c(-2,2)) +
  scale_y_continuous(breaks=seq(0,10,1),limits=c(0,8))+
  geom_point(data=Scer_Sbay_TEs,aes(x=log2FoldChange, y=-log10(padj), shape=gene), color='#8117CE', size=3.1) +
  geom_vline(xintercept=c(0), col="black", size=0.7, linetype="dashed")+
  geom_hline(yintercept=c(-log10(0.01),-log10(0.05)), col="black", size=0.7, linetype="dashed")+
  theme(legend.position = "right",legend.box = "vertical")+
  guides(shape=guide_legend(override.aes = list(size=2.5),nrow=5, ncol=1,byrow=TRUE, title = ""))+
  scale_shape_manual(name = "Legend",
                     labels = c("Tsu4_uva","Ty1_cer","Ty2_cer" ,"Ty3_cer","Ty4_cer" ),
                     values =  c(8, 15, 16, 17, 18))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))+
  theme(text=element_text(size=9))
volcano_Scer_Sbay2

Scer_Sbay_subplot <-subset(Scer_Sbay, padj>-0.00)
volcano_Scer_Sbay2_subplot <- ggplot(data=Scer_Sbay_subplot, aes(x=log2FoldChange, y=-log10(padj)))+geom_point(size=0.5, color="black")+ theme_classic() +
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks=seq(0,300,300),limits=c(0,350),expand = c(0, 0))+
  theme(panel.background = element_rect(colour = NA, fill = NA), plot.background = element_rect(colour = NA, fill = NA))+
  theme(axis.text.x=element_text(color="black"))+
  theme(axis.text.y=element_text(color="black"))
volcano_Scer_Sbay2_subplot

volcano_Scer_Sbay2_complet<-
  ggdraw() +
  draw_plot(volcano_Scer_Sbay2) +
  draw_plot(volcano_Scer_Sbay2_subplot, x = 0.42, y = .50, width = .3, height = .4, scale = 1)
volcano_Scer_Sbay2_complet

ggsave(filename = "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/DS1_ScxSu.png", 
       plot = volcano_Scer_Sbay2_complet,
       width = 10, 
       height = 12,
       units = "cm",
       dpi = 300)

### The same differential expression analysis using DESeq2 (section 2.6.2) was done with S. cerevisiae x S. paradoxus samples ### 

# Create final figures 
figure_Schraiber <- plot_grid(Graph_Schraiber_Scer_Sbay, Graph_Schraiber_Scer_Spar, align = "h", nrow = 1, rel_widths = c(6.5/14, 9/14),labels = c('A', 'B'))
figure_Schraiber
ggsave(filename = "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/figure_Schraiber.jpg", 
       plot = figure_Schraiber,
       width = 20, 
       height = 10,
       units = "cm",
       dpi = 300)

figure_volcanoplot3 <- plot_grid(volcano_Scer_Sbay2_complet, volcano_Scer_Spar2_complet, align = "h", nrow = 1, rel_widths = c(1/2, 1/2),labels = c('A', 'B'))
figure_volcanoplot3
ggsave(filename = "C:/Users/marik/OneDrive/Documents/Université/Maitrise/RNAseq/2e_dataset/bowtie2_filtered_conc/Figures/figure_volcanoplot3.jpg", 
       plot = figure_volcanoplot3,
       width = 20, 
       height = 7,
       units = "cm",
       dpi = 300)

### 2.7 TPM calculation ###
gene_lenght_fc.cer.bay4 <-cbind(fc.cer.bay4$annotation[1],fc.cer.bay4$annotation[6])
colnames(gene_lenght_fc.cer.bay4) <- c("gene","length")
#add lenght column to raw counts file
raw_counts_cer_bay_total <-merge(raw_counts, gene_lenght_fc.cer.bay4, by="gene",all = T)
write.table(raw_counts_cer_bay_total, "/home/madro269/schraiber2013/index_bowtie2_polyT_cer_bay_conc/results_FC2/polyT_trim_counts_cer_bay.mm.4.length.txt", sep="\t")

#TPM normalization
raw_counts_total_cer_bay_total4 <- read.table("/home/madro269/schraiber2013/index_bowtie2_polyT_cer_bay_conc/results_FC2/polyT_trim_counts_cer_bay.mm.4.length.txt", header = TRUE)
raw_counts_tot_cer_bay_TPM<-as.data.frame(raw_counts_total_cer_bay_total4)

for(i in colnames(raw_counts_tot_cer_par_TPM[2:7])){
  raw_counts_tot_cer_bay_TPM[,i]<-raw_counts_tot_cer_bay_TPM[,i]/raw_counts_tot_cer_bay_TPM[,'length']
  raw_counts_tot_cer_bay_TPM[,i]<-raw_counts_tot_cer_bay_TPM[,i]/sum(raw_counts_tot_cer_par_TPM[,i], na.rm=TRUE)*1e6
}
write.table(raw_counts_tot_cer_par_TPM, "/home/madro269/schraiber2013/index_bowtie2_polyT_cer_bay_conc/results_FC2/raw_counts_tot_cer_bay_TPM.txt", sep="\t")
