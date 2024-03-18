# Script for running power analyses on RNA-seq data
# Originally written by Alex Godfrey and Updated by Adrianna San Roman
# Adapted by Laura Blanton to run on CD4+ T cell and monocyte data, 11/04/22
# Bsub and run in terminal...takes several hours.

message("Starting analysis")
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("parallel"))

#Bring in functions
source(file = "/sca-immune/Rfunctions/Function_geneAnno2anno_20221103.R")
source(file = "/sca-immune/Rfunctions/Function_powerAnalysis_20221103.R")
source(file="/sca-immune/Rfunctions/Function_calculatePower_20221103.R")
message("loaded packages and functions")

#### Files needed - adjust to your files #### 
#Gene annotation file that matches your mapped transcriptome - needs "Gene", "chr", "start" columns:
geneAnno <- read.delim("/sca-immune/Annotations/geneAnno_proteinCoding_lncRNA_v107_20221021.txt", stringsAsFactors = FALSE)

#Set output directory for files
OUT.DIR <- "/sca-immune/Power_analysis/cd4_output"  #your output directory

#Load tximport files
# load and pre-process data
load("/sca-immune/Autosomal_linear_regressions/cd4/cd4.txi.rda")
data.cd4 <- cd4.txi$counts


#Load metadata
meta.cd4 <- read.csv(file = "CD4_metadata.csv", sep = "\t", row.names = "SampleID")
names(meta.cd4)[names(meta.cd4) == 'X_count'] <- 'x_count'
names(meta.cd4)[names(meta.cd4) == 'Y_count'] <- 'y_count'
names(meta.cd4)[names(meta.cd4) == 'GTC_Run_ID'] <- 'batch_libprep'
message("loaded metadata")

#Load expressed genes
load("/sca-immune/Autosomal_linear_regressions/cd4/cd4_expressed.rda")
expgenes.cd4 <- cd4_expressed$expressedGenes
message("loaded data")


#### Prep files and metadata ####
#Get gene lists from annotation file 
#Main "Gene" name will be ensembl84
colnames(geneAnno)[c(2,4:6)] <- c("Gene","chr","start","stop")
#Only include genes in both v84 and v107
geneAnno <- geneAnno[! is.na(geneAnno$Gene),]
geneAnno <- geneAnno[! duplicated(geneAnno$gene_name.107) & ! geneAnno$gene_name.107 %in% c("") & ! duplicated(geneAnno$Gene),]

geneAnno_x <- geneAnno[geneAnno$chr == "chrX",]
X_genes_all <-geneAnno_x$Gene
Y_genes_all <- geneAnno[geneAnno$chr == "chrY", "Gene"]
sexChrom_genes <- c(X_genes_all,Y_genes_all)
PAR1_genes <- geneAnno_x[geneAnno_x$start < 2691188,]$Gene
PAR2_genes <- c("SPRY3","VAMP7","IL9R","WASH6P")
PAR_genes_all <- c(PAR1_genes,PAR2_genes)
NPX_genes <- X_genes_all[! X_genes_all %in% PAR_genes_all]
NPY_genes <- Y_genes_all[! Y_genes_all %in% PAR_genes_all]
autosome_genes <- geneAnno[! geneAnno$chr %in% c("chrX","chrY"),]
anno <- na.omit(geneAnno2anno(geneAnno))


#### Do power analysis calculations ####
numCores <- detectCores()

N.REPS <- 20
N.GENES <- 100
BETAS <- c(0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

# CD4s
message("working on cd4s...")
res.cd4 <- powerAnalysis(data.cd4, meta.cd4, betas=BETAS, n.genes=N.GENES, 
                         n.reps=N.REPS, anno=anno, expgenes=expgenes.cd4) 
write.table(res.cd4, file=file.path(OUT.DIR, "power_results.cd4.txt"),
            quote=FALSE, sep="\t")
