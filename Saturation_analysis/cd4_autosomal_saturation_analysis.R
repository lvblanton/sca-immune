#Adapted from Adrianna San Roman and Alex Godfrey
#This is a script to determine whether the analysis 
#has saturated for detection of autosomal genes that respond 
#to Chr X or Chr Y dosage in CD4+ T cells. 

inputValue <- commandArgs(trailingOnly = TRUE)

SAMPLE.SIZE <- inputValue
ITERATIONS <- 100

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("BiocParallel"))
register(MulticoreParam(10))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library("limma"))

#read in the metadata table and restrict to the appropriate samples
metadata <- read.delim("/sca-immune/Autosomal_linear_regressions/cd4/CD4_metadata.csv", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$SampleID
metadata.cd4$genotype <- as.factor(metadata.cd4$Karyotype)

#bring in expressed genes
load(file="/sca-immune/Autosomal_linear_regressions/cd4/cd4_expressed.rda")
expressedGenes_cd4_autosomal <- cd4_expressed$expAutoGenes

#bring in counts
load(file="/sca-immune/Autosomal_linear_regressions/cd4/cd4.txi.rda")

#create results table
saturationResults <- NULL
# resultSampleTables <- vector("list")

set.seed(seed=1)

#for each sample size in the table, choose that many samples and iterate over the 
for(size in c(1:length(SAMPLE.SIZE))){
  # size <- 1
  print(paste0("Sample size: ",SAMPLE.SIZE[size]))
  i <- SAMPLE.SIZE[size]
  
  #do the samplings 
  j <- 1
  repeat{
    print(paste0("Round: ",j))
    
    #randomly select samples
    mySamples <- sample(1:dim(metadata.cd4)[1], i, replace = FALSE)
    print(mySamples)
    #select subset of samples
    cd4.txi.random <- cd4.txi
    cd4.txi.random$counts <- cd4.txi.random$counts[,mySamples]
    cd4.txi.random$abundance <- cd4.txi.random$abundance[,mySamples]
    cd4.txi.random$length <- cd4.txi.random$length[,mySamples]
    
    cd4.txi.random$counts <- cd4.txi.random$counts[,order(colnames(cd4.txi.random$counts))]
    cd4.txi.random$abundance <- cd4.txi.random$abundance[,order(colnames(cd4.txi.random$abundance))]
    cd4.txi.random$length <- cd4.txi.random$length[,order(colnames(cd4.txi.random$length))]

    cd4.metadata.random <- metadata.cd4[metadata.cd4$SampleID %in% colnames(cd4.txi.random$counts),]
    cd4.metadata.random <- cd4.metadata.random[order(rownames(cd4.metadata.random)),]
    
    #Test whether # of samples > coefficients
    num_distinct <- length(unique(cd4.metadata.random$GTC_Run_ID))
    i != (num_distinct - 1 + 3)
    if (i != (num_distinct - 1 + 3)) {
      
      #Test whether matrix is full rank
      mod <- model.matrix( ~ GTC_Run_ID + Y_count + X_count, cd4.metadata.random)
      is.fullrank(mod)
      if(is.fullrank(mod)){
        #If the matrix is full rank
        
        #Perform DESEq
        mydds <- DESeqDataSetFromTximport(cd4.txi.random, colData = cd4.metadata.random, design = ~ GTC_Run_ID + Y_count + X_count)
        mydds <- DESeq(mydds, parallel = TRUE)
        
        ##Do analysis of X genes##
        x_res <- na.omit(results(mydds, name = "X_count", alpha = 0.05))
        x_res_exp_auto_sig <- x_res[rownames(x_res) %in% expressedGenes_cd4_autosomal & x_res$padj < 0.05,]
        
        
        #Do analysis of Y genes
        y_res <- na.omit(results(mydds, name = "Y_count", alpha = 0.05))
        y_res_exp_auto_sig <- y_res[rownames(y_res) %in% expressedGenes_cd4_autosomal & y_res$padj < 0.05,]
        
        
        #Record the number of significant X or Y responsive genes in a table
        sigGenes <- c(i, dim(x_res_exp_auto_sig)[1], dim(y_res_exp_auto_sig)[1])
        names(sigGenes) <- c("Sample_num","x_auto_sig","y_auto_sig")
        saturationResults <- rbind(saturationResults, c(sigGenes , summary(cd4.metadata.random$genotype)))
        
        # resultSampleTables[paste0("round_",as.character(i),"_",as.character(j))] <- list(lcl.metadata.random)
        
        #print the results out
        print(paste0("Finished round ",as.character(i),",",as.character(j),". # of X responsive genes: ", as.character(dim(x_res_exp_auto_sig)[1]),"; # of Y-responsive genes: ",as.character(dim(y_res_exp_auto_sig)[1])))
        
        #advance j
        j <- j + 1
        
      } else{
        print("Model not full rank, drawing samples again.")
      }
      
    } else{
      print("Number of samples equals number of coefficients, drawing samples again.")
    }
    
    if(j == ITERATIONS + 1){
      break
    }
    
  }
  
  ##for each sample size
  saturationResults <- as.data.frame(saturationResults)
  write.table(saturationResults, file=paste0("/sca-immune/Saturation_analysis/cd4_autosomal_output/cd4_saturationResults_100iterations_", as.character(i),".txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
  
}

sessionInfo()