#Adapted from Adrianna San Roman and Alex Godfrey
#This is a script to determine whether the analysis 
#has saturated for detection of sex chromosomal genes 
#that respond to Chr X or Chr Y dosage in CD4+ T cells. 

inputValue <- commandArgs(trailingOnly = TRUE)

SAMPLE.SIZE <- inputValue
ITERATIONS <- 100

suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("BiocParallel"))
register(MulticoreParam(10))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("broom"))

#read in the metadata table and restrict to the appropriate samples
metadata <- read.delim("/sca-immune/Autosomal_linear_regressions/cd4/CD4_metadata.csv", stringsAsFactors = FALSE)
rownames(metadata) <- metadata$SampleID
metadata.cd4$genotype <- as.factor(metadata.cd4$Karyotype)

#bring in expressed genes
load(file="/sca-immune/Autosomal_linear_regressions/cd4/cd4_expressed.rda")
expressedGenes_cd4_sexlinked <- cd4_expressed$expSexChromGenes

#bring in counts
normCounts_expXgenes <- read.delim(file="/sca-immune/Sex_chromosomal_linear_regressions/cd4/CD4_normCounts_expXgenes.txt", stringsAsFactors = FALSE, check.names = FALSE)

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
    normCounts_expXgenes.random <- normCounts_expXgenes[,mySamples]
    normCounts_expXgenes.random <- normCounts_expXgenes.random[,order(colnames(normCounts_expXgenes.random))]

    cd4.metadata.random <- metadata.cd4[metadata.cd4$SampleID %in% colnames(normCounts_expXgenes.random),]
    cd4.metadata.random <- cd4.metadata.random[order(rownames(cd4.metadata.random)),]
    sample_table_par_npx <- cd4.metadata.random[,c("X_count", "Y_count", "GTC_Run_ID")]
    
    #Test whether # of samples > coefficients
    num_distinct <- length(unique(cd4.metadata.random$GTC_Run_ID))
    i != (num_distinct - 1 + 3)
    if (i != (num_distinct - 1 + 3)) {
      
      #Test whether matrix is full rank
      mod <- model.matrix( ~ GTC_Run_ID + Y_count + X_count, cd4.metadata.random)
      is.fullrank(mod)
      if(is.fullrank(mod)){
        #If the matrix is full rank
        myResults_par_npx <- data.frame(intercept=numeric(),intErr = numeric(),xbeta=numeric(), xpval = numeric(), xErr = numeric(),ybeta=numeric(), ypval=numeric(),yErr=numeric(), adj_rsq=numeric(), pval=numeric())
        region_res <- character()
        
        for(gene in rownames(normCounts_expXgenes.random)){
          options(scipen = 0, digits = 3) 
          #pull out data
          myGeneData <- t(normCounts_expXgenes.random[gene,])
          all_Data <- cbind(sample_table_par_npx, myGeneData)
          colnames(all_Data) <- c("X_count", "Y_count", "GTC_Run_ID", "gene_expression")
          #do a linear regression
          myFormula <- formula(gene_expression ~ X_count + Y_count + GTC_Run_ID)
          myFormula_int <- formula(gene_expression ~ X_count + Y_count )
          mylm <- lm(myFormula, data=all_Data)
          mylm_int <- lm(myFormula_int, data=all_Data)
          
          
          #pull out intercept, betas, p-values and add to table
          intercept <- summary(mylm_int)$coefficients["(Intercept)",1]
          interr <- summary(mylm_int)$coefficients["(Intercept)",2]
          xbeta <- summary(mylm)$coefficients["X_count",1]
          xpval <- summary(mylm)$coefficients["X_count",4]
          xerr <- summary(mylm)$coefficients["X_count",2]
          ybeta <- summary(mylm)$coefficients["Y_count",1]
          ypval <- summary(mylm)$coefficients["Y_count",4]
          yerr <- summary(mylm)$coefficients["Y_count",2]
          adj_rsq <- glance(mylm)$adj.r.squared
          pval <- glance(mylm)$p.value
          
          #build results into table  
          myResults_par_npx[gene,] <- c(intercept,interr,xbeta, xpval, xerr,ybeta, ypval, yerr, adj_rsq,pval)

        }
  
        ##Significant genes##
        bh_adj <- p.adjust(myResults_par_npx$xpval, method="BH")
        myResults_par_npx$x_adj_pval <- bh_adj
        x_sig <- myResults_par_npx[myResults_par_npx$x_adj_pval < 0.05,]
        
        
        #Record the number of significant X or Y responsive genes in a table
        sigGenes <- c(i, dim(x_sig)[1])
        names(sigGenes) <- c("Sample_num","x_linked_sig")
        saturationResults <- rbind(saturationResults, c(sigGenes , summary(cd4.metadata.random$genotype)))
        
        # resultSampleTables[paste0("round_",as.character(i),"_",as.character(j))] <- list(lcl.metadata.random)
        
        #print the results out
        print(paste0("Finished round ",as.character(i),",",as.character(j),". # of X responsive genes: ", as.character(dim(x_sig)[1])))
        
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
  write.table(saturationResults, file=paste0("/sca-immune/Saturation_analysis/cd4_sex_chr_output/cd4_saturationResults_100iterations_", as.character(i),".txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
  
}

sessionInfo()