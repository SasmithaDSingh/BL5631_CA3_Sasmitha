  #Installing all the necessary packages
  BiocManager::install("DESeq2")
  BiocManager::install("biomaRt")
  BiocManager::install("EnhancedVolcano")
  install.packages("pheatmap")
  
  #Loading the packages
  library(DESeq2)
  library(biomaRt)
  library(EnhancedVolcano)
  library(pheatmap)
  library(here)
  
  #setting the working directory using the "here" command
  here()
  
  #Reading the metadata from the excel file containing sample details
  metadata <- read.csv("samples.csv")
  rownames(metadata) <- metadata$sampleID
  
  #Creating a list object of the different sample files
  Samples <- lapply(metadata$filename,
              read.table,
              header=TRUE)
  names(Samples) <- metadata$sampleID
  class(Samples)
  length(Samples)
  colnames(Samples[[1]])
  
  #Plotting for effective length Vs length of Day 0 data
  par(mfrow=c(1,2)) 
  plot(y=Samples$Day0$effective_length,
       x=Samples$Day0$length, xlab="Day0 length", ylab="Day0 effective length")
  plot(y=log10(Samples$Day0$effective_length+1),
       x=log10(Samples$Day0$length+1), xlab="log10(Day0 length)", ylab="log10(Day0 effective length)")
  lines(y=1:5, x=1:5, col="red")
  
  
  # Combining all the expected counts of different samples
  countdS <- cbind(Samples[[1]]$expected_count,
                   Samples[[2]]$expected_count,
                   Samples[[3]]$expected_count,
                   Samples[[4]]$expected_count,
                   Samples[[5]]$expected_count,
                   Samples[[6]]$expected_count,
                   Samples[[7]]$expected_count)
  rownames(countdS) <- Samples[[1]]$gene_id
  colnames(countdS) <- names(Samples)
  
  #Creating matrix array and filtering out the zero values
  countmat <- as.matrix(countdS)
  countmat.int <- apply(countmat, 2, function(x) round(x))
  dim(countmat.int)
  barplot(table(rowSums(countmat.int>0)))
  countmat.int2 <- countmat.int[apply(countmat.int,1, function(x) !any(x==0)),]
  dim(countmat.int2)
  
  #Conducting differential expression using DESeq2
  dds <- DESeqDataSetFromMatrix(countData = countmat.int2,
    colData = metadata,
    design = ~ treatment)
  dds <- DESeq(dds)
  res <- results(dds)
  
  head(counts(dds))
  head(counts(dds,normalized=TRUE))
  colData(dds)$sizeFactor
  
  #Finding the Gmean by considering only the non-zero values
  find.gmean <- function(myvec) {
    gmean <- exp(mean(log(myvec[myvec>0])))
    if (all(myvec==0)) {
      return(0)
    } else {
      return(gmean)
    }
  }
  
  gene.gmean <- apply(countmat.int2,
                      1,
                      find.gmean)
  
  #Determining the gene ratios and median ratios (using the size factors command)
  gene.ratio <- t(sapply(1:nrow(countmat.int2),
                         function(x) 
                           countmat.int2[x,]/gene.gmean[x]))
  head(gene.ratio)
  sizeFactors(dds)
  
  #Plotting sizefactor/median ratios vs the total read counts
  par(mfrow=c(1,1)) 
  plot(sizeFactors(dds), colSums(counts(dds)), xlab="Median ratios or size factors(dds)", ylab="Total read counts(dds)")
  
  #Plotting the dispersion estimates
  plotDispEsts(dds)
  
  #Examining the read counts for a single gene across the groups and cross-checking for normalization
  counts(dds, normalized=FALSE)[5,]/sizeFactors(dds) == counts(dds, normalized=TRUE)[5,]
  plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")
  
  #Comparing plots of non-normalized and normalized read counts
  par(mfrow=c(1,2)) 
  boxplot(log2(counts(dds)+1), notch=TRUE,
          main = "Non-normalized read counts",
          ylab="log2(read counts)", cex = 0.75)
  boxplot(log2(counts(dds, normalize= TRUE) +1), notch=TRUE,
          main = "Size-factor-normalized read counts",
          ylab="log2(read counts)", cex=0.75)
  
  #Checking for the number of up-regulated and down regulated genes 
  #Adjusting the p value for significance according our requirement and organizing according to their values
  res.sig <- res[(res$padj<0.01)&(!is.na(res$padj)),]
  res.sig.up <- res.sig[res.sig$log2FoldChange > 1,]
  res.sig.down <- res.sig[res.sig$log2FoldChange < (-1),]
  res.sig.up <- res.sig.up[order(res.sig.up$padj, decreasing=FALSE),]
  res.sig.down <- res.sig.down[order(res.sig.down$padj, decreasing=FALSE),]
  summary(res.sig)
  summary(res.sig.up)
  summary(res.sig.down)
  
  #Deriving rlog transformed, log2 transformed and normalized data sets for further visualization
  log.norm.counts <- log2(counts(dds, normalized=TRUE))
  rlog.norm.counts <- assay(rlog(dds))
  
  #Plotting Day1 Vs Day0 transformed data (Similarly the rest of the days can also be plotted)
  plot(log.norm.counts[,1:2], cex=0.5, main="size factor and log2-transformed")
  plot(rlog.norm.counts[,1:2], cex=0.5, main="rlog transformed")
  
  #Comparing similarity between the different day samples using a PCA plot
  plotPCA(rlog(dds), intgroup="treatment")
  
  #Bocplot representation of the top 20 up-regulated and down regulated genes
  boxplot(rlog.norm.counts[rownames(res.sig.up)[1:20],], main="Top 20 upregulated genes")
  boxplot(rlog.norm.counts[rownames(res.sig.down)[1:20],], main="Top 20 downregulated genes")
  
  #Loading the derived genes into a matrix format
  norm.mat <- rlog.norm.counts[c(rownames(res.sig.up)[1:20],
                                 rownames(res.sig.down)[1:20]),]
  upregulated <- norm.mat[1:20,]
  downregulated <- norm.mat[21:40,]
  
  #Generating a MA-plot/scatter plot of log2 fold changes Vs the mean of normalized counts 
  par(mfrow=c(1,1))
  plotMA(res)
  
  #Generating a basic Volcano plot
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = "log2FoldChange",
                  y="padj")
  
  #Annotating the Ensembl gene IDs to the gene symbols for norm.mat as well as upregulated and downregulated
  diffgenes <- data.frame(ensembl_gene_id_version=rownames(norm.mat),
                          ensembl_gene_id=NA)
  diffgenes$ensembl_gene_id <- sapply(diffgenes$ensembl_gene_id_version,
                                      function(x)
                                        strsplit(x,"[.]")[[1]][1])
  mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
  GID <- getBM(c("ensembl_gene_id","hgnc_symbol"),
             "ensembl_gene_id",
             diffgenes$ensembl_gene_id,
             mart)
  diffgenes.symbol <- merge(diffgenes, GID)
  symbols <- diffgenes.symbol[match(rownames(norm.mat),
                                    diffgenes.symbol$ensembl_gene_id_version),"hgnc_symbol"]
  symbols1 <- diffgenes.symbol[match(rownames(upregulated),
                                     diffgenes.symbol$ensembl_gene_id_version),"hgnc_symbol"]
  symbols2 <- diffgenes.symbol[match(rownames(downregulated),
                                     diffgenes.symbol$ensembl_gene_id_version),"hgnc_symbol"]
  rownames(norm.mat) <- symbols
  rownames(upregulated) <- symbols1
  rownames(downregulated) <- symbols2 
  
  #Generating pheatmaps
  pheatmap(norm.mat, scale="row")
  pheatmap(upregulated, scale = "row", main="Day 0 -> 6 upregulated genes")
  pheatmap(downregulated, scale = "row", main="Day 0 -> 6 down regulated genes")
