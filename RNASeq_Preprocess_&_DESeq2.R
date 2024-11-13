library(DESeq2) # rna-seq
library(stringr)
library(rafalib) # nice plot arrangement
library(dplyr) # data wrangling
library(tibble) # data wrangling
library(tidyr) # data wrangling
library(ggplot2) # plotting
library(edgeR) # rna-seq
library(pheatmap)
library(biomaRt)

rafalib::mypar(mar=c(6,2.5,2.5,1))

# Load data - gene count table
co <- read.delim("gene_counts_original.txt", sep="\t", header=TRUE, stringsAsFactors=F, comment.char="#")


# ***** OUTLIERS REMOVED *****

co<-subset(co, select=-c(Cd0d, Cd2_4c))
# Load metadata table - outliers removed (Cd0d + Cd2_4c)
mo <- read.csv2("metadata_outliers.csv", header=TRUE)
# We need to keep only the counts, gene IDs and gene names/symbol
cr <- co[,1:27]


# ***** ALL DATA *****

# Load metadata table - all data
mo <- read.csv2("metadata.csv",header=TRUE)
# We need to keep only the counts, gene IDs and gene names/symbol
cr <- co[,1:29]

# Also we set row names of count table to the gene IDs, column names to sample name
# and rownames of metadata table to sample name.
rownames(cr) <- cr$gene_id
cr <- subset(cr, select = -gene_id)
colnames(cr) <- substr(colnames(cr),42,43)
colnames(cr) <- mo$SampleName
rownames(mo) <- mo$SampleName

# Check to see if columns in count table are the same as the rows in the metadata table
# should be TRUE if everything worked out
all.equal(colnames(cr),rownames(mo))

# Save prepared data
write.csv(cr,file="output/gene_counts_raw.csv",quote=FALSE)
write.csv(mo,file="output/metadata_raw.csv",quote=FALSE)


# ***** FILTERING *****

# Visualize distribution of counts in raw data using a box plot and density plot.
rafalib::mypar(1,2,mar=c(6,3,3,2))
boxplot(log2(as.matrix(cr)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Raw data")
hist(log2(as.matrix(cr)+1),ylab="",las=2,main="Raw data")
par(mfrow=c(1,1))

# Visualize the number of detected genes per sample
{barplot(colSums(cr>3),main="Number of detected genes",las=2)
  abline(h=median(colSums(cr>3)))}

# Visualize detection rate across genes for the different samples
{barplot(rowSums(cr>3),xlab="Genes",ylab="Number of samples",names.arg="")
  abline(h=median(rowSums(cr>3)),col="red")}
# Or this:
hist(rowSums(cr>3), xlab = "detection rate (>3)")

# Remove genes with low counts
# cr > 5 means that genes with fewer than 5 reads will be removed
# >= 3 means that the gene must've been detected in at least 3 replicates
keep_genes <- rowSums( cr > 5 ) >= 3
cf <- cr[keep_genes,]


# Visualize the new distribution of counts (filtered)
boxplot(log2(as.matrix(cf)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Filtered data")
# and the detection rate of genes for the different samples
hist(rowSums(cf>3), xlab = "detection rate (>3)")

# No samples should've been discarded and the number of rows-columns should remain the same
# should be TRUE
all.equal(colnames(cf),rownames(mo))

# Save the filtered data
write.csv(cf,"output/counts_filtered.csv",quote=F)


# ***** PCA and Clustering *****

# Normalization of data using TPM, for PCA and clustering only
tpm <- function(counts,len) {
  x <- counts/(len/1000)
  return(t(t(x)*1e6/colSums(x)))
}

# read data
metadata <- read.csv("output/metadata_raw.csv",header=TRUE,stringsAsFactors=TRUE,row.names=1)
co <- read.delim("gene_counts_original.txt",
                 sep="\t",
                 header=TRUE,
                 stringsAsFactors=F,
                 comment.char="#")
co<-subset(co, select=-c(Cd0d, Cd2_4c))
g <- data.frame( ensembl_gene_id = co$gene_id , 
                 transcript_length = co$gene_length,
                 stringsAsFactors = F, row.names = co$gene_id)
g <- g[!duplicated(g$ensembl_gene_id),]

igenes <- intersect(rownames(cf),g$ensembl_gene_id)
g1 <- g[igenes,]
cf1 <- cf[igenes,]
all.equal(rownames(cf1),g1$ensembl_gene_id)

# TPM normalization
ct <- tpm(cf1,g1$transcript_length)
log2_ct <- log2( ct + 1 )
boxplot(log2_ct,ylab=expression('Log'[2]~'TPM'),las=2,main="Log2 TPM")
write.csv(ct,"output/counts_tpm.csv",quote=F)

#read data for further analysis
data <- read.csv("output/counts_tpm.csv",header=TRUE,stringsAsFactors=FALSE,row.names=1)
data <- log2( data + 1 )

#convert ensembl ZF gene IDs to hgnc symbols
mart1 <- useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl", host = "https://dec2021.archive.ensembl.org")
mart2 <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

ids <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id",
              values = rownames(data), mart = mart1,
              attributesL = c("ensembl_gene_id", "hgnc_symbol"), martL = mart2, )
colnames(ids) <- c("gene_id", "human_gene_id", "gene_symbol")
data <- rownames_to_column(data, "gene_id")
data <- left_join(data, subset(ids, select = -human_gene_id), by = "gene_id")
data <- data[!(data$gene_symbol == ""),]
data <- na.omit(data)
data <- data[!duplicated(data[,"gene_symbol"]),]
rownames(data) <- NULL
data <- column_to_rownames(data, "gene_symbol")
data <- subset(data, select = -gene_id)
write.csv(data ,"output/gene_tpm_counts.csv",row.names=T)

cf <- read.csv("output/counts_filtered.csv",header=TRUE,stringsAsFactors=FALSE,row.names=1)
ids <- getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id",
              values = rownames(cf), mart = mart1,
              attributesL = c("ensembl_gene_id", "hgnc_symbol"), martL = mart2, )
colnames(ids) <- c("gene_id", "human_gene_id", "gene_symbol")

cf <- rownames_to_column(cf, "gene_id")
cf <- left_join(cf, subset(ids, select = -human_gene_id), by = "gene_id")
cf <- cf[!(cf$gene_symbol == ""),]
cf <- na.omit(cf)
cf <- cf[!duplicated(cf[,"gene_symbol"]),]
id_and_name <- subset(cf, select = c(gene_id, gene_symbol))
rownames(cf) <- NULL
cf <- column_to_rownames(cf, "gene_symbol")
cf <- subset(cf, select = -gene_id)
write.csv(cf ,"output/gene_counts.csv",row.names=T)

# Perform z-score normalization for PCA and clustering analysis
Znorm <- t(apply(data,1,function(x) scale(x,center=T,scale=T)))
colnames(Znorm) <- colnames(data)

# raw data vs normalized data example of top 30 genes
{
  mypar(1,2,mar=c(5,3,2,1))
  boxplot(t(data[1:30,]),ylim=c(0,12),las=2,col="grey",main="data raw_data",cex=.2)
  boxplot(t(Znorm[1:30,]),ylim=c(-4,4),las=2,col="grey",main="Z-score data",cex=.2)
  abline(h=0,col="red",lty=2)
}

#Compute the mean, variance and cv for each gene and sort them in decreasing order
gene_stats <- data.frame(row.names=rownames(data))
means <- apply(data,1,mean)
vars <- apply(data,1,var)

#Plot the row means versus variance
{
  mypar(1,2,mar=c(5,3,2,1))
  plot(means,vars,cex=.1)
}

#Sort and select the top 500 highly variable genes from the data
vars <- sort(vars,decreasing=T)
top_var <- names(vars)[1:500]
boxplot(t(data[top_var[1:15],]),ylim=c(0,12),las=2,col="grey", main="data (top var genes)",cex=.2)

# Plot PCA, PC1 vs PC2 and PC3 vs PC4
{
  mypar(1,2,mar = c(3,3,2,1))
  #PC <-  prcomp( t( Znorm[ top_var, ]) ) #Method1
  PC <-  prcomp( t( data[ top_var, ]), center = TRUE, scale. = TRUE) #Method2
  
  mypar(1,2)
  plot(PC$x[,1],PC$x[,2],cex=2,col=factor(metadata$Group),xlab="PC1",ylab="PC2",pch=16,main="PCA",las=1)
  text(PC$x[,1],PC$x[,2],cex=.7,labels = paste0(metadata$SampleName),pos=3)
  
  plot(PC$x[,3],PC$x[,4],cex=2,col=factor(metadata$Group),xlab="PC3",ylab="PC4",pch=16,main="PCA",las=1)
  text(PC$x[,3],PC$x[,4],cex=.7,labels = paste0(metadata$SampleName),pos=3)
}

# compute and plot PC variance
PC_sd <- setNames(PC$sdev,paste0("PC",1:length(PC$sdev)))
PC_var_expl <- (PC_sd^2)/sum(PC_sd^2)*100

{
  mypar(1,1)
  barplot(PC_var_expl,las=2,ylab="% variance explained")
  abline(h=10,lty=2)
}

leading_genes <- PC$rotation
head(leading_genes)

leading_PC1 <- sort(leading_genes[,1],decreasing=T)
leading_PC2 <- sort(leading_genes[,2],decreasing=T)

#plot leading genes in PC1 and PC2
{
  mypar(1,2,mar=c(3,4,2,2))
  barplot(leading_PC1[15:1],las=2,horiz=T,cex.names=.8,cex.axis=0.8,yaxs="i")
  abline(v=0,lwd=2)
  barplot(leading_PC2[15:1],las=2,horiz=T,cex.names=.8,cex.axis=0.8,yaxs="i")
  abline(v=0,lwd=2)
}

d <- dist( t(data) , method="euclidean")
d

#Compute sample correlations
sample_cor <- cor( data )
round(sample_cor,4)
pheatmap(sample_cor)

#Transform the scale from correlations
cor_distance <- -(sample_cor-1)/2
round(cor_distance,4)
pheatmap(cor_distance)

#Convert it to a distance object
d2 <- as.dist(cor_distance)
d2

#Clustering using euclidean distance and correlation
{
  mypar(1,2,mar=c(6,4,2,1))
  h <- hclust(d,method="complete")
  plot( as.dendrogram(h),las=1,main="d=euclidean\nh=complete")
  points(1:ncol(data),rep(0,ncol(data)),pch=16,cex=2,col=metadata$Group[h$order])
}

h2 <- hclust(d2,method="complete")
{
  plot( as.dendrogram(h2),las=1, main="d=correlation\nh=complete")
  points(1:ncol(data),rep(0,ncol(data)),pch=16,cex=2, col=metadata$Group[h2$order])
}

gene_cor  <- cor(t(Znorm[top_var, ]))
gene_dist <- as.dist(-(gene_cor-1)/2)
gene_clus <- hclust(gene_dist,method="complete")

HEIGHT <- 0.9
gene_clusters <- cutree(gene_clus,h=HEIGHT)
gene_clusters

{
  mypar(1,1,mar=c(6,4,2,1))
  plot( as.dendrogram(gene_clus),las=1,main="d=correlation\nh=complete")
  
  rect.hclust(gene_clus,h=HEIGHT)
  abline(h=HEIGHT,col="red",lty=2)
  points(1:length(gene_clusters),rep(0,length(gene_clusters)),pch=16,cex=2, col=factor(gene_clusters)[gene_clus$order])
  legend("topright",levels(factor(gene_clusters)),pch=16,col=  factor(levels(factor(gene_clusters))))
}

pheatmap( data[top_var,] , scale="row" , color = colorRampPalette(c("navy","white","firebrick"))(90),
          border_color=NA, cluster_cols=F)



# ***** Differential Gene Expression analysis *****

# Data is converted into a DESeq2 object.
# The "design" parameter is the generalized linear model, and "~Group" means
# we only compare by that variable.
mo$Group <- factor(mo$Group)
d <- DESeqDataSetFromMatrix(countData=cf, colData=mo, design=~Group)

# Estimate normalization factors.The median-of-ratios method is used.
d <- DESeq2::estimateSizeFactors(d,type="ratio")
gene_counts_norm <- counts(d, normalized=TRUE)
saveRDS(gene_counts_norm, "output/outliers_removed_QC/DESeq2_normalized_gene_counts.rds")

# Step 2. Gene dispersion
# Measures of variability like variance and standard deviation are not fit for RNASeq data.
# Instead, dispersion is used, as it works well for negative bionomially distributed data.
# DESeq2 estsimates the dispersion for each gene based on counts and fits a curve through the estimates.
# Then, shrinks the dispersion towards the fitted curve.
d <- DESeq2::estimateDispersions(d)
# Dispersions visualized:
plotDispEsts(d)
# Black points are the maximum dispersion estimate for each gene.
# Blue points are the new gene dispersion estimates after shrinkage toward the curve.
# Blue circled black points are estimates that were not shrunk, too far away from the curve.

#To investigate additional groups:
dg <- DESeq(d)

# Repeat for other comparisons and insert data into new table
res <- id_and_name

t <- results(dg, contrast=c("Group", "PCB0_08", "Cd0"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".PCB0_08_vs_Cd0"))

t <- results(dg, contrast=c("Group", "PCB0_4", "Cd0"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".PCB0_4_vs_Cd0"))

t <- results(dg, contrast=c("Group", "PCB0_8", "Cd0"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".PCB0_8_vs_Cd0"))

t <- results(dg, contrast=c("Group", "Cd2_4", "Cd0"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".Cd2_4_vs_Cd0"))

t <- results(dg, contrast=c("Group", "Cd12", "Cd0"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".Cd12_vs_Cd0"))

t <- results(dg, contrast=c("Group", "Cd24", "Cd0"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".Cd24_vs_Cd0"))

t <- results(dg, contrast=c("Group", "Cd24", "PCB0_8"), alpha=0.05)
t <- as.data.frame(t)
t <- rownames_to_column(t, var="gene_symbol")
t$padj[is.na(t$padj)] <- 1
res <- left_join(res, t, by="gene_symbol", suffix=c("", ".Cd24_vs_PCB0_8"))
colnames(res)[2:8] <- c("gene_symbol", "baseMean.PCB0_08_vs_Cd0",
                        "log2FoldChange.PCB0_08_vs_Cd0", "lfcSE.PCB0_08_vs_Cd0",
                        "stat.PCB0_08_vs_Cd0", "pvalue.PCB0_08_vs_Cd0", "padj.PCB0_08_vs_Cd0")

#res <- read.csv("output/outliers_removed_QC/dge_results.csv")
A <- res[res$padj.PCB0_08_vs_Cd0 < 0.05 ,2]
B <- res[res$padj.PCB0_4_vs_Cd0 < 0.05 ,2]
C <- res[res$padj.PCB0_8_vs_Cd0 < 0.05 ,2]
D <- res[res$padj.Cd2_4_vs_Cd0 < 0.05 ,2]
E <- res[res$padj.Cd12_vs_Cd0 < 0.05 ,2]
F <- res[res$padj.Cd24_vs_Cd0 < 0.05 ,2]

# Generate venn diagrams with the number of DEGs
VennDiagram::venn.diagram(list("PCB 0.08" = A, "PCB 0.4" = B, "PCB 0.8" = C), filename = "Venn_PCB_padj0.05_DEGS.png", disable.logging = TRUE, 
                          imagetype = "png", fontfamily = "sans", cat.fontfamily = "sans", fill = c("dodgerblue2", "orchid3", "seagreen3"),
                          cex = 1.5, cat.cex = 1.25, alpha = 0.4, resolution = 300, width = 1600, height = 1600, na = "remove")

VennDiagram::venn.diagram(list("Cd 2.4" = D, "Cd 12" = E, "Cd 24" = F), filename = "Venn_Cd_padj0.05_DEGS.png", disable.logging = TRUE, 
                          imagetype = "png", fontfamily = "sans", cat.fontfamily = "sans", fill = c("gold3", "darkorange2", "red3"),
                          cex = 1.5, cat.cex = 1.25, alpha = 0.4, resolution = 300, width = 1600, height = 1600, na = "remove")


# Save results as a CSV file
write.csv(res,"output/dge_results.csv",row.names=F)


# Some useful visualizations for the DEG output:

# To generate a PCA plot to visualize the principal components
dgvst<-varianceStabilizingTransformation(dg)
plotPCA(DESeqTransform(dgvst), intgroup="Group")

#MA plot (mean expression vs log FC for all genes)
DESeq2::plotMA(dg)

# Volcano plot (log FC vs padj)
# Grey line denotes the significance threshold line, smallest p-values are at the top.
ggplot()+
  geom_point(data=as.data.frame(res),aes(x=log2FoldChange.PCB0_08_vs_Cd0,y=-log10(padj.PCB0_08_vs_Cd0)),col="grey80",alpha=0.5)+
  geom_point(data=filter(as.data.frame(res),padj.PCB0_08_vs_Cd0<0.05),aes(x=log2FoldChange.PCB0_08_vs_Cd0,y=-log10(padj.PCB0_08_vs_Cd0)),col="red",alpha=0.7)+
  geom_hline(aes(yintercept=-log10(0.05)),alpha=0.5)+
  coord_cartesian(ylim=c(0, 20))+
  theme_bw()

ggplot()+
  geom_point(data=as.data.frame(res),aes(x=log2FoldChange.PCB0_4_vs_Cd0,y=-log10(padj.PCB0_4_vs_Cd0)),col="grey80",alpha=0.5)+
  geom_point(data=filter(as.data.frame(res),padj.PCB0_4_vs_Cd0<0.05),aes(x=log2FoldChange.PCB0_4_vs_Cd0,y=-log10(padj.PCB0_4_vs_Cd0)),col="red",alpha=0.7)+
  geom_hline(aes(yintercept=-log10(0.05)),alpha=0.5)+
  theme_bw()

ggplot()+
  geom_point(data=as.data.frame(res),aes(x=log2FoldChange.PCB0_8_vs_Cd0,y=-log10(padj.PCB0_8_vs_Cd0)),col="grey80",alpha=0.5)+
  geom_point(data=filter(as.data.frame(res),padj.PCB0_8_vs_Cd0<0.05),aes(x=log2FoldChange.PCB0_8_vs_Cd0,y=-log10(padj.PCB0_8_vs_Cd0)),col="red",alpha=0.7)+
  geom_hline(aes(yintercept=-log10(0.05)),alpha=0.5)+
  theme_bw()

ggplot()+
  geom_point(data=as.data.frame(res),aes(x=log2FoldChange.Cd2_4_vs_Cd0,y=-log10(padj.Cd2_4_vs_Cd0)),col="grey80",alpha=0.5)+
  geom_point(data=filter(as.data.frame(res),padj.Cd2_4_vs_Cd0<0.05),aes(x=log2FoldChange.Cd2_4_vs_Cd0,y=-log10(padj.Cd2_4_vs_Cd0)),col="red",alpha=0.7)+
  geom_hline(aes(yintercept=-log10(0.05)),alpha=0.5)+
  theme_bw()

ggplot()+
  geom_point(data=as.data.frame(res),aes(x=log2FoldChange.Cd12_vs_Cd0,y=-log10(padj.Cd12_vs_Cd0)),col="grey80",alpha=0.5)+
  geom_point(data=filter(as.data.frame(res),padj.Cd12_vs_Cd0<0.05),aes(x=log2FoldChange.Cd12_vs_Cd0,y=-log10(padj.Cd12_vs_Cd0)),col="red",alpha=0.7)+
  geom_hline(aes(yintercept=-log10(0.05)),alpha=0.5)+
  theme_bw()

ggplot()+
  geom_point(data=as.data.frame(res),aes(x=log2FoldChange.Cd24_vs_Cd0,y=-log10(padj.Cd24_vs_Cd0)),col="grey80",alpha=0.5)+
  geom_point(data=filter(as.data.frame(res),padj.Cd24_vs_Cd0<0.05),aes(x=log2FoldChange.Cd24_vs_Cd0,y=-log10(padj.Cd24_vs_Cd0)),col="red",alpha=0.7)+
  geom_hline(aes(yintercept=-log10(0.05)),alpha=0.5)+
  theme_bw()

