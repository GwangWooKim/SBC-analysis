library("DESeq2")
library(glue)
library(ggplot2)
library("vsn")
library("pheatmap")
library(plyr)
library("RColorBrewer")

cts <- as.matrix(merged_count_matrix_3[,c(2:138)])
rownames(cts) <- merged_count_matrix_3[[1]]

coldata <- data.frame("cell_type" = condition_3$condition)
rownames(coldata) <- colnames(cts)
coldata[sapply(coldata, is.character)] <- lapply(coldata[sapply(coldata, is.character)], as.factor)

coldata$V1 <- condition_3$V1
coldata$V2 <- condition_3$V2

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ V1 + V2 + cell_type)
                              
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

## pairwise inference

cell <- "CCC" # possible choices in [CCC, CRC, CRT]

# default analysis
res <- results(dds, contrast=c("cell_type", glue("{cell}"), 'BNEW'), alpha = 0.05)
summary(res)

# shrinkage method
resLFC <- lfcShrink(dds, coef=glue("cell_type_{cell}_vs_BNEW"), res = res)
summary(resLFC)

# result save
resLFCOrdered <- resLFC[order(resLFC$padj),]
write.csv(as.data.frame(resLFCOrdered), 
          file=glue("resLFC_cell_type_{cell}_vs_BNEW.csv"))
          
# plotting to compare two results
plotMA(res, ylim=c(-2,2), main=glue("res_{cell}_vs_BNEW"))
plotMA(resLFC, ylim=c(-2,2), main=glue("resLFC_{cell}_vs_BNEW"))

# other two options of shrinkage
resNorm <- lfcShrink(dds, glue("cell_type_{cell}_vs_BNEW"), type="normal", res = res)
summary(resNorm)
resAsh <- lfcShrink(dds, glue("cell_type_{cell}_vs_BNEW"), type="ashr", res = res)
summary(resAsh)

# comparisons to each other
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

# report filtering result 1
plot(metadata(resLFC)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(resLFC)$lo.fit, col="red")
abline(v=metadata(resLFC)$filterTheta)

# report filtering result 2
use <- resLFC$baseMean > metadata(resLFC)$filterThreshold
h1 <- hist(resLFC$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resLFC$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "resLFC", ylab="frequency", xlab="p-value")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

# report genes with the lowest p-value
# par(mfrow=c(1,4), mar=c(4,4,2,1))
# 
# for (i in 1:3) {
#   num = which(rownames(dds) == rownames(resOrdered[resOrdered$log2FoldChange<0,])[i])
#   plotCounts(dds, gene=num, intgroup="cell_type")
# }
# 
# num = which(rownames(dds) == "Gm4204")
# plotCounts(dds, gene=num, intgroup="cell_type")

# report genes with the lowest p-value
par(mfrow=c(1,4), mar=c(4,4,2,1))
lst <- c("Sparcl1", "Rufy3", "Pfn2", "Gm4204")

for (i in 1:4) {
  d <- plotCounts(dds, gene=lst[i], intgroup="cell_type", 
                  returnData=TRUE)
  d$cell_type <- as.factor(d$cell_type)
  p <- ggplot(data=d, aes(x=cell_type, y=count)) + 
    geom_boxplot()
  p
}

d <- plotCounts(dds, gene=lst[3], intgroup="cell_type", 
                returnData=TRUE)
d$cell_type <- as.factor(d$cell_type)
ggplot(data=d, aes(x=cell_type, y=count)) + 
  scale_y_log10(breaks=c(1e0, 1e1, 1e2, 1e3, 1e4)) +
  ylab("normalized count") +
  ggtitle(glue("{lst[3]}")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_boxplot()

# mean-variance trend 
vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- data.frame("cell_type" = condition_3$condition_2)
rownames(df) <- colnames(cts)
# df <- arrange(df, cell_type)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE, show_colnames = FALSE,
         cluster_cols=FALSE, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
        
# distance matrix among samples
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell_type, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
        
# 2d visualization (note: there is a batch effect)
# plotPCA(vsd, intgroup=c("cell_type"))

# perhaps, the following code provides an adjusted 2d visualization
# mat <- assay(vsd)
# mm <- model.matrix(~condition[[2]], colData(vsd))
# mat <- limma::removeBatchEffect(mat, batch=, design=mm)
# assay(vsd) <- mat
# plotPCA(vsd, intgroup=c("cell_type"))
