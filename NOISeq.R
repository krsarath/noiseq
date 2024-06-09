library(BiocManager)
library(NOISeq)
library(ggplot2)
library(VennDiagram)
library(ggVennDiagram)
setwd("/home/sarath/R/rwd")

mycounts <- read.csv("matrixcounts.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
head(mycounts)

mycounts <- mycounts[, c(1:4)]
head(mycounts)

myfactors <- data.frame(group = c("control","treated","control","treated"))
head(myfactors)

mydata <- readData(data=mycounts, factors=myfactors)
head(mydata)

exprs <- assayData(mydata)$exprs
exprs <- apply(exprs, 2, function(x) {
  if (any(is.na(x))) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
  }
  return(x)
})
myTMM <- tmm(exprs, long = 1000, lc = 0)
head(myTMM)

myfilt = filtered.data(mycounts, factor = myfactors, norm = T, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")

mydata2 <- NOISeq::readData(data=myfilt, factors=myfactors)

myresults <- noiseq(mydata2, factor = "group", k = 0, norm = "tmm", condition = c("control","treatment"), pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "technical")
myresults.deg = degenes(myresults, q = 0.8, M = NULL)
myresults.deg1 = degenes(myresults, q = 0.8, M = "up")
myresults.deg2 = degenes(myresults, q = 0.8, M = "down")

up <- rownames(myresults.deg1)
down <- rownames(myresults.deg2)
all <- rownames(exprs)
no_expression_genes <- setdiff(all_genes, union(up_genes, down_genes))
gene_sets <- list(
  UP = up,
  DOWN = down,
  ALL = no_expression_genes
)
ggVennDiagram(gene_sets, label_alpha = 0)

N <- 20 
top_genes <- myresults.deg[order(abs(myresults.deg$ranking), decreasing = TRUE), ]
top_genes <- head(top_genes, N)
top_gene_names <- rownames(top_genes)
top_gene_exprs <- exprs[top_gene_names, ]
heatmaply(top_gene_exprs, 
          main = "Heatmap of Top Significant Genes",
          xlab = "Samples", 
          ylab = "Genes",
          scale = "row",  
          colors = viridis::viridis(256))

head(myresults.deg2)
head(myresults.deg1)
head(myresults.deg)
DE.plot(myresults, q = 0.8, graphic = "expr", log.scale = TRUE)
DE.plot(myresults, q = 0.8, graphic = "MD")

save(myresults, file = "NoiseqResults.RData")
write.table(myresults@results[[1]], file="noiseqresults.txt", sep="\t", row.names=T, col.names=T, quote=F)
