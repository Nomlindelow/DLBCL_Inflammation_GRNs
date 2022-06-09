library(data.table)
library(edgeR)
library(DESeq2)

# Raw-counts
data=read.table("../data/s1s2-rawcounts.csv", row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

ensembl <- row.names(data)
ensembl <- sub("\\.\\d+$","",ensembl)
row.names(data) <- ensembl

group <- factor(c(rep("S1",8), rep("S2",17)))
dge <- DGEList(counts = data, group = group)
keep <- filterByExpr(dge)
summary(keep)
dge <- dge[keep, , keep.lib.sizes=FALSE]
genes <- rownames(dge$counts)

# pData
pData <- read.csv("../data/s1s2-pData.tsv",row.names = 1, sep = ",", header = T)
all(rownames(pData)==colnames(m))
metadata <- data.frame(labelDescription=c("Condition Name"), row.names = c("sampleId"))
phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata=metadata)
sampleId <- phenoData$sampleId
eset <- ExpressionSet(assayData = m, phenoData = phenoData)
colData <- DataFrame(pData(eset))
se <- SummarizedExperiment(exprs(eset), colData=colData)

mode(assay(se)) <- "integer"

dds <- DESeqDataSet(se, design = ~ sampleId)
dds <- DESeq(dds)

res <- results(dds, alpha = 0.05, altHypothesis = "greaterAbs")
summary(res)
resultsNames(dds)
mcols(res, use.names = T)

res <- res[order(res$padj),]
resSig <- res[res$padj < 0.05 & !is.na(res$padj), ]
summary(resSig)
write.csv(res, "../outputs/Stage1and2_DEGs.csv")


