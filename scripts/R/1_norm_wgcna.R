# prepare folders
path <- '/home/hocine/my-projects/nomli/project/new/'
setwd(path)

# import libraries 
library(DESeq2)
library(magrittr)
library(WGCNA)
library(ggplot2)
library(edgeR)

#import data files
data_dir <- "data"
deg_file <- file.path(data_dir, "s1s4-degs.csv")
data_file <- file.path(data_dir, "s1s4-rawcounts.csv")
metada_file <- file.path(data_dir, "s1s4-pData.tsv")
anno_file <- file.path(data_dir, "stages_pheno.csv")
file.exists(deg_file)
file.exists(data_file)
file.exists(metada_file)
file.exists(anno_file)

raw_counts <- readr::read_csv(data_file) %>%
  tibble::column_to_rownames("transcriptID")
degs <- readr::read_csv(deg_file) 
pData <- readr::read_tsv(metada_file) %>%
  tibble::column_to_rownames("id")
traitData = readr::read_csv(anno_file)

group <- factor(c(rep("S1",8), rep("S4",12)))
dge_list <- DGEList(counts = raw_counts, group = group)
keep <- filterByExpr(dge_list)
summary(keep)
dge_list <- dge_list[keep, , keep.lib.sizes=FALSE]
counts <- dge_list$counts %>%
  as.data.frame()

all.equal(rownames(pData), colnames(counts))
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = pData,
  design = ~1
)

# normalization
dds_norm <- varianceStabilizingTransformation(dds)

# Uncomment in case you want to save the normalized matrix in a file.
#norm_matrix <- assay(dds_norm) %>% 
#  as.data.frame() %>% 
#  tibble::rownames_to_column('transcriptID') %>%
#  dplyr::mutate(transcriptID = sub("\\.\\d+$", "", transcriptID))
#readr::write_csv(norm_matrix, file = file.path("results", "normalized_s1s2.tsv"))

#datExpr <- as.data.frame(t(subset(assay(dds_norm), rownames(assay(dds_norm)) %in% ens)))
datExpr <- assay(dds_norm) %>%
  as.data.frame() %>%
  tibble::rownames_to_column('transcriptID') %>%
  dplyr::mutate(transcriptID = sub("\\.\\d+$", "", transcriptID)) %>%
  tibble::column_to_rownames("transcriptID") %>%
  t()

# Use real sample names instead of stage1x or stage2x
datExpr <- datExpr %>%
  as.data.frame() %>%
  tibble::rownames_to_column('sampleID') %>%
  dplyr::mutate(sampleID = pData$sampleId) %>%
  tibble::column_to_rownames('sampleID') %>%
  as.matrix()

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

samples <- rownames(datExpr)
traitRows <- match(samples, traitData$sample_id)
datTraits <- traitData[traitRows, -1] %>%
  dplyr::mutate(traitData[traitRows, 1]) %>%
  tibble::column_to_rownames('sample_id')

# convert gender to numeric
gender <- character(0)
j <- 1
for(g in datTraits$gender){
  if(g == "male")
    gender[j] <- 0
  else if (g == "female")
    gender[j] <- 1
  j <- j + 1
}
datTraits$gender <- as.numeric(gender)

# convert clinical stage to numeric
stage <- character(0)
j <- 1
for(s in datTraits$clinical_stage){
  if(s == "Stage I")
    stage[j] <- 1
  else if(s == "Stage II")
    stage[j] <- 2
  else if(s == "Stage III")
    stage[j] <- 3
  else if(s == "Stage IV")
    stage[j] <- 4
  j <- j + 1
}
datTraits$clinical_stage <- as.numeric(stage)

collectGarbage()

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "DLBCL-dataInput.RData")
