# prepare folders
path <- '/home/hocine/my-projects/nomli/project/new'
setwd(path)

# import libraries 
library(DESeq2)
library(magrittr)
library(GENIE3)
library(ggplot2)
library(edgeR)

data_dir <- "data"
if(!dir.exists(data_dir)) dir.create(data_dir)
plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir)
results_dir <- "results"
if(!dir.exists(results_dir)) dir.create(results_dir)

#import data files
deg_file <- file.path(data_dir, "s1s4-degs.csv")
data_file <- file.path(data_dir, "s1s4-rawcounts.csv")
metada_file <- file.path(data_dir, "s1s4-pData.tsv")
file.exists(deg_file)
file.exists(data_file)
file.exists(metada_file)

raw_counts <- readr::read_csv(data_file) %>%
  tibble::column_to_rownames("transcriptID")
degs <- readr::read_csv(deg_file) 
metadata <- readr::read_tsv(metada_file) %>%
  tibble::column_to_rownames("id")

ens <- sub("\\.\\d+$", "", degs$transcriptID) 

group <- factor(c(rep("S1",8), rep("S4",12)))
dge_list <- DGEList(counts = raw_counts, group = group)
keep <- filterByExpr(dge_list)
summary(keep)
dge_list <- dge_list[keep, , keep.lib.sizes=FALSE]
counts <- dge_list$counts %>%
  as.data.frame()

all.equal(rownames(metadata), colnames(counts))
#meta_cols <- data.frame(labelDescription=c("Condition Name", "Sample ID"), row.names = c("conditionName", "sampleId"))
#pData <- new("AnnotatedDataFrame", data = metadata, varMetadata=meta_cols)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~1
)

# normalization
dds_norm <- varianceStabilizingTransformation(dds)
#dds_norm <- vst(dds)
normalized_counts <- assay(dds_norm)[,9:20] %>%
  as.data.frame() %>%
  tibble::rownames_to_column('transcriptID') %>%
  dplyr::mutate(transcriptID = sub("\\.\\d+$", "", transcriptID)) %>%
  tibble::column_to_rownames("transcriptID") %>%
  as.matrix()

# run genie3. Change the target based on the compared stages
set.seed(123)
res <- GENIE3(normalized_counts, regulators = NULL, targets = ens, nCores=1)
w <- getLinkList(res, threshold = 0.00211) # (s2) 0.0015 (s4) 0.00211 
readr::write_tsv(w, file.path("results", "s1s4-genie3.tsv"))
dim(w)
