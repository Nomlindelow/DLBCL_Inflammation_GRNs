library("org.Hs.eg.db")
library(magrittr)


data_dir <- "data"
results_dir <- "results"
data_file <- file.path(results_dir, "s1s4-genie3.tsv")
deg_file <- file.path(data_dir, "s1s4-degs.csv")
anno_file <- file.path(data_dir, "GeneAnnotation.csv")

data <- readr::read_tsv(data_file)
degs <- readr::read_csv(deg_file)
anno <- readr::read_csv(anno_file)

ens <- sub("\\.\\d+$", "", degs$transcriptID) 
data <- data[data$regulatoryGene %in% ens & data$targetGene %in% ens, ] %>%
  dplyr::mutate(regulatoryGene = anno$genename[match(regulatoryGene, anno$ensembl_Id)]) %>%
  dplyr::mutate(targetGene = anno$genename[match(targetGene, anno$ensembl_Id)])
head(data)
readr::write_tsv(data, file.path(results_dir, "filtered_genie3-s4.tsv"))

