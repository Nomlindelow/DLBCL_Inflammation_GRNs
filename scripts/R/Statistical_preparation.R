genie3_genes <- c("MTF1", "MNDA", "IGF2R", "DLL1", "CPA6", "CHMP4C", "GRIN3A", "LIPM", "ANKRD22", "MYOF", "CASP5", "CD163",
"IRG1", "STXBP6", "MEFV", "RP11-77H9.8", "MT1H", "CCL8","PSTPIP2", "FFAR2", "FPR2", "LILRA1", "SIGLEC1", "KCNE1",
"MB")

wgcna_genes <- c("MTF1", "CXCL2", "PDE8B", "DLL1", "GRIN3A", "CASP5", "CD163", "IRG1", "RNASE2",
"MEFV", "RP11-77H9.8", "CPD", "CCL8", "KCNJ2", "PSTPIP2", "FFAR2", "LILRA1", "SIGLEC1")

c25_genes <- c('C11orf53', 'KCNQ1OT1', 'TSPAN1', 'FFAR2', 'BARX2', 'CASP5', 'CACNA1E', 'RP11-830F9.7', 
                    'PTBP1P', 'CXCL2', 'BCL2L10', 'RASGRP4', 'MYOCD', 'ADAMTS6', 'LINC01415', 'VWCE', 'CCL8',
                    'RP4-781K5.6', 'CABLES1', 'MEFV', 'CTSE', 'SPTBN2', 'TMEM132B', 'SHANK1', 'SOX7')

anno <- read.csv("../../project/new/data/GeneAnnotation.csv", sep = ',', header = TRUE)
ens <- anno$ensembl_Id[anno$genename %in% c25_genes]

vst <- read.csv("./files/s4s1/all_vst_counts.csv", sep = ',', header = TRUE)
df <- vst[vst$X %in% ens,]
write.csv(df, "./files/c25/Inf_vst_s4.csv", row.names = FALSE)

library(dplyr)
#df <- read.csv("./files/c25/Inf_vst_s2.csv", sep=",", header = TRUE)
s1 <- df %>% select(matches("^stage4")) %>%
    unlist() %>% data.frame()

write.csv(s1, "./files/c25/s4.csv")
