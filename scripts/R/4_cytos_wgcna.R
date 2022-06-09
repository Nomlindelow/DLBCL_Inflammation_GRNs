setwd('/home/hocine/my-projects/nomli/project/new/')

library(WGCNA)

options(stringsAsFactors = FALSE);
lnames = load(file = "DLBCL-dataInput.RData");
lnames = load(file = "DLBCL-networkConstruction.RData");

deg_file = read.csv(file = "./data/s1s4-degs.csv")
degs <- deg_file$ensembl
annot = read.csv(file = "./data/GeneAnnotation.csv");

TOM = TOMsimilarityFromExpr(datExpr, power = 7, nThreads = 6);
modules = c("brown");
probes = colnames(datExpr)
#probes <- sub("\\.\\d+$","",probes)
inModule = is.finite(match(moduleColors, modules));
inDegs = probes %in% degs
modProbes = probes[inDegs & inModule]
modGenes = annot$genename[match(modProbes, annot$ensembl_Id)];
modTOM = TOM[inDegs & inModule, inDegs & inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inDegs & inModule]);
