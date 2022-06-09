setwd('/home/hocine/my-projects/nomli/project/new/')

library(WGCNA)

load(file = "DLBCL-dataInput.RData");
load(file = "DLBCL-networkConstruction.RData")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

clinical_stage = as.data.frame(datTraits$clinical_stage)
names(clinical_stage) = "stage"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
geneTraitSignificance = as.data.frame(cor(datExpr, clinical_stage, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(clinical_stage), sep = "")
names(GSPvalue) = paste("p.GS.", names(clinical_stage), sep = "")

module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for clinical stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


ensembl <- names(datExpr)[moduleColors=="blue"]
annot = read.csv(file = "./data/GeneAnnotation.csv")
names(annot)
probes = colnames(datExpr)
#probes <- sub("\\.\\d+$","",probes)
probes2annot = match(probes, annot$ensembl.Id)
sum(is.na(probes2annot))

geneInfo0 = data.frame(ensembl = probes,
                       geneSymbol = annot$genename[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
modOrder = order(-abs(cor(MEs, clinical_stage, use = "p")))
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                         MMPvalue[,modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                       paste("p.MM", modNames[modOrder[mod]], sep = ""))
}

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.stage))
geneInfo = geneInfo0[geneOrder,]
write.csv(geneInfo, file = "./results/wgcna/geneInfo-S4-N.csv")
