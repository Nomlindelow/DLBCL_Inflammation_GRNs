setwd('/home/hocine/my-projects/nomli/project/new/')

library(WGCNA)

#enableWGCNAThreads()
lnames = load(file = "DLBCL-dataInput.RData")

powers = c(c(1:10))
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft <- pickSoftThreshold(datExpr,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "unsigned")

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1.0;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

bwnet <- blockwiseModules(datExpr,
                          maxBlockSize = 4000,
                          TOMType = "unsigned",
                          power = 7,
                          numericLabels = TRUE,
                          randomSeed = 1234)

table(bwnet$colors)
mergedColors = labels2colors(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "DLBCL-networkConstruction.RData")
