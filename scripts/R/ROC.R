library(OptimalCutpoints)
ROCData <- read.csv("s1s2ROC.csv", sep=",", header = TRUE)
ROCData <- factor(ROCData$ labels=c("Absent", "Present"))
results<- optimal.cutpoints(X= CCL8,
                            status="status",
                            tag.healthy="Absent",
                            methods="Youden", data=CCL8)
summary(results)S