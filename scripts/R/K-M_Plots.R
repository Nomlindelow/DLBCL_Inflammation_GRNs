library(survival)
library(ranger)
library(dplyr)
library(ggfortify)
library(survminer)

data <- read.table("../test-data/km-os.csv", row.names=1, header=TRUE, sep=",", stringsAsFactors=FALSE)

data$time <- as.numeric(data$time)

data <- mutate(data, status = ifelse((score < 2.16), 0, 1), status = as.numeric(status)) 
#data$status <- as.numeric(data$status)

fit <- survfit(Surv(time, status) ~ strata(Risk), data = data)

#autoplot(fit, conf.int = F)

ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           palette = c("#E7B800", "#2E9FDF"),data = data)

