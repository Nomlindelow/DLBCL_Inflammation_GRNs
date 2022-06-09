# packages
library(DescTools)  
library(dplyr)
library(car)
library(agricolae)

# Set working directory
setwd("/home/hocine/my-projects/nomli/thesis/stats/")

# Read csv data
df <- read.csv("files/c25/stats_c25.csv", sep=",", header = TRUE)

# Convert Grp variable to factor
df$Grp <- factor(df$Grp)

#### Levene test ####
lv <- leveneTest(Exp ~ Grp , data = df)

# Print results
lv

# Extract p value
lv[1,3]

#### ANOVA ####
av <- aov(Exp ~ Grp, data = df)

# Print results
summary(av)

#Extract p-val
summary(Anova_S1vs2)[[1]][1,5]

## No need for a posthoc analysis

#### Kruskal.test ####
kw <- kruskal.test(Exp ~ Grp, data = df)

# Print results
kw

#### No need for a posthoc test
