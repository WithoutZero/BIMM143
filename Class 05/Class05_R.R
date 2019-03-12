#' ---
#'title: "Class05_R.R"
#'author: "Zhuohang_Wu"
#'output: word_document
#'date: "Thu Jan 24 09:32:32 2019"
#' ---

# Class 05 R grapgic intro

# my first boxplot
x <- rnorm(1000,0)
boxplot(x)
summary(x)
hist(x)
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)
plot(weight, pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="AGE(Months)", 
     ylab="WEIGHT(Kg)", main="Baby Weight with Age")
feature <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")

par(mar=c(5,12,4,3))
barplot(feature$Count, horiz=TRUE, xlab="feature count", names.arg = feature$Feature, 
       main= "Feature Count", las=1, xlim = c(0,80000))

phenotype <- read.table("bimm143_05_rstats/up_down_expression.txt", header= TRUE)
table(phenotype$State)
palette(c("green", "red", "blue"))
plot(phenotype$Condition1, phenotype$Condition2, col=phenotype$State, 
     xlab="expresscondition 1", ylab="express condition 2")

# Lets plot expresion vs gene regulation
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")
plot(meth$gene.meth, meth$expression)

dcols <- densCols(meth$gene.meth, meth$expression)

# Plot changing the plot character ('pch') to a solid circle
plot(meth$gene.meth, meth$expression, col = dcols, pch = 20)
# Find the indices of genes with above 0 expresion
inds <- meth$expression > 0

# Plot just these genes
plot(meth$gene.meth[inds], meth$expression[inds])

## Make a desnisty color vector for these genes and plot
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])

plot(meth$gene.meth[inds], meth$expression[inds], col = dcols, pch = 20)

dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                         colramp = colorRampPalette(c("blue2",
                                                      "green",
                                                      "orange",
                                                      "red2")) )

plot(meth$gene.meth[inds], meth$expression[inds], 
     col = dcols.custom, pch = 20)


