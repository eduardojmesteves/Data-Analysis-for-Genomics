##########
#Preparation

# install Bioconductor and check version
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::version()

# install Bioconductor packages
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install(c("genefilter","geneplotter"))

# load installed packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(genefilter)
library(geneplotter)

# get help through the documentation
help.start()
?mean
help(mean)
help(package="genefilter")

# inspect objects, classes and methods
library(Biobase)    # load one of the core Bioconductor packages
?ExpressionSet
?"ExpressionSet-class"
methods(class = ExpressionSet)

# inspect the source code for functions and methods
read.csv
plotMA
showMethods("plotMA")
getMethod("plotMA","data.frame")

# vignettes teach you how to use various functions in a package
vignette(package="Biobase")
vignette("ExpressionSetIntroduction")
browseVignettes(package="Biobase")

# report key details about your R session
sessionInfo()


BiocManager::install(c("genefu",
                       "COPDSexualDimorphism",
                       "gwascat",
                       "hgu133a.db",
                       "genomicsclass/tissuesGeneExpression"))
###############
## Section 1 ##
## Overview assessment
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]

#Q2
sum(is.na(sig.gene70$NCBI.gene.symbol))

#Q3
library(dplyr)

sig.gene70 %>%
  filter(Description == "cyclin E2") %>%
  select(NCBI.gene.symbol)

#Q4
?grep

kinase_indexes <- grep("kinase", sig.gene70$Description, ignore.case = TRUE)
num_kinase_genes <- length(kinase_indexes)
print(num_kinase_genes)


## Assessment: Phenotypes ## 
library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)

#Q1
table(expr.meta$GENDER)

#Q2
summary_table <- summary(expr.meta)
summary_table

#Q3
qqnorm(expr.meta$pkyrs, pch = 1, frame = FALSE)
qqline(expr.meta$pkyrs, col = "steelblue", lwd = 2)

# # Load necessary library
# library(ggplot2)
# 
# # Create a QQ plot
# ggplot(expr.meta, aes(sample = pkyrs)) +
#   stat_qq() +
#   stat_qq_line() +
#   ggtitle("QQ Plot of Pack-Years Smoked") +
#   xlab("Theoretical Quantiles") +
#   ylab("Sample Quantiles")

#Q4
boxplot(pkyrs~gender, data=expr.meta)


######
#Assessment: Chromosomes and SNPs
#Q2
library(gwascat)
data(ebicat_2020_04_30)
ebicat_2020_04_30

sort(table(ebicat_2020_04_30$CHR_ID),decreasing=TRUE)

#Q3
library(GenomicRanges)
sort(mcols(ebicat_2020_04_30)[,"DISEASE/TRAIT"])
Diseases <- sort(table(ebicat_2020_04_30$`DISEASE/TRAIT`))
table(ebicat_2020_04_30$`DISEASE/TRAIT`)
tail(names(sort(table(ebicat_2020_04_30$`DISEASE/TRAIT`))), 1)


#####
#Verified assessment: Biology background
#Q3
library(tissuesGeneExpression)
data(tissuesGeneExpression)

head(e[,1:5])
table(tissue)

#Q4
library(dplyr)

gene_expression <- e["209169_at", ]
expression_data <- data.frame(expression = as.numeric(gene_expression), 
                              tissue = tissue)
mean_expression_by_tissue <- expression_data %>%
  group_by(tissue) %>%
  summarize(mean_expression = mean(expression, na.rm = TRUE))

top_tissues <- mean_expression_by_tissue %>%
  arrange(desc(mean_expression)) %>%
  slice(1:2)

print(top_tissues)


#Q6
library(hgu133a.db)    # installed in section overview
symbol = mapIds(hgu133a.db, keys=rownames(e), column="SYMBOL", keytype="PROBEID")

print(symbol["209169_at"])

#Q7
sum(symbol == "GAPDH", na.rm=TRUE)
sum(symbol == "H2AX", na.rm=TRUE)

#Q9
boxplot(as.numeric(e["205436_s_at",])~tissue)

#Q10
IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")
library(rafalib)
mypar(3,2)
sapply(IDs, function(x) boxplot(e[x,]~tissue, las=2))
