# install once
#install.packages("bigsnpr")
library(bigsnpr)
library(bigstatsr)

# downloads these into data_1000G/
# 1000G_phase3_common_norel.bed/.bim/.fam/.fam2
bigsnpr::download_1000G("data_1000G")

bedfile <- "data_1000G/1000G_phase3_common_norel.bed"

# IMPORTANT:
rds <-  snp_readBed(bedfile)
obj <- snp_attach(rds)

G   <- obj$genotypes
map <- obj$map
fam <- obj$fam

# read population labels
fam2 <- read.table("data_1000G/1000G_phase3_common_norel.fam2",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# quick sanity checks
stopifnot(all(c("sample.ID","Population","Super.Population") %in% names(fam2)))
