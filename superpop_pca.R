# ----- PARAMETERS -----
start_bp <- 41500000
end_bp   <- 43500000

# ----- SNP SUBSET (locus window) -----
ind.locus <- which(map$chromosome == 22 &
                     map$physical.pos >= start_bp &
                     map$physical.pos <= end_bp)

# MAF filter
maf <- snp_MAF(G, ind.col = ind.locus)
ind.locus2 <- ind.locus[maf >= 0.05 & maf <= 0.95]

cat("SNPs after MAF:", length(ind.locus2), "\n")

# ----- PCA ON INDIVIDUALS -----
X <- as.matrix(G[, ind.locus2])

# drop zero-variance SNPs
keep <- apply(X, 2, var) > 0
X <- X[, keep, drop = FALSE]

# scale
X <- scale(X)

# PCA
pca2 <- prcomp(X, center = FALSE, scale. = FALSE)

# ----- MERGE LABELS -----
pca2_df <- data.frame(
  sample.ID = fam$sample.ID,
  PC1 = pca2$x[,1],
  PC2 = pca2$x[,2],
  stringsAsFactors = FALSE
)

pca2_df <- merge(pca2_df, fam2, by = "sample.ID", all.x = TRUE)

# ----- PLOT (colored by superpop) -----
sp <- factor(pca2_df$Super.Population)

plot(pca2_df$PC1, pca2_df$PC2,
     col = sp, pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = paste0("chr22:", start_bp, "-", end_bp, " (CYP2D6 locus region)"))

legend("topright", legend = levels(sp),
       col = seq_along(levels(sp)), pch = 19, cex = 0.8)
