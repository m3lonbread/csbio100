# ----- SUBSET TO AFR + EAS INDIVIDUALS -----
pops_keep <- c("AFR", "EAS")
sub_ids <- fam2$sample.ID[fam2$Super.Population %in% pops_keep]
ind_rows <- which(fam$sample.ID %in% sub_ids)

# ----- PCA RECOMPUTED ON ONLY AFR+EAS INDIVIDUALS -----
# reuse the SAME SNP set ind.locus2 from Script 2 (important)
X_sub <- as.matrix(G[ind_rows, ind.locus2])

keep <- apply(X_sub, 2, var) > 0
X_sub <- X_sub[, keep, drop = FALSE]
X_sub <- scale(X_sub)

pca_sub <- prcomp(X_sub, center = FALSE, scale. = FALSE)

# dataframe for individuals in AFR+EAS
pca_sub_df <- data.frame(
  sample.ID = fam$sample.ID[ind_rows],
  PC1 = pca_sub$x[,1],
  PC2 = pca_sub$x[,2],
  stringsAsFactors = FALSE
)

pca_sub_df <- merge(pca_sub_df, fam2, by = "sample.ID")

# ----- POPULATION CENTROIDS (GBR/FIN/etc — but here only AFR+EAS pops exist) -----
centroids_pop <- aggregate(cbind(PC1, PC2) ~ Population,
                           data = pca_sub_df,
                           mean)

# add superpopulation label for each Population
centroids_pop <- merge(centroids_pop,
                       unique(fam2[, c("Population", "Super.Population")]),
                       by = "Population")

# ----- PLOT CENTROIDS COLORED BY SUPERPOP -----
cols <- ifelse(centroids_pop$Super.Population == "AFR", "black", "red")

plot(centroids_pop$PC1, centroids_pop$PC2,
     col = cols, pch = 19, cex = 2,
     xlab = "PC1", ylab = "PC2",
     main = "Population-Level PCA (AFR vs EAS)")

text(centroids_pop$PC1, centroids_pop$PC2,
     labels = centroids_pop$Population,
     pos = 3, cex = 0.9)

legend("topright", legend = c("AFR", "EAS"),
       col = c("black", "red"), pch = 19)

# ----- MOST SEPARATED AFR vs EAS POPULATION PAIR (centroid distance) -----
afr_cent <- subset(centroids_pop, Super.Population == "AFR")
eas_cent <- subset(centroids_pop, Super.Population == "EAS")

dist_matrix <- matrix(NA, nrow = nrow(afr_cent), ncol = nrow(eas_cent))
rownames(dist_matrix) <- afr_cent$Population
colnames(dist_matrix) <- eas_cent$Population

for (i in 1:nrow(afr_cent)) {
  for (j in 1:nrow(eas_cent)) {
    dist_matrix[i, j] <- sqrt(
      (afr_cent$PC1[i] - eas_cent$PC1[j])^2 +
        (afr_cent$PC2[i] - eas_cent$PC2[j])^2
    )
  }
}

max_idx <- which(dist_matrix == max(dist_matrix), arr.ind = TRUE)
afr_pop <- rownames(dist_matrix)[max_idx[1]]
eas_pop <- colnames(dist_matrix)[max_idx[2]]

cat("Most separated AFR–EAS pair:", afr_pop, "vs", eas_pop, "\n")
cat("Distance =", max(dist_matrix), "\n")


