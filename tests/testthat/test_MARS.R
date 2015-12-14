# Count data of MARS (353 * 228, matrix)
data(MARS)
expect_equivalent(dim(MARS), c(353, 228))
expect_true(is.matrix(MARS))

# Cell label vector of MARS (353 * 228, matrix)
data(label.MARS)
expect_equivalent(length(label.MARS), 228)
expect_true(is.vector(label.MARS))

# Result of PCA against MARS (353 * 228, matrix)
data(result.pca.MARS)
expect_equivalent(dim(result.pca.MARS), c(228, 6))
expect_true(is.matrix(result.pca.MARS))

# Result of t-SNE against MARS (353 * 228, matrix)
data(result.tsne.MARS)
expect_equivalent(dim(result.tsne.MARS), c(228, 2))
expect_true(is.matrix(result.tsne.MARS))

# Result of DiffusionMap against MARS (353 * 228, matrix)
data(result.dmap.MARS)
expect_equivalent(dim(result.dmap.MARS), c(228, 6))
expect_true(is.matrix(result.dmap.MARS))
