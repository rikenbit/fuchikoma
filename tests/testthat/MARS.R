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
expect_equivalent(names(result.pca.MARS), c("sdev", "rotation", "center", "scale", "x"))
expect_is(result.pca.MARS, "prcomp")

# Result of DiffusionMap against MARS (353 * 228, matrix)
data(result.destiny.MARS)
expect_is(result.pca.MARS, "DiffusionMap")
