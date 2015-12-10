data(label.MARS)

# 2-class
two_class <- label.MARS[union(which(names(label.MARS) == "B cell"), which(names(label.MARS) == "NK_cell"))
]

res.simple2 <- CatKernel(names(two_class), type="simple")
expect_equivalent(dim(res.simple2), c(93, 93))

res.two2 <- CatKernel(names(two_class), type="two")
expect_equivalent(dim(res.two2), c(93, 93))

expect_warning(res.one_vs_rest2 <- CatKernel(names(two_class), type="one_vs_rest"))

expect_warning(res.each2 <- CatKernel(names(two_class), type="each"))

# 4-class
res.simple4 <- CatKernel(names(label.MARS), type="simple")
expect_equivalent(dim(res.simple4), c(228, 228))

expect_warning(res.two4 <- CatKernel(names(label.MARS), type="two"))

res.one_vs_rest4 <- CatKernel(names(label.MARS), type="one_vs_rest")
expect_equivalent(dim(res.one_vs_rest4), c(228, 228))
res.each4 <- CatKernel(names(label.MARS), type="each")
expect_equivalent(dim(res.each4), c(228, 228))
