#! /bin/sh

rm -rf FUCHIKOMA
R CMD BATCH createPackage.R log.txt

# ##### Editing ######
# echo "import (NOISeq)
# import (snow)
# import (Rcpp)

# export(
# 	meta.readData,
# 	meta.oneside.noiseq,
# 	Fisher.test,
# 	Stouffer.test,
# 	other.oneside.pvalues,
# 	Accelerate.NOISeq,
# 	Reset.Accelerate.NOISeq
# )" > ./FUCHIKOMA/NAMESPACE

rm ./FUCHIKOMA/Read-and-delete-me

##### Data #####
cp ????.Rdata FUCHIKOMA/data

##### DESCRIPTION #####
cp DESCRIPTION FUCHIKOMA

##### vignettes #####
mkdir ./FUCHIKOMA/vignettes
mkdir ./FUCHIKOMA/inst
mkdir ./FUCHIKOMA/inst/doc
cp ./vignettes/FUCHIKOMA.pdf ./FUCHIKOMA/inst/doc
cp ./vignettes/FUCHIKOMA.Rnw ./FUCHIKOMA/inst/doc
# cp ./vignettes/Fig1.jpeg ./FUCHIKOMA/inst/doc
# cp ./vignettes/Fig2.png ./FUCHIKOMA/inst/doc

cp ./vignettes/FUCHIKOMA.Rnw ./FUCHIKOMA/vignettes
# cp ./vignettes/Fig1.jpeg ./FUCHIKOMA/vignettes
# cp ./vignettes/Fig2.png ./FUCHIKOMA/vignettes

##### テストコード #####
cp -rf tests/ ./FUCHIKOMA/inst

##### マニュアル #####
rm -rf ./FUCHIKOMA/man
cp -rf man/ ./FUCHIKOMA/man

# Then edit man/.Rd and inst/doc/.Rnw (vignettes/.Rnw) by hand ...

R CMD BUILD --resave-data FUCHIKOMA
R CMD INSTALL FUCHIKOMA_0.0.1.tar.gz
# R CMD CHECK FUCHIKOMA_0.0.1.tar.gz

rm log.txt
