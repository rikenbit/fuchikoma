#! /bin/sh

rm -rf FUCHIKOMA
R CMD BATCH createPackage.R log.txt

##### Editing ######
echo "import (destiny)
import (foreach)
import (doParallel)

export(
	custom.DiffusionMap,
	CatKernel,
	HSIC,
	FUCHIKOMA
)" > ./FUCHIKOMA/NAMESPACE

rm ./FUCHIKOMA/Read-and-delete-me

##### Data #####
cp MARS.Rdata FUCHIKOMA/data

##### DESCRIPTION #####
cp DESCRIPTION FUCHIKOMA

##### NEWS #####
mkdir ./FUCHIKOMA/inst
cp NEWS ./FUCHIKOMA/inst

# ##### vignettes #####
# mkdir ./FUCHIKOMA/vignettes
# mkdir ./FUCHIKOMA/inst/doc

# cp ./vignettes/FUCHIKOMA.pdf ./FUCHIKOMA/inst/doc
# cp ./vignettes/FUCHIKOMA.Rnw ./FUCHIKOMA/inst/doc

# cp ./vignettes/FUCHIKOMA.pdf ./FUCHIKOMA/vignettes
# cp ./vignettes/FUCHIKOMA.Rnw ./FUCHIKOMA/vignettes

# ##### テストコード #####
# cp -rf tests/ ./FUCHIKOMA

##### マニュアル #####
rm -rf ./FUCHIKOMA/man
cp -rf man/ ./FUCHIKOMA/man

# Then edit man/.Rd and inst/doc/.Rnw (vignettes/.Rnw) by hand ...

R CMD BUILD --resave-data FUCHIKOMA
R CMD INSTALL FUCHIKOMA_0.0.1.tar.gz
R CMD CHECK FUCHIKOMA_0.0.1.tar.gz

rm log.txt
