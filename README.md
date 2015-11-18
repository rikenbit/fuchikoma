# FUCHIKOMA
FUCHIKOMA（ふちこま、天乃斑駒が由来）は、BAHSICを利用し、非線形なPCAを解析時に主成分に寄与した遺伝子を特定する手法です。

"F"eat"u"re Sele"c"tion method based on BA"H"S"I"C and "K"ernel's Multivariate Analysis（BAHSICとカーネルの多変量解析に基づく特徴量抽出法）の略です。

# 特徴・メリット
1. 各主成分に寄与する遺伝子をBAHSICによる特徴選択で特定する
2. 複数の主成分が組み合わさった方向に寄与する遺伝子も選択できる
3. あらかじめクラス数の判定をする必要がない（教師無し）
4. ***による並列化
5. single-cell RNA-Seqに固有な***の補正
6. ***によるp-value、その後のFDR補正の対応
7. RNA-Seq, DNAマイクロアレイどちらでも利用

# インストールの仕方（予定）
```bash
git clone https://github.com/rikenbit/FUCHIKOMA
R CMD INSTALL FUCHIKOMA
```

# インストールの仕方（予定）
```r
library(devtools)
devtools::install_github("rikenbit/FUCHIKOMA")
```

# 依存パッケージ（予定）
- kernlab（v-X.Y.Z）
- destiny (v-X.Y.Z)

# 使い方（予定）
```r
# destinyを読み込む
library(destiny)

# 遺伝子発現量データ読み込み
data("testdata")

# 細胞周期ラベル読み込み
data("cellcycle")

# オブジェクト化
testdata.obj <- as.ExpressionSet(testdata)

# Diffusion Mapを実行（第10主成分までを見る）
dif <- DiffusionMap(testdata.obj, n.eigs=10)

# 目視でどの主成分に分かれるのかを確認
pairs(eigenvectors(dif), col=cellcycle)

# FUCHIKOMA実行（教師あり、ラベルあり）
result <- FUCHIKOMA(dif, DC=c(1,2))

# FUCHIKOMA実行（教師なし、ラベルなし）
result <- FUCHIKOMA(dif, DC=c(1,2))

################### FUCHIKOMA実行 ###################
# 教師あり（クラスラベルを与える）
result1 <- FUCHIKOMA(data=testdata, mode="Supervised", label=cellcycle, type="each", n.eigs=10)

# 教師なし（指定した主成分を使う）
result2 <- FUCHIKOMA(data=testdata, mode="Unsupervised", Comp=c(1,2), n.eigs=10)

# DEGs
head(result1$DEGs)
head(result2$DEGs)
```
