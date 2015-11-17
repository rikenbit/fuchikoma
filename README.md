# FUCHIKOMA
FUCHIKOMA（ふちこま、天斑駒が由来）は、BAHSICとDiffusion Mapを利用した、非線形なPCAにおける主成分に寄与した遺伝子を特定する手法です。

# 特徴・メリット
1. Diffusion Mapを利用した非線形なPCA
2. 各主成分に寄与する遺伝子をBAHSICによる特徴選択で特定する
3. 複数の主成分が組み合わさった方向に寄与する遺伝子も選択できる
4. あらかじめクラス数の判定をする必要がない（教師無し）
5. ***による並列化
6. single-cell RNA-Seqに固有な***の補正
7. ***による$p$-value、その後のFDR補正の対応

# 使い方
```r
input <- read.table("input.txt")
result <- FUCHIKOMA(input)
summary(result)
```
