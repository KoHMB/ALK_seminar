########################################################################
########################################################################
## 
## 体長ビンと成長モデルがALK推定精度に与える影響評価シミュレーション
## 
########################################################################
########################################################################
library(tidyverse)
library(magrittr)
library(FSA)


# パラメータ設定
## 基本、某資源瀬戸内系群を想定種にしている

## 資源年齢に対しての漁業の選択性selectivity
FAA <- c(0.31, 0.26, 0.39, 0.64, 0.64) # 2010年のVPAのFAA参考
selectivity <- prop.table(FAA)

## 成長モデル（von Bertalanffy）
L_inf_true <- 560
K_true <- 0.6
a0_true <- -0.144
VB <- function(x, L_inf=L_inf_true, K=K_true, a0=a0_true)L_inf*(1-exp(-K*(x-a0)))
### 成長モデルの不確実性
sd_VB <- 0.5


##-------------------------------------------------------------------##
#  年齢別漁獲尾数の生成
## Catch_total <- 10000


## 超幾何分布で乱数生成するのは、また後日
## ここは値を固定して指定したほうがいい(割合はCAA準拠)

CAA_true <- c(5000,2500,2500,2300,2000) # 2016年のCAA
age_vec_true <- list()
for(age in 1:length(CAA_true))age_vec_true[[age]] <- rep(age-1, CAA_true[age])
age_vec_true <- unlist(age_vec_true)
paCAA_true <- prop.table(CAA_true)


##-------------------------------------------------------------------##
#  サンプリング

## length frequencyデータの抽出
### サンプリング尾数
n_sample1 <- 300
age_sample1_true <- sample(age_vec_true, n_sample1)

### 体長の乱数生成 from VB(場所ここでいいのかな?)
length_sample1 <- purrr::map(age_sample1_true, 
                             function(x){VB(x)*exp(rnorm(1, -.5*sd_VB^2, sd_VB))}
                             ) %>% unlist()
LF_mat_true <- data.frame(id = 1:n_sample1,
                          age = age_sample1_true,
                          length = length_sample1
                          )
LF_mat_true %<>% mutate(bin_L = lencat(length, w=50))


## age-lengthアンプルの抽出
### サンプリング尾数
n_sample2 <- 100
sample2_label <- sample(1:n_sample1, n_sample2) %>% sort()
AL_mat <- data.frame(id = sample2_label,
                     age = LF_mat_true[sample2_label,]$age,
                     bin_L = LF_mat_true[sample2_label,]$bin_L)






### ALKの計算（FSAパッケージ利用）
#ALK_freq <- xtabs(~bin_L + age, data = AL_mat)
#ALK_est <- prop.table(ALK_freq, margin = 1)
ALK_freq <- xtabs(~bin_L + age, data = LF_mat_true)
ALK_est <- prop.table(ALK_freq, margin = 1)
alkPlot(ALK_est, type = "area", pal = "gray", showLegend = TRUE,
        leg.cex = .7, xlab = "Total Length(mm)")


length_n <- xtabs(~bin_L, data = LF_mat_true)
alkAgeDist(ALK_est, lenA.n = rowSums(ALK_freq), len.n = length_n)
## length freqとALKのビンの構成が合っていないとダメ？？





# ビン数の最適化
bin_sim <- c(10,100)
#bin_sim <- seq(10,100,10)
res <- list()
for(i in 1:length(bin_sim)){
  LF_mat_true <- data.frame(id = 1:n_sample1,
                            age = age_sample1_true,
                            length = length_sample1
  )
  LF_mat_true %<>% mutate(bin_L = lencat(length, w=bin_sim[i]))
  ALK_freq <- xtabs(~bin_L + age, data = LF_mat_true)
  ALK_est <- prop.table(ALK_freq, margin = 1)
  length_n <- xtabs(~bin_L, data = LF_mat_true)
  res[[i]] <- alkAgeDist(ALK_est, lenA.n = rowSums(ALK_freq), len.n = length_n)
}


a <- numeric()
for(i in 1:length(res))a[i] <- res[[i]]$prop[1]


