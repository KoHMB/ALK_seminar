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
#library(DLMtool)


# パラメータ設定
## 基本、某資源瀬戸内系群を想定種にしている

## 資源年齢に対しての漁業の選択性selectivity
FAA <- c(0.31, 0.26, 0.39, 0.64, 0.64) # 2010年のVPAのFAA参考
selectivity <- prop.table(FAA)

## 成長モデル（von Bertalanffy）
L_inf_true <- 560
K_true <- 0.6
a0_true <- -0.144
VB <- function(x, L_inf=L_inf_true, K=K_true, a0=a0_true, 
               a_sd=0.05, b_sd=0.1, deterministic = FALSE){
  if(isTRUE(deterministic)){
    L_inf*(1-exp(-K*(x-a0)))
  } else {
    sd <- sqrt(a_sd + b_sd*x)
    L_inf*(1-exp(-K*(x-a0)))*exp(rnorm(1, -.5*sd^2, sd))
  }
}

a <- numeric()
for(i in 1:1000){a[i] <- VB(5)}
hist(a)

## 最大体長の閾値（これ以上を＋グループにする）
max_length <-700
max_age <-4




##-------------------------------------------------------------------##
#  年齢別漁獲尾数の生成
## Catch_total <- 10000


## 超幾何分布で乱数生成するのは、また後日
## ここは値を固定して指定したほうがいい(割合はCAA準拠)


CAA_true <- c(5000,2500,2500,2300,2000)*5 # 2016年のCAA
age1_dist <- rnorm(CAA_true[1])+1     # 1歳魚の乱数
age2_dist <- rnorm(CAA_true[2])+2     # 2歳魚の乱数
age3_dist <- rnorm(CAA_true[3])+3     # 3歳魚の乱数
age4_dist <- rnorm(CAA_true[4])+4     # 4歳魚の乱数
age5_dist <- rnorm(CAA_true[5])+5     # 5歳魚の乱数



age_vec_true <- c(age1_dist, age2_dist, age3_dist, age4_dist, age5_dist)
age_vec_true <- age_vec_true[-which(age_vec_true <= 0.5)]

CAAprob_true <- floor(age_vec_true) %>% table() %>% prop.table() %>% round(digits = 3)

Sample_true <- data.frame(id = 1:length(age_vec_true), 
                         age = age_vec_true,
                         length = purrr::map(age_vec_true, VB) %>% unlist())

##-------------------------------------------------------------------##
#  Selectivity アイデアないので、とりあえず体長を指定して抽出

#double_norm <- function(x, b1=300, b2=0.02, b3=400, b4=0.02){
#  (1/(1+exp(-b2*(x-b1)))*(1-(1/(1+exp(-b4*(x-b3))))))/1
#}
#curve(double_norm(x), xlim = c(0,700), ylim=c(0,1))


sel_id <- which(Sample_true$length>=300 & Sample_true$length<=600)
random_label <- sample(sel_id, 10000) %>% sort()
Catch_true <- Sample_true[random_label,]

hist(Catch_true$age)
prob_Catch_true <- Catch_true$age %>% floor() %>% table()%>% prop.table()


##-------------------------------------------------------------------##
#  サンプリング

## length frequencyデータの抽出
### サンプリング尾数
n_sample1 <- 200
sample1_true <- sample(1:length(Catch_true$id), n_sample1)

LF_mat_true <- data.frame(id = Catch_true$id[sample1_true],
                          age = Catch_true$age[sample1_true],
                          length = Catch_true$length[sample1_true]
                          )
LF_mat_true %<>% mutate(bin_L = lencat(length, w=20))
LF_mat_true$bin_L[(LF_mat_true$bin_L>=max_length)] <- max_length


## age-lengthサンプルの抽出
### 各体長階級から10尾ずつサンプリング
### 10尾いない場合は全て持ってくる
id_tmp <- age_tmp <- length_tmp <- bin_tmp <- NA
bL_list <- unique(LF_mat_true$bin_L) %>% sort()
for(bL in 1:length(bL_list)){
  subset_data <- LF_mat_true[LF_mat_true$bin_L==bL_list[bL],]
  if(length(subset_data[,1]) >= 10){
    sample2_label <- sample(1:length(subset_data[,1]), 10) %>% sort()
    id_tmp <- c(id_tmp, subset_data$id[sample2_label])    
    age_tmp <- c(age_tmp, floor(subset_data$age[sample2_label]))
    length_tmp <- c(length_tmp, subset_data$length[sample2_label])
    bin_tmp <- c(bin_tmp, subset_data$bin_L[sample2_label])
  } else {
    id_tmp <- c(id_tmp, subset_data$id)    
    age_tmp <- c(age_tmp, floor(subset_data$age))
    length_tmp <- c(length_tmp, subset_data$length)
    bin_tmp <- c(bin_tmp, subset_data$bin_L)
  }
}

plus_tmp <- cbind(age_tmp, rep(NA, length(age_tmp)))
for (i in 2:length(age_tmp)) {
  plus_tmp[i,2] <- if(plus_tmp[i,1]>=max_age)max_age else plus_tmp[i,1]
}
AL_mat <- data.frame(id = na.omit(id_tmp),
                     age = na.omit(plus_tmp[,2]),
                     length = na.omit(length_tmp),
                     bin_L = na.omit(bin_tmp))

table(AL_mat$bin_L)







### ALKの計算（FSAパッケージ利用）
ALK_freq <- xtabs(~bin_L + age, data = AL_mat)
ALK_est <- prop.table(ALK_freq, margin = 1)
alkPlot(ALK_est, type = "area", pal = "gray", showLegend = TRUE,
        leg.cex = .7, xlab = "Total Length(mm)")
alkPlot(ALK_est, type = "bubble", pal = "gray", showLegend = TRUE,
        leg.cex = .7, xlab = "Total Length(mm)")


length_n <- xtabs(~bin_L, data = LF_mat_true)
alkAgeDist(ALK_est, lenA.n = rowSums(ALK_freq), len.n = length_n)
## length freqとALKのビンの構成が合っていないとダメ？？





# ビン数の最適化
bin_sim <- seq(10,100,10)
iteration <- 100
n_sample1 <- 200
res_ALKest <- res_ALKse <- LF_mat_sim <- AL_mat_sim <- list()
est_mat <- se_mat <- matrix(NA, ncol = 5, nrow = iteration)
#res_ALKse <- rep(NA, 5)
 
for(bb in 1:length(bin_sim)){
  for (ite in 1:iteration) {
    sample1_true <- sample(1:length(Catch_true$id), n_sample1)
  
    ## length frequency sample ----------------------------- ##
    LF_mat_true <- data.frame(id = Catch_true$id[sample1_true],
                              age = Catch_true$age[sample1_true],
                              length = Catch_true$length[sample1_true]
    )
    LF_mat_true %<>% mutate(bin_L = lencat(length, w=bin_sim[bb]))
    LF_mat_true$bin_L[(LF_mat_true$bin_L>=max_length)] <- max_length
    
    ## age-length sample ----------------------------- ##
    id_tmp <- age_tmp <- length_tmp <- bin_tmp <- NA
    bL_list <- unique(LF_mat_true$bin_L) %>% sort()
    for(bL in 1:length(bL_list)){
      subset_data <- LF_mat_true[LF_mat_true$bin_L==bL_list[bL],]
      if(length(subset_data[,1]) >= 10){
        sample2_label <- sample(1:length(subset_data[,1]), 10) %>% sort()
        id_tmp <- c(id_tmp, subset_data$id[sample2_label])    
        age_tmp <- c(age_tmp, floor(subset_data$age[sample2_label]))
        length_tmp <- c(length_tmp, subset_data$length[sample2_label])
        bin_tmp <- c(bin_tmp, subset_data$bin_L[sample2_label])
      } else {
        id_tmp <- c(id_tmp, subset_data$id)    
        age_tmp <- c(age_tmp, floor(subset_data$age))
        length_tmp <- c(length_tmp, subset_data$length)
        bin_tmp <- c(bin_tmp, subset_data$bin_L)
      }
    }
    plus_tmp <- cbind(age_tmp, rep(NA, length(age_tmp)))
    for (pp in 2:length(age_tmp)) {
      plus_tmp[pp,2] <- if(plus_tmp[pp,1]>=max_age)max_age else plus_tmp[pp,1]
    }
    AL_mat <- data.frame(id = na.omit(id_tmp),
                         age = na.omit(plus_tmp[,2])+1,
                         length = na.omit(length_tmp),
                         bin_L = na.omit(bin_tmp))
    
    ## calculation ALK
    ALK_freq <- xtabs(~bin_L + age, data = AL_mat)
    ALK_est <- prop.table(ALK_freq, margin = 1)
    length_n <- xtabs(~bin_L, data = LF_mat_true)
    res_ALK <- alkAgeDist(ALK_est, lenA.n = rowSums(ALK_freq), len.n = length_n)
    
    if(!length(res_ALK$prop)==5){
      est_mat[ite,] <- c(rep(0, (5-length(res_ALK$prop))), res_ALK$prop)  ## ここでは０歳魚の選択率が             ##
      se_mat[ite,] <- c(rep(0, (5-length(res_ALK$prop))), res_ALK$se)     ## 明らかに低いから、この項が入っている ##
    } else {
      est_mat[ite,] <- res_ALK$prop
      se_mat[ite,] <- res_ALK$se
    }
    
    ## back-up of samples
#    res_ALKest[[(ite-1)*10+bb]] <- res_ALK
    LF_mat_sim[[(ite-1)*10+bb]] <- LF_mat_true
    AL_mat_sim[[(ite-1)*10+bb]] <- AL_mat
  }#for(ite)
  res_ALKest[[bb]] <- est_mat
  res_ALKse[[bb]] <- se_mat
}#for(bb)

boxplot(res_ALKest[[1]])
points(1:5, prob_Catch_true[1:5], col="red", pch=16, cex=1.6)






