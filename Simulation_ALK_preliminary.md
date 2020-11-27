-   [目的](#目的)
-   [解析の流れ](#解析の流れ)
-   [成長モデル（von Bertalanffy）](#成長モデルvon-bertalanffy)
    -   [最大ビン幅と年齢+グループの設定](#最大ビン幅と年齢グループの設定)
-   [シミュレーション1回の流れ](#シミュレーション1回の流れ)
    -   [資源集団の生成と体長データの生成](#資源集団の生成と体長データの生成)
    -   [漁獲物の抽出](#漁獲物の抽出)
    -   [Length
        frequencyサンプルの抽出](#length-frequencyサンプルの抽出)
    -   [Age-lengthサブサンプルの抽出](#age-lengthサブサンプルの抽出)
    -   [ALKの計算](#alkの計算)
-   [ビン幅を変えての計算とその結果](#ビン幅を変えての計算とその結果)
    -   [結果](#結果)
-   [メモ](#メモ)

目的
====

-   ALKを得る際に、どのくらいの体長ビン幅が推定上最適なのか知ること
    -   (正直、細かければ細かいほうがいいという辺りだとは思う)
    -   例えばCV10%で推定結果を得るには?みたいにするのがベターかな
-   成長関数の各パラメータの値と最適ビン幅の関係の把握
    -   L\_infの10%体長とかそんな風に得られたら面白いかな

解析の流れ
==========

1.  5万匹くらいの資源の全体を設定

2.  年齢データがあるから、体長に変換する

3.  体長についての選択曲線から漁獲物を抽出(上手いアイデアなかったのでとりあえず300-500mmをランダム抽出)

4.  1万匹の漁獲サンプルからlength frequencyサンプル200を抽出

5.  さらに各体長ビンから10尾ずつ年齢査定をしてALKの計算

6.  4と5を繰り返してシミュレーション

成長モデル（von Bertalanffy）
=============================

トラフグ瀬戸内系群を参考にしてみた。

雌雄差が少なさそうなので、資源評価書の値を雄雌それぞれほどほどに見て扱った

    L_inf_true <- 560
    K_true <- 0.6
    a0_true <- -0.144
    VB <- function(x, L_inf=L_inf_true, K=K_true, a0=a0_true, 
                   a_sd=0.05, b_sd=0.1, deterministic = FALSE){
      if(isTRUE(deterministic)){
        L_inf*(1-exp(-K*(x-a0)))
      } else {
        sd <- a_sd + b_sd*x
        L_inf*(1-exp(-K*(x-a0)))*exp(rnorm(1, -.5*sd^2, sd))
      }
    }

    curve(VB(x, deterministic = TRUE), xlim = c(0,5), ylim = c(0,1000),
          ylab = "Length (mm)", xlab = "Age", lwd = 3)
    for(i in 1:100){
      points(1,VB(1))
      points(2,VB(2))
      points(3,VB(3))
      points(4,VB(4))
      points(5,VB(5))
    }

![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-1-1.png)

<span
style="color: red; ">分散に線形構造加えると、対数正規の補正ずれる</span>
<span style="color: red; ">**あとでチェック**</span>

#### 最大ビン幅と年齢+グループの設定

    max_length <-700
    max_age <-4

シミュレーション1回の流れ
=========================

資源集団の生成と体長データの生成
--------------------------------

-   なんだが漁獲組成から反映させてしまった
-   適当に指数減少的にしてもいいのかも

<!-- -->

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

漁獲物の抽出
------------

アイデアわかなかったので、300mm-600mmがランダムに1万匹漁獲される

選択曲線がいいのだけど、どう匹数データに換算するか知恵が足りなかった

だって、選択曲線ってあくまで確率だし、選択曲線全てかけて最終的に10000とれるように設定するのも筋が違う気がするし

    sel_id <- which(Sample_true$length>=300 & Sample_true$length<=600)
    random_label <- sample(sel_id, 10000) %>% sort()
    Catch_true <- Sample_true[random_label,]

    hist(Catch_true$age)

![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-4-1.png)

    (prob_Catch_true <- Catch_true$age %>% floor() %>% table()%>% prop.table())

    ## .
    ##      0      1      2      3      4      5      6      7      8 
    ## 0.0125 0.2909 0.2828 0.1993 0.1318 0.0641 0.0169 0.0016 0.0001

Length frequencyサンプルの抽出
------------------------------

    n_sample1 <- 200
    sample1_true <- sample(1:length(Catch_true$id), n_sample1)

    LF_mat_true <- data.frame(id = Catch_true$id[sample1_true],
                              age = Catch_true$age[sample1_true],
                              length = Catch_true$length[sample1_true]
                              )
    LF_mat_true %<>% mutate(bin_L = lencat(length, w=20))
    LF_mat_true$bin_L[(LF_mat_true$bin_L>=max_length)] <- max_length
    head(LF_mat_true)

    ##      id      age   length bin_L
    ## 1 13089 2.338734 427.5456   420
    ## 2  1513 1.617907 373.1175   360
    ## 3  5625 1.598441 410.1616   400
    ## 4 19204 2.417651 488.4532   480
    ## 5 44633 1.630511 358.0441   340
    ## 6 22104 2.227879 394.7670   380

Age-lengthサブサンプルの抽出
----------------------------

各体長幅ビンから10尾ずつ年齢査定してるんだけど、これって細かいほうが自ずと有利だよね

年齢査定の総尾数は固定したほうがいい？？

目的あってのシミュレーションなわけだし、実際の現場に即している方法を真鍋さんに伺う

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

    ## 
    ## 300 320 340 360 380 400 420 440 460 480 500 520 540 560 580 
    ##  10  10  10  10  10  10  10  10  10  10   8   6   7   6   6

ALKの計算
---------

    ALK_freq <- xtabs(~bin_L + age, data = AL_mat)
    ALK_est <- prop.table(ALK_freq, margin = 1)
    alkPlot(ALK_est, type = "area", pal = "gray", showLegend = TRUE,
            leg.cex = .7, xlab = "Total Length(mm)")

![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    alkPlot(ALK_est, type = "bubble", pal = "gray", showLegend = TRUE,
            leg.cex = .7, xlab = "Total Length(mm)")

![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-7-2.png)

    length_n <- xtabs(~bin_L, data = LF_mat_true)
    alkAgeDist(ALK_est, lenA.n = rowSums(ALK_freq), len.n = length_n)

    ##   age   prop          se
    ## 1   0 0.0080 0.008226786
    ## 2   1 0.2555 0.042031115
    ## 3   2 0.2330 0.037873608
    ## 4   3 0.2685 0.040873657
    ## 5   4 0.2350 0.038989646

ざっとこんな感じです

ビン幅を変えての計算とその結果
==============================

    bin_sim <- seq(10,100,10)
    iteration <- 100
    n_sample1 <- 200
    res_ALKest <- res_ALKse <- LF_mat_sim <- AL_mat_sim <- list()
    est_mat <- se_mat <- matrix(NA, ncol = 5, nrow = iteration)

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
        LF_mat_sim[[(ite-1)*10+bb]] <- LF_mat_true
        AL_mat_sim[[(ite-1)*10+bb]] <- AL_mat
      }#for(ite)
      res_ALKest[[bb]] <- est_mat
      res_ALKse[[bb]] <- se_mat
    }#for(bb)

結果
----

パワープレイで恐縮ですが…

    for(i in 1:10){
      boxplot(res_ALKest[[i]])
      points(1:5, prob_Catch_true[1:5], col="red", pch=16, cex=1.6)
    }

![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-1.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-2.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-3.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-4.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-5.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-6.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-7.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-8.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-9.png)![](Simulation_ALK_preliminary_files/figure-markdown_strict/unnamed-chunk-9-10.png)

メモ
====

目的を見失っていないかチェック

なんのためのシミュレーションなのか

現実に即しているようにもっと勉強する

選択率のところ考える

選択率は結局外れ値的にとれる、あるいは歯抜けでとれる体長ビン対策

シミュレーションで自動的に計算するのと、データ見てビン幅とか調整していくのは割と解析上の相性が悪い
