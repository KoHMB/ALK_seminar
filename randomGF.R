VB <- function(x, L_inf=L_inf_true, K=K_true, a0=a0_true, 
               a_sd=0.05, b_sd=0.1, deterministic = FALSE){
  if(isTRUE(deterministic)){
    L_inf*(1-exp(-K*(x-a0)))
  } else {
    sd <- 0.2
    L_inf*(1-exp(-K*(x-a0)))*exp(rnorm(1, -.5*sd^2, sd))
  }
}

VB2 <- function(x, L_inf=L_inf_true, K=K_true, a0=a0_true, 
                a_sd=.3, b_sd=.1, deterministic = FALSE){
  if(isTRUE(deterministic)){
    L_inf*(1-exp(-K*(x-a0)))
  } else {
    sd <- a_sd + b_sd*x
    exp(log(L_inf*(1-exp(-K*(x-a0))))+rnorm(1, -.5*sd^2, sd))
  }
}

curve(VB2(x, K = K_true_vec[7], deterministic = TRUE), xlim = c(0,5), ylim = c(0,1000),
      ylab = "Length (mm)", xlab = "Age", lwd = 3)
for(i in 1:100){
  points(1,VB2(1, K = K_true_vec[7]))
  points(2,VB2(2, K = K_true_vec[7]))
  points(3,VB2(3, K = K_true_vec[7]))
  points(4,VB2(4, K = K_true_vec[7]))
  points(5,VB2(5, K = K_true_vec[7]))
}

a <- numeric()
for(i in 1:10000){
  a[i] <- VB2(5, K = K_true_vec[7])
}
mean(a)
VB2(5, K = K_true_vec[7], deterministic = TRUE)
