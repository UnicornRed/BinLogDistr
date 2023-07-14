p_0 <- 0.2

log_0_prepering <- function(p = 0.5){
  tdistr <- (1 - p_0) * (-p) / log(1 - p)
  sum_distr <- p_0
  sum <- c(p_0)
  k <- 2
  
  while(tdistr > 1e-10) {
    sum_distr <- sum_distr + tdistr
    sum <- c(sum, sum_distr)
    tdistr <- tdistr * p / k * (k - 1)
    k <- k + 1
  }
  
  return(sum)
}

rlog_0 <- function(p = 0.5, sum) {
  x <- runif(1)
  j <- 1
  
  while (x > sum[j]) {
    j <- j + 1
  }
  
  return(j - 1)
}

rnlog_0 <- function(n = 1, p = 0.5, sum){
  v <- replicate(n, rlog_0(p, sum))
  
  return(v)
}

table_log <- log_0_prepering(0.75)
v <- rnlog_0(10000, 0.75, table_log)
hist(v[v < 9], probability = TRUE, breaks = seq(-0.5, max(v) + 0.5, 1), xlim = c(-0.5, 8.5), xlab = "Sample", main = "")


log_prepering <- function(p = 0.5){
  tdistr <- (-p) / log(1 - p)
  sum_distr <- tdistr
  sum <- c(sum_distr)
  k <- 2
  
  while(tdistr > 1e-10) {
    tdistr <- tdistr * p / k * (k - 1)
    sum_distr <- sum_distr + tdistr
    sum <- c(sum, sum_distr)
    k <- k + 1
  }
  
  return(sum)
}

rlog <- function(sum) {
  x <- runif(1)
  j <- 1
  
  while (x > sum[j]) {
    j <- j + 1
  }
  
  return(j)
}

rnlog <- function(n = 1, sum){
  v <- replicate(n, rlog(sum))
  
  if(n == 0)
    v <- 0;
  
  return(v)
}

rlogbinom <- function(num = 1, sumlog, n = 1, p = 0.5) {
  res <- c()
  
  for(i in rbinom(num, n, p)){
    res <- c(res, sum(0, rnlog(i, sumlog)))
  }
  
  return(res)
}

rbinomlog <- function(num = 1, sumlog, n = 1, p = 0.5) {
  res <- c()
  
  for(i in rnlog(num, sumlog)){
    res <- c(res, sum(0, rbinom(i, n, p)))
  }
  
  return(res)
}

get_prob_0 <- function(p, n, q){
  return ((1 - p)^n)
}

get_prob_1 <- function(p, n, q){
  return (n * p * (1 - p)^(n - 1) * (-q) / log(1 - q))
}

get_prob_2 <- function(p, n, q){
  return (n * p * (1 - p)^(n - 2) * (q^2) / log(1 - q) * ((n - 1) * p / log(1 - q) - (1 - p)) / 2)
}

get_prob_3 <- function(p, n, q){
  return (n * p * (1 - p)^(n - 3) * (-q^3) / log(1 - q) * ((n - 1) * (n - 2) * p^2  / (log(1 - q))^2 - 3 * (n - 1) * p * (1 - p) / log(1 - q) + 2 * (1 - p)^2) / 6)
}

get_prob_4 <- function(p, n, q){
  return (n * p * (1 - p)^(n - 4) * (q^4) / log(1 - q) * ((n - 1) * (n - 2) * (n - 3) * p^3  / (log(1 - q))^3 - 6 * (n - 1) * (n - 2) * p^2 * (1 - p) / (log(1 - q))^2 + 11 * (n - 1) * p * (1 - p)^2 / log(1 - q) - 6 * (1 - p)^3) / 24)
}

get_prob <- function(p, n, q, k){
  if (k == 0){
    return(get_prob_0(p, n, q))
  }
  
  if (k == 1){
    return(get_prob_1(p, n, q))
  }
  
  if (k == 2){
    return(get_prob_2(p, n, q))
  }
  
  if (k == 3){
    return(get_prob_3(p, n, q))
  }
  
  if (k == 4){
    return(get_prob_4(p, n, q))
  }
  
  if (k >= 5){
    return(1 - get_prob_0(p, n, q) - get_prob_1(p, n, q) - get_prob_2(p, n, q) - get_prob_3(p, n, q) - get_prob_4(p, n, q))
  }
}

m_log_lik <- function(params, x.in, E, n) {
  q = params[1]
  res <- 0
  p <- - log(1 - q) * (1 - q) / n / q * E
  
  for (i in x.in){
    prob <- get_prob(p, n, q, i)
    
    if (prob < 0)
      prob <- get_prob(0.99, n, q, i)
    
    res <- res + log(prob)
  }
  
  return(-res)
}

generate_sample <- function(num_k){
  sam <- c()
  
  for (i in 1:length(num_k)){
    sam <- c(sam, rep.int(i - 1, num_k[i]))
  }
  
  return(sam)
}

my_chisq <- function(exp_prob, prob){
  res <- 0
  
  for (i in 1:length(prob)){
    res <- res + (exp_prob[i] - prob[i])^2/prob[i]
  }
  
  return(res)
}

hist_make <- function (n, exp_prob, get_hist){
  # q <- 0.2
  # p <- 0.25
  # sum <- log_prepering(q)
  # 
  # v <-rlogbinom(10000, sum, n, p)
  
  res <- c(n)
  sam <- generate_sample(exp_prob)
  exp_prob <- exp_prob / 100
  
  opt.res <- optimise(function(q) m_log_lik(q, sam, mean(sam), n),lower = 0, upper = 1)
  q <-opt.res$minimum
  res <- c(res, q)
  
  p <- - log(1 - q) * (1 - q) / n / q * mean(sam)
  res <- c(res, p)
  
  prob <- c()
  
  for (i in 0:4){
    prob <- c(prob, get_prob(p, n, q, i))
  }
  
  prob <- c(prob, 1 - sum(prob))
  
  max_el_x = max(sam)
  max_el_y = max(exp_prob, prob)
  
  if (get_hist == TRUE){
    hist.default(sam, probability = TRUE, breaks = seq(-0.5, max_el_x + 0.5, 1), xlim = c(-0.5,6.5), ylim = c(0, max_el_y), xlab = "Sample", main = "")
    points(0:5, prob)
  }
  
  for (i in 1:length(prob)){
    while (100 * prob[i] < 5 && i < length(prob)){
      prob[i] <- prob[i] + prob[i + 1]
      prob <- prob[-(i + 1)]
      exp_prob[i] <- exp_prob[i] + exp_prob[i + 1]
      exp_prob <- exp_prob[-(i + 1)]
    }
  }
  
  if (100 * prob[length(prob)] < 5){
    prob[length(prob) - 1] <- prob[length(prob) - 1] + prob[length(prob)]
    prob <- prob[-length(prob)]
    exp_prob[length(exp_prob) - 1] <- exp_prob[length(exp_prob) - 1] + exp_prob[length(exp_prob)]
    exp_prob <- exp_prob[-length(exp_prob)]
  }
  
  df <- length(prob) - 3

  if (df < 1){
    df <- 1
  }
  
  chi <- 100 * my_chisq(exp_prob, prob)
  res <- c(res, 1 - pchisq(chi, df))
  
  return(res)
}

vec <- c(66, 31, 2, 1, 0, 0)
vec <- c(vec, 50, 35, 13, 2, 0, 0)
vec <- c(vec, 41, 39, 17, 2, 1, 0)
vec <- c(vec, 27, 39, 29, 3, 1, 1)
vec <- c(vec, 22, 22, 32, 15, 6, 3)
vec <- c(vec, 33, 39, 18, 8, 1, 1)
vec <- c(vec, 21, 29, 21, 14, 10, 5)
vec <- c(vec, 17, 24, 32, 11, 10, 6)
vec <- c(vec, 15, 17, 29, 18, 11, 10)
vec <- c(vec, 13, 14, 24, 24, 12, 13)
vec <- c(vec, 68, 25, 7, 0, 0, 0)
vec <- c(vec, 74, 19, 5, 2, 0, 0)
vec <- c(vec, 59, 24, 16, 1, 0, 0)
vec <- c(vec, 48, 33, 11, 3, 2, 3)
vec <- c(vec, 59, 31, 5, 4, 1, 0)
vec <- c(vec, 37, 31, 22, 6, 2, 2)
vec <- c(vec, 35, 37, 17, 5, 3, 3)
vec <- c(vec, 26, 36, 19, 10, 5, 4)
vec <- c(vec, 19, 33, 25, 11, 7, 5)

m <- matrix(vec, 19, 6, byrow = TRUE)

df <- data.frame(n = c(0), q = c(0), p = c(0), p_value = c(0))

for (i in 1:19){
  max <- 0
  n <- 1
  for (j in 1:50){
    res <- hist_make(j, m[i,], FALSE)
    
    if (res[4] > max){
      max <- res[4]
      n <- j
    }
  }
  res <- hist_make(n, m[i,], FALSE)
  df[nrow(df) + 1, ] <- res
  print(res)
}

df <- df[-1, ]

pq_n <-function(n){
  p_n <- c()
  q_n <- c()
  
  for (j in 1:10){
    res <- hist_make(j, m[n,], FALSE)
    
    p_n <- c(p_n, res[3])
    q_n <- c(q_n, res[2])
  }
  
  plot(1:10, p_n, type = "l", col = "red", ylim = c(0, max(p_n, q_n)))
  lines(1:10, q_n, type = "l", col = "blue")
}

pq_n(8)

plot(seq(0, 45, 5), df$q[1:10], type = "l", col = "red", ylim = c(0, max(df$q, df$p)))
lines(seq(0, 40, 5), df$q[11:19], type = "l", col = "blue")

print(df)

print(hist_make(2.46, m[3,], TRUE))

p_il_invitro <- c(0.32, 0.20, 0.21, 0.13, 0.05, 0.04, 0.02, 0.01, 0.01)
q_il_invitro <- c(0.32, 0.39, 0.25, 0.41, 0.22, 0.19, 0.22, 0.19, 0.09)

plot(seq(0, 40, 5), p_il_invitro, type = "l", col = "red", ylim = c(0, max(p_il_invitro, q_il_invitro)), main = "Red - P; Blue - q", xlab = "Dose, Gy", ylab = "Probability")
lines(seq(0, 40, 5), q_il_invitro, type = "l", col = "blue")

p_il_invivo  <- c(0.27, 0.32, 0.28, 0.21, 0.14, 0.04, 0.02, 0.01, 0.01, 0.00)
q_il_invivo  <- c(0.47, 0.55, 0.28, 0.12, 0.17, 0.13, 0.34, 0.19, 0.19, 0.12)

plot(seq(0, 45, 5), p_il_invivo, type = "l", col = "red", ylim = c(0, max(p_il_invivo, q_il_invivo)), main = "Red - P; Blue - q", xlab = "Dose, Gy", ylab = "Probability")
lines(seq(0, 45, 5), q_il_invivo, type = "l", col = "blue")

n_vect <- c(1, 1.36, 2.46, 5.02, 10.92, 24.74, 57.63, 137.08, 331.22, 810.31, 1, 1.36, 2.46, 5.02, 10.92, 24.74, 57.63, 137.08, 331.22)

df <- data.frame(n = c(0), q = c(0), p = c(0), p_value = c(0))

for (i in 1:19){
  res <- hist_make(n_vect[i], m[i,], FALSE)
  df[nrow(df) + 1, ] <- res
}

df <- df[-1, ]

plot(seq(0, 45, 5), df$p[1:10] * df$n[1:10], type = "l", col = "blue", ylab = "Average", xlab = "Dose, Gy", main = "Red - in vitro; Blue - in vivo")
plot(seq(0, 40, 5), df$p[11:19] * df$n[11:19], type = "l", col = "black", ylab = "Average", xlab = "Dose, Gy", main = "in vitro")
lines(seq(0, 40, 5), df$p[11:19] * df$n[11:19], type = "l", col = "black")

n_vect <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 1, 5, 10, 15, 20, 25, 30, 35, 40)

df <- data.frame(n = c(0), q = c(0), p = c(0), p_value = c(0))

for (i in 1:19){
  res <- hist_make(n_vect[i], m[i,], FALSE)
  df[nrow(df) + 1, ] <- res
}

df <- df[-1, ]

lines(seq(0, 45, 5), df$p[1:10] * df$n[1:10], type = "b", col = "blue")
lines(seq(0, 40, 5), df$p[11:19] * df$n[11:19], type = "b", col = "black")

mid <- c(0.38, 0.67, 0.83, 1.15, 1.70, 1.08, 1.83, 1.91, 2.23, 2.47, 0.39, 0.35, 0.59, 0.87, 0.57, 1.11, 1.13, 1.44, 1.69)

lines(seq(0, 45, 5), mid[1:10], type = "l", lty = 2, col = "green")
lines(seq(0, 40, 5), mid[11:19], type = "l", lty = 2, col = "black")
