#Experiment in Supplementary Section H.1

library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Parameter Setup ---

alpha <- 0.1     
B <- 100        # Simulations for Adaptive Quantile
N_sims <- 50   # Simulations for r_t estimation
iterations <- 500

# --- 2. Core Logic Functions ---

# CUSUM Detector 
run_cusum <- function(data, threshold) {
  S <- 0
  for (i in seq_along(data)) {
    lr_inc <- dnorm(data[i], 1, 1, log = TRUE) - dnorm(data[i], 0, 1, log = TRUE)
    S <- max(0, S + lr_inc)
    if (exp(S) >= threshold) return(i)
  }
  return(NA)
}

# Universal Test Statistic Mt (Equation 9) 
calc_Mt <- function(data, t, tau) {
  if (t > tau) return(-Inf)
  fwd_data <- data[t:tau]
  fwd_lrs <- cumsum(dnorm(fwd_data, 0, 1, log = TRUE) - dnorm(fwd_data, 1, 1, log = TRUE))
  R_max <- exp(max(fwd_lrs))
  if (t > 1) {
    bwd_data <- data[1:(t-1)]
    bwd_lrs <- cumsum(rev(dnorm(bwd_data, 1, 1, log = TRUE) - dnorm(bwd_data, 0, 1, log = TRUE)))
    S_max <- exp(max(bwd_lrs))
  } else { S_max <- 0 }
  return(max(R_max, S_max))
}

# --- 3. Simulation Wrapper ---
# This function handles a single run and returns stats + plot data
run_single_sim <- function(true_T, A) {
  obs <- c(rnorm(true_T - 1, 0, 1), rnorm(500, 1, 1))
  tau <- run_cusum(obs, A)
  if (is.na(tau) || tau < true_T) return(NULL) 
  
  data_at_tau <- obs[1:tau]
  null_stops <- replicate(N_sims, {
    s <- run_cusum(rnorm(tau + 2, 0, 1), A)
    ifelse(is.na(s), tau + 2, s)
  })
  
  univ_set <- c()
  
  for (t in 1:tau) {
    rt <- (1 + sum(null_stops >= t)) / (N_sims + 1)
    mt <- calc_Mt(data_at_tau, t, tau)
    if (mt < 2 / (alpha * rt)) 
      univ_set <- c(univ_set, t)
  }
  
  res <- list(stats = data.frame(T = true_T, Delay = tau - true_T, 
                                 Univ_Size = length(univ_set), 
                                 Univ_Cover = (true_T %in% univ_set)))
  return(res)
}

# --- 4. Parallel Execution for Table ---
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores); registerDoParallel(cl)

table_results <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  r1 <- run_single_sim(50,100)$stats
  r2 <- run_single_sim(1000,100)$stats
  r3 <- run_single_sim(2000,100)$stats
  r4 <- run_single_sim(500,1000)$stats
  r5 <- run_single_sim(10000,1000)$stats
  r6 <- run_single_sim(20000,1000)$stats
  rbind(r1, r2, r3, r4, r5, r6)
}
stopCluster(cl)

# --- 6. Printing Results ---

# Table Summary 
print("--- REPRODUCED TABLE 2: SETTING I ---")
summary_tab <- table_results %>% group_by(T) %>%
  summarise(Coverage_Univ = mean(Univ_Cover), 
            Size_Univ = mean(Univ_Size), Delay = mean(Delay))
print(summary_tab)







