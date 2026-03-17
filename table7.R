library(dplyr)
library(foreach)
library(doParallel)

run_cusum <- function(data, threshold) {
  Tn <- 0
  for (i in seq_along(data)) {
    Tn <- max(0, Tn + data[i])
    if (Tn >= threshold) return(i)
  }
  return(NA)
}

# Test Statistic Mt (Equation 9) 
calc_Mt <- function(data, t, tau, mu_val) {
  if (t > tau) return(-Inf)
  fwd_data <- data[t:tau]
  fwd_lrs <- cumsum(dnorm(fwd_data, -mu_val, 1, log = TRUE) - dnorm(fwd_data, mu_val, 1, log = TRUE))
  R_max <- exp(max(fwd_lrs))
  if (t > 1) {
    bwd_data <- data[1:(t-1)]
    bwd_lrs <- cumsum(rev(dnorm(bwd_data, mu_val, 1, log = TRUE) - dnorm(bwd_data, -mu_val, 1, log = TRUE)))
    S_max <- exp(max(bwd_lrs))
  } else { S_max <- 0 }
  return(max(R_max, S_max))
}

# --- 3. Simulation Wrapper ---
# This function handles a single run and returns stats + plot data
run_adaptive_sim <- function(obs, tau, d, mu_val, alpha, N_sims) {
  
  data_at_tau <- obs[1:tau]
  null_stops <- replicate(N_sims, {
    s <- run_cusum(rnorm(tau + 2, -mu_val, 1), d)
    ifelse(is.na(s), tau + 2, s)
  })
  
  univ_set <- c()
  
  for (t in 1:tau) {
    rt <- (1 + sum(null_stops >= t)) / (N_sims + 1)
    mt <- calc_Mt(data_at_tau, t, tau, mu_val)
    if (mt < 2 / (alpha * rt)) univ_set <- c(univ_set, t)
  }
  return(univ_set)
}

run_table7_expt <- function(T_true, mu_val, alpha, iter = 100) {
  d <- if(mu_val == 0.25) 8.59 else 7.56 # Thresholds from Page 35
  mu_pre <- -mu_val
  
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(i = 1:iter, .combine = rbind, .packages = c("stats", "dplyr"),
            .export = c("get_wu_ci", "run_adaptive_ci")) %dopar% {
              
              # Generate Setting I Data 
              obs <- c(rnorm(T_true-1, mu_pre, 1), rnorm(500, mu_val, 1))
              
              # CUSUM Detector
              Tn <- 0; tau <- 0
              for (k in 1:length(obs)) {
                Tn <- max(0, Tn + obs[k])
                if (Tn >= d) { tau <- k; break }
              }
              
              # CONDITIONAL REQUIREMENT: Discard iterations where tau < T or tau is NA
              if (tau < T_true || tau == 0) return(NULL)
              
              # Apply both methods to same data
              wu_ci <- get_wu_ci(obs, tau, mu_val, alpha)
              our_ci <- run_adaptive_ci(obs, tau, d, mu_val, alpha, 100)
              
              data.frame(
                T_val = T_true,
                Our_Cov = T_true %in% our_ci,
                Wu_Cov = T_true %in% wu_ci,
                Our_Size = length(our_ci),
                Wu_Size = length(wu_ci)
              )
            }
  }, finally = { stopCluster(cl) })
  return(results)
}

# --- 4. Final Summary ---
res_100 <- run_table7_expt(100, 0.25, 0.1, iter = 500)
print("--- TABLE 7 REPRODUCTION (T=100) ---")
print(res_100 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

res_500 <- run_table7_expt(500, 0.25, 0.1, iter = 500)
print("--- TABLE 7 REPRODUCTION (T=500) ---")
print(res_500 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

res_100 <- run_table7_expt(100, 0.3, 0.1, iter = 500)
print("--- TABLE 7 REPRODUCTION (T=100) ---")
print(res_100 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

res_500 <- run_table7_expt(500, 0.3, 0.1, iter = 500)
print("--- TABLE 7 REPRODUCTION (T=100) ---")
print(res_500 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

