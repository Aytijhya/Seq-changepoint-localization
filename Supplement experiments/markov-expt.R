#Experiment in Supplementary Section H.2
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Global Setup ---
alpha_val <- 0.1
A_thresh <- 1000
N_sims <- 100 # For estimating r_t
B_sims <- 100  # For estimating adaptive M_t quantile

# Transition probabilities from State 1 (as per the formulas)
p0 <- 0.2  # P(X=0 | X_prev=1) under P0
p1 <- 0.6  # P(X=0 | X_prev=1) under P1

# Full matrices (Transitions from 0 are identical so they provide 0 log-evidence)
P0 <- matrix(c(0.5, 0.5, p0, 1 - p0), nrow = 2, byrow = TRUE)
P1 <- matrix(c(0.5, 0.5, p1, 1 - p1), nrow = 2, byrow = TRUE)

# --- 2. High-Performance Core Functions ---

# Fast Markov Generator
sim_markov_fast <- function(n, P, start_state = 1) {
  X <- integer(n)
  X[1] <- start_state
  if (n > 1) {
    U <- runif(n - 1)
    for (i in 2:n) {
      X[i] <- if (U[i-1] < P[X[i-1] + 1, 1]) 0 else 1
    }
  }
  return(X)
}

# Fast Detector using exact log-LR logic for n10 and n11
run_markov_detector_user <- function(data, A, p0_val, p1_val) {
  n <- length(data)
  if (n < 2) return(NA)
  
  log_A <- log(A)
  log_e <- 0
  
  log_lr_10 <- log(p1_val / p0_val)
  log_lr_11 <- log((1 - p1_val) / (1 - p0_val))
  
  for (i in 2:n) {
    # Only accumulate evidence when transitioning FROM state 1
    if (data[i-1] == 1 && data[i] == 0) {
      lr_val <- log_lr_10
    } else if (data[i-1] == 1 && data[i] == 1) {
      lr_val <- log_lr_11
    } else {
      lr_val <- 0 # P_1(x|0) = P_0(x|0), cancels out
    }
    
    log_e <- max(0, log_e) + lr_val
    if (log_e >= log_A) return(i)
  }
  return(NA)
}


calc_Mt_user <- function(data, t, tau, p0_val, p1_val) {
  if (t > tau) return(-Inf)
  
  # Forward R^{(t)}_n: Collection t <= n <= tau
  max_R <- 0.5 # Default empty product is 1/2
  if (t <= tau && tau > 1) {
    idx_start <- max(1, t - 1)
    prev_fwd <- data[idx_start:(tau - 1)]
    curr_fwd <- data[(idx_start + 1):tau]
    
    n10_fwd <- cumsum((prev_fwd == 1) & (curr_fwd == 0))
    n11_fwd <- cumsum((prev_fwd == 1) & (curr_fwd == 1))
    
    R_vals <- 0.5 * (p0_val / p1_val)^n10_fwd * ((1 - p0_val) / (1 - p1_val))^n11_fwd
    max_R <- max(max(R_vals), 0.5)
  }
  
  # Backward S^{(t)}_n: Collection 1 <= n <= t-1
  max_S <- 0.5 # Default empty product is 1/2
  if (t > 2) {
    # To compute for segments n to t-1, we read backwards from t-1 down to 2
    prev_bwd <- data[(t - 2):1]
    curr_bwd <- data[(t - 1):2]
    
    n10_bwd <- cumsum((prev_bwd == 1) & (curr_bwd == 0))
    n11_bwd <- cumsum((prev_bwd == 1) & (curr_bwd == 1))
    
    # Notice the flipped ratios for the backward hypothesis testing
    S_vals <- 0.5 * (p1_val / p0_val)^n10_bwd * ((1 - p1_val) / (1 - p0_val))^n11_bwd
    max_S <- max(max(S_vals), 0.5)
  }
  
  return(max(max_R, max_S))
}

# --- 3. Parallel Adaptive Runner ---

run_adaptive_markov <- function(true_T, iter = 100) {
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  export_vars <- c("sim_markov_fast", "run_markov_detector_user", 
                   "calc_Mt_user", "P0", "P1", "p0", "p1", "A_thresh", 
                   "alpha_val", "N_sims", "B_sims")
  
  results <- tryCatch({
    foreach(i = 1:iter, .combine = rbind, .packages = c("stats", "dplyr"), .export = export_vars) %dopar% {
      
      # 1. Simulate Stream
      pre_data <- sim_markov_fast(true_T, P0)
      post_data <- sim_markov_fast(300, P1, start_state = pre_data[true_T])
      obs <- c(pre_data, post_data[-1])
      
      # 2. Detect
      tau <- run_markov_detector_user(obs, A_thresh, p0, p1)
      if (is.na(tau) || tau < true_T) return(NULL) 
      
      # 3. Simulate r_t
      null_stops <- replicate(N_sims, {
        s <- run_markov_detector_user(sim_markov_fast(tau + 1, P0), A_thresh, p0, p1)
        if(is.na(s)) tau + 1 else s
      })
      
      in_ci <- logical(tau)
      obs_tau_data <- obs[1:tau]
      
      # 4. Invert Test Statistic
      for (t_cand in 1:tau) {
        rt <- (1 + sum(null_stops >= t_cand)) / (N_sims + 1)
        mt_obs <- calc_Mt_user(obs_tau_data, t_cand, tau, p0, p1)
        
        # Simulate adaptive threshold under H_{0,t}
        mt_sims <- numeric(B_sims)
        for (b in 1:B_sims) {
          sim_pre <- if(t_cand > 1) sim_markov_fast(t_cand, P0) else sim_markov_fast(1, P0)
          sim_post <- sim_markov_fast(200, P1, start_state = sim_pre[length(sim_pre)])
          sim_obs <- c(sim_pre, sim_post[-1])
          
          sim_tau <- run_markov_detector_user(sim_obs, A_thresh, p0, p1)
          if (!is.na(sim_tau) && sim_tau >= t_cand) {
            mt_sims[b] <- calc_Mt_user(sim_obs[1:sim_tau], t_cand, sim_tau, p0, p1)
          } else {
            mt_sims[b] <- 0
          }
        }
        
        threshold <- quantile(mt_sims, max(0, min(1, 1 - alpha_val * rt)), na.rm = TRUE)
        if (mt_obs <= threshold) in_ci[t_cand] <- TRUE
      }
      
      ci_idx <- which(in_ci)
      data.frame(T_val = true_T, Covered = (true_T %in% ci_idx), Size = length(ci_idx), Delay = tau - true_T)
    }
  }, finally = { stopCluster(cl) })
  
  return(results)
}

# --- 4. Execution ---
cat("Running Fast Adaptive Markovian Simulation...\n")
res_100 <- run_adaptive_markov(100, iter = 100) 

print("--- ADAPTIVE MARKOV METHOD RESULTS (T=100) ---")
print(res_100 %>% summarise(
  T = 100,
  Cond_Coverage = mean(Covered), 
  Avg_Size = mean(Size),
  Avg_Delay = mean(Delay)
))

cat("Running Fast Adaptive Markovian Simulation...\n")
res_500 <- run_adaptive_markov(500, iter = 100) 

print("--- ADAPTIVE MARKOV METHOD RESULTS (T=500) ---")
print(res_500 %>% summarise(
  T = 500,
  Cond_Coverage = mean(Covered), 
  Avg_Size = mean(Size),
  Avg_Delay = mean(Delay)
))
