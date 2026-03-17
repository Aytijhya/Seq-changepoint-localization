library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Wu (2007) Asymptotic CI (Source: Page 35, Section 6) ---
get_wu_ci <- function(obs, tau_stop, theta_1, alpha) {
  Tn <- numeric(tau_stop + 1)
  for (i in 1:tau_stop) Tn[i+1] <- max(0, Tn[i] + obs[i])
  
  # Asymptotic constants from Wu (2007)
  v_hat <- max(which(Tn[1:tau_stop] == 0)) 
  s_approx <- log(alpha) * (1 / (sqrt(2) * theta_1) + 0.088)
  c_approx <- -log(1 - sqrt(1 - alpha)) / (2 * theta_1) - 0.583
  
  zeros <- which(Tn[1:v_hat] == 0)
  s_int <- ceiling(abs(s_approx))
  Ls <- if(length(zeros) >= s_int) rev(zeros)[s_int] else 1
  
  Vc_indices <- which(Tn[(v_hat + 1):(tau_stop + 1)] <= c_approx) + v_hat
  return(unique(c(Vc_indices, Ls:(v_hat - 1))))
}

# --- 2. Corrected Adaptive CI (Source: Algorithm 1 & 2) ---
run_adaptive_ci <- function(obs, tau_stop, A, theta_1, alpha, N_sims) {
  # Estimate rt with strict positivity safeguard 
  null_stops <- replicate(N_sims, {
    Tn <- 0; i <- 0
    while(Tn < A && i < (tau_stop + 1)) {
      i <- i + 1
      Tn <- max(0, Tn + rnorm(1, -theta_1, 1)) 
    }
    i
  })
  
  adapt_set <- c()
  for (t in 1:tau_stop) {
    # Strictly positive rt estimator 
    rt <- (1 + sum(null_stops >= t)) / (N_sims + 1)
    
    # Mt calculation (Equation 9) [cite: 201]
    fwd_data <- obs[t:tau_stop]
    # Log-space sum to prevent overflow, then exp [cite: 187]
    log_fwd <- cumsum(-2 * theta_1 * fwd_data)
    max_fwd <- exp((max(log_fwd)))
    
    max_bwd <- 0
    if (t > 1) {
      bwd_data <- obs[1:(t-1)]
      log_bwd <- cumsum(rev(2 * theta_1 * bwd_data))
      max_bwd <- exp(max(log_bwd))
    }
    mt_obs <- max(max_fwd, max_bwd)
    
    # Simulation-based quantile [cite: 315]
    m_t_sims <- replicate(50, {
      sim_obs <- c(rnorm(t - 1, -theta_1, 1), rnorm(300, theta_1, 1))
      Tn_sim <- 0; tau_sim <- 0
      for (k in 1:length(sim_obs)) {
        Tn_sim <- max(0, Tn_sim + sim_obs[k])
        if (Tn_sim >= A) { tau_sim <- k; break }
      }
      if (tau_sim < t) return(0)
      
      sim_seg <- sim_obs[1:tau_sim]
      s_fwd <- exp(max(cumsum(-2 * theta_1 * sim_seg[t:tau_sim])))
      s_bwd <- if(t > 1) exp(max(cumsum(rev(2 * theta_1 * sim_seg[1:(t-1)])))) else 0
      max(s_fwd, s_bwd)
    })
    
    # Avoid NA/NaN in quantile 
    if (mt_obs <= quantile(c(mt_obs,m_t_sims),  1 - alpha* rt, na.rm = TRUE)) {
      adapt_set <- c(adapt_set, t)
    }
  }
  return(adapt_set)
}

# --- 3. Parallel Comparison Runner ---
run_table8_expt <- function(T_true, mu_val, iter = 100) {
  d <- if(mu_val == 0.25) 8.59 else 7.56 # Thresholds from Page 35
  mu_pre <- -mu_val
  
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(i = 1:iter, .combine = rbind, .packages = c("stats", "dplyr"),
            .export = c("get_wu_ci", "run_adaptive_ci")) %dopar% {
              
              # Generate Setting I Data [cite: 524]
              obs <- c(rnorm(T_true, mu_pre, 1), rnorm(500, mu_val, 1))
              
              # CUSUM Detector [cite: 670]
              Tn <- 0; tau <- 0
              for (k in 1:length(obs)) {
                Tn <- max(0, Tn + obs[k])
                if (Tn >= d) { tau <- k; break }
              }
              
              # Requirement: Conditional on tau >= T 
              if (tau < T_true) return(NULL) 
              
              # Apply both methods to same data
              wu_ci <- get_wu_ci(obs, tau, mu_val, 0.05)
              our_ci <- run_adaptive_ci(obs, tau, d, mu_val, 0.05, 100)
              
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
res_100 <- run_table8_expt(100, 0.25, iter = 500)
print("--- TABLE 8 REPRODUCTION (T=100) ---")
print(res_100 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

res_500 <- run_table8_expt(500, 0.25, iter = 200)
print("--- TABLE 8 REPRODUCTION (T=100) ---")
print(res_500 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

res_100 <- run_table8_expt(100, 0.3, iter = 200)
print("--- TABLE 8 REPRODUCTION (T=100) ---")
print(res_100 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

res_500 <- run_table8_expt(500, 0.3, iter = 200)
print("--- TABLE 8 REPRODUCTION (T=100) ---")
print(res_500 %>% summarise(
  Our_Cond_Coverage = mean(Our_Cov), 
  Wu_Cond_Coverage = mean(Wu_Cov),  
  Our_Size = mean(Our_Size),        
  Wu_Size = mean(Wu_Size)            
))

