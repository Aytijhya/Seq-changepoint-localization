# Reproduction of Figure 2(a) and Table 3 (Setting II)
# Paper: Post-detection inference for sequential changepoint localization
library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Core Experiment Function ---

run_setting2_expt <- function(true_T, theta_min, alpha_univ, alpha_adapt, beta_adapt, iterations, is_viz = FALSE) {
  
  # Parameters for Weighted CUSUM
  A <- 1000
  B <- 100        # Adaptive simulations
  N_sims <- 50   # rt estimation simulations
  theta_grid <- seq(theta_min, theta_min + 9 * 0.2, by = 0.2) 
  weights <- c(exp(-(0:8)/2) - exp(-(1:9)/2), exp(-4.5))
  
  # Helper: Fast Weighted CUSUM
  run_cusum <- function(data, threshold) {
    S <- 0
    for (i in seq_along(data)) {
      lr_inc <- dnorm(data[i], 0.75, 1, log = TRUE) - dnorm(data[i], 0, 1, log = TRUE)
      S <- max(0, S + lr_inc)
      if (exp(S) >= threshold) return(i)
    }
    return(NA)
  }
  
  # Helper: Mt Statistic for Setting II
  calc_Mt <- function(data, t, tau, grid, w, th_min) {
    # Forward: Evidence for P0 against P1 (closest element th_min)
    lr_fwd <- dnorm(data[t:tau], 0, 1, log = TRUE) - dnorm(data[t:tau], th_min, 1, log = TRUE)
    R_t <- exp(max(cumsum(lr_fwd)))
    
    # Backward: Evidence for P1 (mixture) against P0
    if (t > 1) {
      bwd_seg <- data[1:(t-1)]
      S_t <- sum(sapply(seq_along(grid), function(idx) {
        theta <- grid[idx]
        exp(max(cumsum(rev(dnorm(bwd_seg, theta, 1, log = TRUE) - dnorm(bwd_seg, 0, 1, log = TRUE))))) * w[idx]
      }))
    } else { S_t <- 0 }
    return(max(R_t, S_t))
  }
  
  # --- Simulation Loop (Setting II: Adaptive Method Algorithm 2) ---
  
  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
      
      # 1. Data Generation (Setting I & II parameters) [cite: 524, 549]
      obs <- c(rnorm(true_T - 1, 0, 1), rnorm(400, 1, 1))
      tau <- run_cusum(obs, A) 
      if (is.na(tau) || tau < true_T) return(NULL) 
      
      data_tau <- obs[1:tau]
      
      # 2. Estimate rt (Unbiased estimator using Monte Carlo) [cite: 214, 222, 225]
      # Using the (1 + count)/(N + 1) version for strict positivity [cite: 225]
      null_stops <- replicate(N_sims, {
        s <- run_cusum(rnorm(tau + 2, 0, 1), A)
        if(is.na(s)) tau + 2 else s
      })
      
      univ_set <- c(); adapt_set <- c()
      
      # 3. Candidate Inversion Loop [cite: 152, 313, 316]
      for (t in 1:tau) {
        rt_val <- (1 + sum(null_stops >= t)) / (N_sims + 1)
        # Safeguard rt_val from being too small
        rt_val <- max(rt_val, 1/(N_sims + 1)) 
        
        mt_obs <- calc_Mt(data_tau, t, tau, theta_grid, weights, theta_min)
        
        # --- UNIVERSAL METHOD (Equation 10/14) [cite: 215, 271] ---
        if (!is.na(mt_obs) && mt_obs < 2 / (alpha_univ * rt_val)) {
          univ_set <- c(univ_set, t)
        }
        
        # --- ADAPTIVE METHOD (Algorithm 2) [cite: 317, 390] ---
        # a. Build Confidence Sequence (CS) for theta_1 [cite: 355, 359]
        n_t <- tau - t + 1
        seg_mean <- mean(data_tau[t:tau])
        
        # Target coverage 1 - beta*rt [cite: 359]
        cs_err <- max(beta_adapt * rt_val, 1e-6) 
        
        # Howard et al. boundary: nonasymptotic and time-uniform [cite: 355, 356, 785]
        # Radius with numerical safety
        radius <- sqrt( (2 * (n_t + 1) * log(sqrt(n_t + 1) / cs_err)) / max(n_t^2, 1) )
        theta_lower <- max(theta_min, seg_mean - radius)
        
        # b. Simulate B streams under H_{0,t} [cite: 311, 415]
        m_t_sims <- replicate(B, {
          # X_1...X_{t-1} ~ F_theta0, X_t... ~ F_theta_lower [cite: 325, 414]
          sim_obs <- c(rnorm(t - 1, 0, 1), rnorm(250, theta_lower, 1))
          sim_tau <- run_cusum(sim_obs, A)
          
          if (is.na(sim_tau) || sim_tau < t) return(0)
          
          # Compute Mt for simulated sequence [cite: 312, 339, 416]
          calc_Mt(sim_obs[1:sim_tau], t, sim_tau, theta_grid, weights, theta_min)
        })
        
        # c. Inversion check using (1 - alpha*rt) quantile [cite: 315, 341, 433]
        q_val <- quantile(m_t_sims, 1 - (alpha_adapt * rt_val), na.rm = TRUE)
        if (!is.na(mt_obs) && mt_obs <= q_val) {
          adapt_set <- c(adapt_set, t)
        }
      }
      
      # 4. Return results [cite: 530]
      if (is_viz) {
        return(data.frame(Time = 1:tau, Data = data_tau, In_CI = (1:tau) %in% univ_set, Run = i))
      } else {
        return(data.frame(T = true_T, Theta_Min = theta_min, Delay = tau - true_T, 
                          Univ_Sz = length(univ_set), Adapt_Sz = length(adapt_set),
                          Univ_Cov = (true_T %in% univ_set), Adapt_Cov = (true_T %in% adapt_set)))
      }
    }
  }, finally = { stopCluster(cl) })
  
  return(results)
}

# --- 2. Execute Table 3 (Sample of Rows) ---

cat("Running Table 3 simulations (500 iterations, this may take a few minutes/hours)...\n")

table_rows <- list(
  run_setting2_expt(100, 0.75, 0.075, 0.1, 0.025, 100),
  run_setting2_expt(100, 0.9,  0.075, 0.1, 0.025, 100),
  run_setting2_expt(500, 0.75, 0.075, 0.1, 0.025, 100),
  run_setting2_expt(500, 0.9,  0.075, 0.1, 0.025, 100)
)

table_3_summary <- bind_rows(table_rows) %>%
  group_by(T, Theta_Min) %>%
  summarise(Univ_Coverage = mean(Univ_Cov), Adapt_Coverage = mean(Adapt_Cov),
            Univ_Size = mean(Univ_Sz), Adapt_Size = mean(Adapt_Sz), Avg_Delay = mean(Delay))

print("--- REPRODUCED TABLE 3 ---")
print(table_3_summary)

# --- 3. Execute Figure 2(a) Visualization ---

cat("\nGenerating Figure 2(a) runs...\n")
plot_data <- run_setting2_expt(100, 0.9, 0.1, 0.1, 0, 5, is_viz = TRUE)

p=ggplot(plot_data, aes(x = Time, y = Data)) +
  geom_point(size = 0.6) +
  geom_vline(xintercept = 100, color = "black", size = 1) + # True T
  geom_point(data = subset(plot_data, In_CI), aes(x = Time, y = 0), color = "red", shape = 16, size = 1.5) +
  facet_wrap(~Run, ncol = 1, scales = "fixed") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), strip.text = element_blank()) +
  labs(title = "",
       x = "Time", y = "Data")
ggsave("fig2a.png", plot = p, width = 8, height = 5, dpi = 300)
