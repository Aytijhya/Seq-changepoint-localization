library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Parameter Setup ---

alpha_univ <- 0.1     # Target for Setting III Universal 
alpha_adapt <- 0.05   # Adaptive alpha
beta_adapt <- 0.025   # Adaptive beta 
gamma_adapt <- 0.025  # Adaptive gamma 
A <- 1000             # Detection Threshold
B <- 100              # Adaptive simulations 
N_sims <- 50         # rt estimation simulations 
iterations <- 500     # For Table 4 

# Weight setup for Theta1 as in Setting II 
theta_grid_post <- seq(0.9, 0.9 + 9 * 0.2, by = 0.2) 
w_post <- c(exp(-(0:8)/2) - exp(-(1:9)/2), exp(-4.5))

# Weight setup for Theta0 (Pre-change) 
theta_grid_pre <- seq(0.1, 0.1 - 9 * 0.2, by = -0.2)
w_pre <- c(exp(-(0:8)/2) - exp(-(1:9)/2), exp(-4.5))

# --- 2. Core Functions ---

# Weighted CUSUM-type detector (Equation 23) 
run_wcusum_composite <- function(data, threshold, grid_post, weights_post, theta0_star) {
  S <- numeric(length(grid_post))
  for (j in seq_along(data)) {
    # LR of P1(theta) vs P0(theta0_star) 
    lr_inc <- dnorm(data[j], grid_post, 1, log = TRUE) - dnorm(data[j], theta0_star, 1, log = TRUE)
    S <- pmax(0, S + lr_inc)
    if (sum(exp(S) * weights_post) >= threshold) return(j)
  }
  return(NA)
}

# Universal Test Statistic Mt for Setting III
calc_Mt_setting3 <- function(data, t, tau, grid_pre, w_pre, grid_post, w_post, theta0_star, theta1_star) {
  if (t > tau) return(-Inf)
  
  # Forward Segment (t to tau): Mixture over Theta0 vs Theta1_star
  fwd_data <- data[t:tau]
  R_t <- sum(sapply(seq_along(grid_pre), function(i) {
    exp(max(cumsum(dnorm(fwd_data, grid_pre[i], 1, log = TRUE) - dnorm(fwd_data, theta1_star, 1, log = TRUE)))) * w_pre[i]
  }))
  
  # Backward Segment (1 to t-1): Mixture over Theta1 vs Theta0_star
  if (t > 1) {
    bwd_data <- data[1:(t-1)]
    S_t <- sum(sapply(seq_along(grid_post), function(i) {
      exp(max(cumsum(rev(dnorm(bwd_data, grid_post[i], 1, log = TRUE) - dnorm(bwd_data, theta0_star, 1, log = TRUE))))) * w_post[i]
    }))
  } else { S_t <- 0 }
  
  return(max(R_t, S_t))
}

# --- 3. Experiment Runner ---

run_setting3_expt <- function(true_T, theta0_star, theta1_star, alpha_univ, alpha_adapt, beta_adapt, gamma_adapt, is_viz = FALSE) {
  # Generate true data: N(0,1) then N(1,1) 
  obs <- c(rnorm(true_T - 1, 0, 1), rnorm(400, 1, 1))
  tau <- run_wcusum_composite(obs, A, theta_grid_post, w_post, theta0_star)
  if (is.na(tau) || tau < true_T) return(NULL)
  
  data_tau <- obs[1:tau]
  
  # rt* estimation under H0 (using closest element theta0_star)
  null_stops <- replicate(N_sims, {
    s <- run_wcusum_composite(rnorm(tau + 2, theta0_star, 1), A, theta_grid_post, w_post, theta0_star)
    if(is.na(s)) tau + 2 else s
  })
  
  univ_set <- c(); adapt_set <- c()
  
  for (t in 1:tau) {
    rt_star <- (1 + sum(null_stops >= t)) / (N_sims + 1)
    mt_obs <- calc_Mt_setting3(data_tau, t, tau, theta_grid_pre, w_pre, theta_grid_post, w_post, theta0_star, theta1_star)
    
    # 1. Universal Method (Equation 14) 
    if (mt_obs < 2 / (alpha_univ * rt_star)) univ_set <- c(univ_set, t)
    
    # 2. Adaptive Method (Algorithm 3)
    # In practice, discretized CS is used 
    # Optimization: Check boundary elements
    cs_err_post <- beta_adapt * rt_star
    cs_err_pre <- gamma_adapt * rt_star
    
    # Simple CS bounds 
    n_post <- tau - t + 1
    n_pre <- t - 1
    post_mu <- mean(data_tau[t:tau]); pre_mu <- if(t > 1) mean(data_tau[1:(t-1)]) else 0
    rad_post <- sqrt((log(log(2*n_post)) + 0.72 * log(10.4 / cs_err_post)) / max(n_post, 1) )
    rad_pre <- qnorm(1-cs_err_pre/2)/sqrt( max(n_pre, 1))
    
    theta_post_bound <- max(theta1_star, post_mu - rad_post)
    theta_pre_bound <- min(theta0_star, pre_mu + rad_pre)
    
    m_t_sims <- replicate(B, {
      sim_obs <- c(rnorm(t - 1, theta_pre_bound, 1), rnorm(250, theta_post_bound, 1))
      sim_tau <- run_wcusum_composite(sim_obs, A, theta_grid_post, w_post, theta0_star)
      if(is.na(sim_tau) || sim_tau < t) return(0)
      calc_Mt_setting3(sim_obs[1:sim_tau], t, sim_tau, theta_grid_pre, w_pre, theta_grid_post, w_post, theta0_star, theta1_star)
    })
    
    if (mt_obs <= quantile(m_t_sims, 1 - alpha_adapt * rt_star, na.rm = TRUE)) adapt_set <- c(adapt_set, t)
  }
  
  if (is_viz) {
    return(data.frame(Time = 1:tau, Data = data_tau, In_CI = (1:tau) %in% univ_set, Run = 1))
  } else {
    return(data.frame(T = true_T, Delay = tau - true_T, 
                      Univ_Sz = length(univ_set), Adapt_Sz = length(adapt_set),
                      Univ_Cov = (true_T %in% univ_set), Adapt_Cov = (true_T %in% adapt_set),
                      T_hat_error = abs(median(univ_set) - true_T)))
  }
}

# --- 4. Execution ---

# Parallel setup for Table 4 
cl <- makeCluster(parallel::detectCores()-1); registerDoParallel(cl)
table_data_4 <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  # Using first row params: T=100, Theta0=[0, 0.25], Theta1=[0.75, Inf) 
  run_setting3_expt(100, 0.25, 0.75, 0.1,0.05,0.025,0.025)
}
stopCluster(cl)

# Summary Output for Table 4 
print("--- TABLE 4 REPRODUCTION: SETTING III (T=100) ---")
table_4_summary <- table_data_4 %>% 
  summarise(Univ_Coverage = mean(Univ_Cov), Adapt_Coverage = mean(Adapt_Cov),
            Univ_Size = mean(Univ_Sz), Adapt_Size = mean(Adapt_Sz), 
            Avg_Abs_Dev = mean(T_hat_error), Avg_Delay = mean(Delay))
print(table_4_summary)

###

l <- makeCluster(parallel::detectCores()-1); registerDoParallel(cl)
table_data_4 <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  # Using first row params: T=100, Theta0=[0, 0.25], Theta1=[0.75, Inf) 
  run_setting3_expt(100, 0.1, 0.9, 0.1,0.05,0.025,0.025)
}
stopCluster(cl)

# Summary Output for Table 4 
print("--- TABLE 4 REPRODUCTION: SETTING III (T=100) ---")
table_4_summary <- table_data_4 %>% 
  summarise(Univ_Coverage = mean(Univ_Cov), Adapt_Coverage = mean(Adapt_Cov),
            Univ_Size = mean(Univ_Sz), Adapt_Size = mean(Adapt_Sz), 
            Avg_Abs_Dev = mean(T_hat_error), Avg_Delay = mean(Delay))
print(table_4_summary)


# --- Define the Parameter Configurations for Table 4 ---
# Based on the specifications in Table 4 
configs <- list(
  list(T = 100, theta0_star = 0.25, theta1_star = 0.75, name = "Row 1: T=100, [0,0.25], [0.75,inf)"),
  list(T = 100, theta0_star = 0.1,  theta1_star = 0.9,  name = "Row 2: T=100, [0,0.1], [0.9,inf)"),
  list(T = 500, theta0_star = 0.25, theta1_star = 0.75, name = "Row 3: T=500, [0,0.25], [0.75,inf)"),
  list(T = 500, theta0_star = 0.1,  theta1_star = 0.9,  name = "Row 4: T=500, [0,0.1], [0.9,inf)")
)

# --- Parallel Execution across all Configurations ---
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

# We use a nested loop structure or a combined list to run all experiments
# Warning: Running 500 iterations for all 4 rows is computationally heavy 
final_table_4 <- foreach(conf = configs, .combine = rbind) %do% {
  
  message(paste("Processing", conf$name))
  
  row_results <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
    # Calls the primary experiment function for Setting III 
    # Ensure the run_setting3_expt function is defined in the global environment
    run_setting3_expt(true_T = conf$T, 
                      theta0_star = conf$theta0_star, 
                      theta1_star = conf$theta1_star)
  }
  
  # Summarize the 500 runs for this specific configuration 
  row_results %>%
    summarise(
      T = conf$T,
      Theta0_Bound = conf$theta0_star,
      Theta1_Bound = conf$theta1_star,
      Univ_Coverage = mean(Univ_Cov), # Target ~0.9 
      Adapt_Coverage = mean(Adapt_Cov),
      Univ_Size = mean(Univ_Sz),
      Adapt_Size = mean(Adapt_Sz), # Typically smaller than Universal 
      Avg_Abs_Dev = mean(T_hat_error),
      Avg_Delay = mean(Delay)
    )
}

stopCluster(cl)

# --- Final Formatted Output ---
print("--- TABLE 4 REPRODUCTION: SETTING III (ALL ROWS) ---")
# Column names match the paper's metrics 
print(final_table_4)

# Figure 2(b) Visualization
viz_runs_4 <- bind_rows(lapply(1:5, function(i) {
  d <- run_setting3_expt(100, 0.1, 0.9, is_viz = TRUE) # Params for Fig 2b 
  d$Run <- paste("Run", i); d
}))


p=ggplot(viz_runs_4, aes(x = Time, y = Data)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 100, color = "black", size = 1) + # True T
  geom_point(data = subset(viz_runs_4, In_CI), aes(x = Time, y = 0), color = "red", shape = 16, size = 1.5) +
  facet_wrap(~Run, ncol = 1, scales = "fixed") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), strip.text = element_blank()) +
  labs(title = "",
       x = "Time", y = "Data")+
  theme_minimal() +
  xlim(0, 130) +            # Set x-axis limits
  scale_x_continuous(breaks = seq(0, 130, by = 10), 
                     minor_breaks = seq(0, 130, by = 2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank()   # Customize facet label size if desired
  )
ggsave("fig2b.png", plot = p, width = 8, height = 5, dpi = 300)



