library(dplyr)
library(foreach)
library(doParallel)

run_aue_kirch_detector <- function(data, threshold, historic_data) {
  # 1. Estimation from historic (training) data
  m <- length(historic_data)
  x_bar_m <- mean(historic_data)
  sigma_hat <- sd(historic_data)
  
  # 2. Monitoring phase
  S <- 0 # This will track the partial sum of deviations
  
  for (i in seq_along(data)) {
    # Update the partial sum Psi(m, i)
    # The paper uses: sum_{j=1 to i} (X_{m+j} - x_bar_m)
    S <- S + (data[i] - x_bar_m)
    
    # Apply the weight function w_gamma with gamma = 0
    # Formula: w(m, i) = (1 + i/m)^-1
    weight <- 1 / (1 + i/m)
    
    # Calculate the detector statistic T(m, i)
    # Note: We use the absolute value for a two-sided change (0 to 1 or 0 to -1)
    # If strictly monitoring for a mean increase, remove the abs()
    statistic <- (1 / (sigma_hat * sqrt(m))) * weight * abs(S)
    
    # Check against the boundary (threshold)
    # In the paper's framework, 'threshold' corresponds to the critical value c_alpha
    if (statistic >= threshold) {
      return(i)
    }
  }
  
  return(NA)
}

theta_grid_post <- seq(0.9, 0.9 + 9 * 0.2, by = 0.2) 
w_post <- c(exp(-(0:8)/2) - exp(-(1:9)/2), exp(-4.5))

# Weight setup for Theta0 (Pre-change) 
theta_grid_pre <- seq(0.1, 0.1 - 9 * 0.2, by = -0.2)
w_pre <- c(exp(-(0:8)/2) - exp(-(1:9)/2), exp(-4.5))

# --- 2. Core Functions ---


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

run_setting3_expt <- function(true_T, theta0_star, theta1_star, alpha_univ, A, is_viz = FALSE) {
  # Generate true data: N(0,1) then N(1,1) 
  hdata <- rnorm(200,0,1)
  obs <- c(rnorm(true_T - 1, 0, 1), rnorm(400, 1, 1))
  tau <- run_aue_kirch_detector(obs, A, hdata)
  if (is.na(tau) || tau < true_T) return(NULL)
  
  data_tau <- obs[1:tau]
  
 
  univ_set <- c()
  
  for (t in 1:tau) {
    mt_obs <- calc_Mt_setting3(data_tau, t, tau, theta_grid_pre, w_pre, theta_grid_post, w_post, theta0_star, theta1_star)
    
    # 1. Universal Method (Equation 14) 
    if (mt_obs < 2 / (alpha_univ))
      univ_set <- c(univ_set, t)
  }
  
  if (is_viz) {
    return(data.frame(Time = 1:tau, Data = data_tau, In_CI = (1:tau) %in% univ_set, Run = 1))
  } else {
    return(data.frame(T = true_T, Delay = tau - true_T, 
                      Univ_Sz = length(univ_set), 
                      Univ_Cov = (true_T %in% univ_set)
                     ))
  }
}

# --- 4. Execution ---

# Parallel setup for Table 4 
cl <- makeCluster(parallel::detectCores()-1); registerDoParallel(cl)
table_data_4 <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  # Using first row params: T=100, Theta0=[0, 0.25], Theta1=[0.75, Inf) 
  run_setting3_expt(100, 0.25, 0.75, 0.1, 3)
}
stopCluster(cl)

# Summary Output for Table 4 
print("--- TABLE 4 REPRODUCTION: SETTING III (T=100) ---")
table_4_summary <- table_data_4 %>% 
  summarise(Univ_Coverage = mean(Univ_Cov),
            Univ_Size = mean(Univ_Sz),
            Avg_Delay = mean(Delay))
print(table_4_summary)

###

cl <- makeCluster(parallel::detectCores()-1); registerDoParallel(cl)
table_data_4 <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  # Using first row params: T=100, Theta0=[0, 0.25], Theta1=[0.75, Inf) 
  run_setting3_expt(100, 0.25, 0.75, 0.1, 2)
}
stopCluster(cl)

# Summary Output for Table 4 
print("--- TABLE 4 REPRODUCTION: SETTING III (T=100) ---")
table_4_summary <- table_data_4 %>% 
  summarise(Univ_Coverage = mean(Univ_Cov),
            Univ_Size = mean(Univ_Sz), 
             Avg_Delay = mean(Delay))
print(table_4_summary)

###

cl <- makeCluster(parallel::detectCores()-1); registerDoParallel(cl)
table_data_4 <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  # Using first row params: T=100, Theta0=[0, 0.25], Theta1=[0.75, Inf) 
  run_setting3_expt(500, 0.25, 0.75, 0.1, 3)
}
stopCluster(cl)

# Summary Output for Table 4 
print("--- TABLE 4 REPRODUCTION: SETTING III (T=100) ---")
table_4_summary <- table_data_4 %>% 
  summarise(Univ_Coverage = mean(Univ_Cov),
            Univ_Size = mean(Univ_Sz),
            Avg_Delay = mean(Delay))
print(table_4_summary)

###

cl <- makeCluster(parallel::detectCores()-1); registerDoParallel(cl)
table_data_4 <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  # Using first row params: T=100, Theta0=[0, 0.25], Theta1=[0.75, Inf) 
  run_setting3_expt(500, 0.25, 0.75, 0.1, 2)
}
stopCluster(cl)

# Summary Output for Table 4 
print("--- TABLE 4 REPRODUCTION: SETTING III (T=100) ---")
table_4_summary <- table_data_4 %>% 
  summarise(Univ_Coverage = mean(Univ_Cov),
            Univ_Size = mean(Univ_Sz),
            Avg_Delay = mean(Delay))
print(table_4_summary)