library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Parameter Setup ---

alpha <- 0.1
A <- 1000
N_sims <- 100
iterations <- 50
mu <- 0.25

# Numerically solve for lambda* (Equation: e^l = (1+l(1-mu))/(1-l*mu))
find_lambda_star <- function(m) {
  root_fn <- function(l) exp(l) - (1 + l * (1 - m)) / (1 - l * m)
  uniroot(root_fn, lower = 0.01, upper = 1/m - 0.01)$root
}
lambda_star <- find_lambda_star(mu) # Approx 1.705

# --- 2. Core Methodological Functions ---

# Histogram-based p_hat estimator (Predictable plug-in)
get_p_hat <- function(segment, x_val) {
  bin_idx <- pmin(pmax(ceiling(x_val * 10), 1), 10)
  k_counts <- if(length(segment) == 0) 0 else sum(pmin(pmax(ceiling(segment * 10), 1), 10) == bin_idx)
  # Formula: (1 + sum 1(x in Bk) * sum 1(Xj in Bk)) / (9 + t - i)
  return((1 + k_counts) / (9 + length(segment)))
}

# e-detector (Shin et al., 2023) using the recursive logic
run_e_detector <- function(data, threshold) {
  n <- length(data)
  for (i in 1:n) {
    # Check max log-e-process over all starting points j
    log_e_vals <- sapply(1:i, function(j) {
      segment <- data[j:i]
      # Predictable histogram estimator p_hat'
      s_log_e <- 0
      if(length(segment) > 1) {
        for(k in 2:length(segment)) {
          p_val <- get_p_hat(segment[1:(k-1)], segment[k])
          s_log_e <- s_log_e + log(p_val / 0.1) # p0(x) = 1/10
        }
      }
      return(s_log_e)
    })
    if (exp(max(log_e_vals)) >= threshold) return(i)
  }
  return(NA)
}

# Universal Mt Statistic for Setting IV (Equation 9)
calc_Mt_IV <- function(data, t, tau) {
  if (t > tau) return(-Inf)
  # Forward t-delay e-process (Numeraire e-variable)
  fwd_data <- data[t:tau]
  log_R_t <- max(cumsum(log(1 + lambda_star * (fwd_data - mu))))
  
  # Backward t-delay e-process (Predictable plug-in)
  log_S_t <- 0
  if (t > 1) {
    bwd_data <- data[1:(t-1)]
    log_S_vals <- numeric(t - 1)
    s_acc <- 0
    for (i in (t-1):2) {
      p_val <- get_p_hat(bwd_data[i:(t-1)], bwd_data[i-1])
      s_acc <- s_acc + log(p_val / 0.1)
      log_S_vals[i-1] <- s_acc
    }
    log_S_t <- max(log_S_vals)
  }
  return(exp(max(log_R_t, log_S_t)))
}

# --- 3. Experiment Runner ---

run_setting4_full <- function(true_T, d_type, is_viz = FALSE) {
  # True post-change density f1(x)
  p1_data <- if(d_type == "beta") rbeta(600, 1, 4) else 
    sample(c(runif(480, 0, 0.2), runif(120, 0.2, 1)))
  
  obs <- c(runif(true_T - 1, 0, 1), p1_data)
  tau <- run_e_detector(obs, A)
  
  if (is.na(tau) || tau < true_T) return(NULL) # Conditional on tau >= T
  
  data_tau <- obs[1:tau]
  null_stops <- replicate(N_sims, {
    s <- run_e_detector(runif(tau + 200, 0, 1), A)
    if(is.na(s)) tau + 200 else s
  })
  
  in_ci <- logical(tau)
  for (t in 1:tau) {
    rt <- (1 + sum(null_stops >= t)) / (N_sims + 1)
    mt <- calc_Mt_IV(data_tau, t, tau)
    if (mt < 2 / (alpha * rt)) in_ci[t] <- TRUE
  }
  
  if (is_viz) {
    return(data.frame(Time = 1:tau, Data = data_tau, In_CI = in_ci))
  } else {
    ci_idx <- which(in_ci)
    return(data.frame(T = true_T, Delay = tau - true_T, Size = length(ci_idx),
                      Covered = (true_T %in% ci_idx), Error = abs(median(ci_idx) - true_T)))
  }
}

# --- 4. Execution ---

cl <- makeCluster(parallel::detectCores() - 1); registerDoParallel(cl)

# Table 5 Metrics (Row 1: T=100, 4(1-x)^3)
table_res <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr"),
                     .export = c("run_e_detector", "get_p_hat", "calc_Mt_IV", "mu", "lambda_star", "A")) %dopar% {
                       run_setting4_full(100, "beta")
                     }

# Visualization for Figure 3(a) (5 runs)
viz_data_list <- list()
while(length(viz_data_list) < 5) {
  res <- run_setting4_full(100, "beta", is_viz = TRUE)
  if(!is.null(res)) viz_data_list[[length(viz_data_list) + 1]] <- res
}
stopCluster(cl)

# --- 5. Output Results ---

print("--- TABLE 5 (SETTING IV PARTIAL REPRODUCTION) ---")
print(table_res %>% summarise(T = 100, Coverage = mean(Covered), Size = mean(Size), 
                              Abs_Dev = mean(Error), Delay = mean(Delay)))

viz_df <- bind_rows(lapply(1:5, function(i) { viz_data_list[[i]]$Run <- paste("Run", i); viz_data_list[[i]] }))
ggplot(viz_df, aes(x = Time, y = Data)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_point(data = filter(viz_df, In_CI), color = "red", size = 0.8) +
  geom_vline(xintercept = 100, color = "black", size = 1) +
  facet_wrap(~Run, ncol = 1, scales = "free_y") + theme_bw() +
  labs(title = "Figure 3(a) Reproduction: Setting IV",
       subtitle = "Red points = Universal Confidence Set (10)")