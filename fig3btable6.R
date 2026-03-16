library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Core Functions ---

# e-detector for Setting V (Equation 612)
run_e_detector_V <-function(data, threshold) {
    S <- 0
    for (i in seq_along(data)) {
      if(i==1)
        mu=0
      else 
        mu=min(0,mean(data[1:(i-1)]))
      lr_inc <- dnorm(data[i], mu, 1, log = TRUE) - dnorm(data[i], 0.75, 1, log = TRUE)
      S <- max(0, S + lr_inc)
      if (exp(S) >= threshold) return(i)
    }
    return(NA)
}

calc_Mt_V <- function(data, t, tau) {
  if (t > tau) return(-Inf)
  
  # Forward t-delay e-process (R_n): lambda_i = max(0.75, mean_past) [cite: 187, 614]
  fwd_data <- data[t:tau]
  log_R_t <- 0
  if (length(fwd_data) > 1) {
    for (i in 2:length(fwd_data)) {
      lambda_i <- max(0.75, mean(fwd_data[1:(i-1)]))
      log_R_t <- log_R_t -(-lambda_i * fwd_data[i-1] + lambda_i^2 / 2)
    }
  }
  
  # Backward t-delay e-process (S_n): mu_i = min(0, mean_future) - 0.75 [cite: 188, 614]
  log_S_t <- 0
  if (t > 2) { # Ensure segment is long enough for predictable plug-in
    bwd_data <- data[1:(t-1)]
    s_vals <- c(0)
    s_acc <- 0
    # Predictable mu_i requires at least one data point to the 'right'
    for (i in (t-1):2) {
      mu_i <- min(0, mean(bwd_data[i:(t-1)])) 
      s_acc <- s_acc + (-(0.75-mu_i) * bwd_data[i-1] + (0.25-mu_i^2) / 2)
      s_vals <- c(s_vals, s_acc)
    }
    log_S_t <- max(s_vals)
  }
  
  return(exp(max(log_R_t, log_S_t)))
}

# --- 2. Simulation Runner ---


run_setting5_combined <- function(true_T, iterations, is_viz = FALSE) {
  alpha <- 0.1
  A <- 100
  N_sims <- 100
  
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  results <- tryCatch({
    # CRITICAL: Exporting functions to workers via .export 
    foreach(i = 1:iterations, .combine = rbind, 
            .packages = c("stats", "dplyr"),
            .export = c("run_e_detector_V", "calc_Mt_V")) %dopar% {
              
              # True Data: N(1,1) -> U[-1.2, 0.8] [cite: 605]
              obs <- c(rnorm(true_T - 1, 1, 1), runif(300, -1.2, 0.8))
              tau <- run_e_detector_V(obs, A)
              print(tau)
              if (is.na(tau) || tau < true_T) return(NULL) # Conditional on tau >= T [cite: 63]
              
              data_tau <- obs[1:tau]
              
              # Estimate rt* using P0* = N(0.75, 1) (Assumption 1) [cite: 613]
              null_stops <- replicate(N_sims, {
                s <- run_e_detector_V(rnorm(tau + 1, 0.75, 1), A)
                if(is.na(s)) tau + 1 else s
              })
              
              in_ci <- logical(tau)
              for (t in 1:tau) {
                rt_star <- (1 + sum(null_stops >= t)) / (N_sims + 1)
                mt <- calc_Mt_V(data_tau, t, tau)
                # Universal threshold check [cite: 271]
                if (mt < 2 / (alpha * rt_star)) in_ci[t] <- TRUE
              }
              
              if (is_viz) {
                return(data.frame(Time = 1:tau, Data = data_tau, In_CI = in_ci, Run = i))
              } else {
                return(data.frame(T = true_T, Delay = tau - true_T, Size = sum(in_ci),
                                  Covered = (true_T %in% which(in_ci))))
              }
            }
  }, finally = { stopCluster(cl) })
  
  return(results)
}

# --- 3. Execution ---

# Table 6 Metrics
cat("Calculating Table 6 metrics...\n")
table_res <- run_setting5_combined(100, iterations = 500)
print("--- TABLE 6 (SETTING V) ---")
print(table_res %>% summarise(T = 100, Coverage = mean(Covered), Size = mean(Size), Delay = mean(Delay)))

table_res <- run_setting5_combined(500, iterations = 500)
print("--- TABLE 6 (SETTING V) ---")
print(table_res %>% summarise(T = 500, Coverage = mean(Covered), Size = mean(Size), Delay = mean(Delay)))

# Figure 3(b) Plot
cat("\nGenerating Figure 3(b)...\n")
viz_runs_3<- run_setting5_combined(100, iterations = 5, is_viz = TRUE)

p=ggplot(viz_runs_3, aes(x = Time, y = Data)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 100, color = "black", size = 1) + # True T
  geom_point(data = subset(viz_runs_3, In_CI), aes(x = Time, y = 0), color = "red", shape = 16, size = 1.5) +
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
ggsave("fig3b.png", plot = p, width = 8, height = 5, dpi = 300)


