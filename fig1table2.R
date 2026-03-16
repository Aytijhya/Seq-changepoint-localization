# Reproduction of Figure 1 & Table 2
# Setting I: Known pre- and post-change distributions
library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Parameter Setup ---

alpha <- 0.1
A <- 1000      # CUSUM threshold
B <- 100        # Simulations for Adaptive Quantile
N_sims <- 50   # Simulations for r_t estimation
iterations <- 500

# --- 2. Core Logic Functions ---

# CUSUM Detector [cite: 526]
run_cusum <- function(data, threshold) {
  S <- 0
  for (i in seq_along(data)) {
    lr_inc <- dnorm(data[i], 1, 1, log = TRUE) - dnorm(data[i], 0, 1, log = TRUE)
    S <- max(0, S + lr_inc)
    if (exp(S) >= threshold) return(i)
  }
  return(NA)
}

# Universal Test Statistic Mt (Equation 9) [cite: 201]
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
run_single_sim <- function(true_T, is_viz = FALSE) {
  obs <- c(rnorm(true_T - 1, 0, 1), rnorm(500, 1, 1))
  tau <- run_cusum(obs, A)
  if (is.na(tau) || tau < true_T) return(NULL) 
  
  data_at_tau <- obs[1:tau]
  null_stops <- replicate(N_sims, {
    s <- run_cusum(rnorm(tau + 2, 0, 1), A)
    ifelse(is.na(s), tau + 2, s)
  })
  
  univ_set <- c(); adapt_set <- c()
  
  for (t in 1:tau) {
    rt <- (1 + sum(null_stops >= t)) / (N_sims + 1)
    mt <- calc_Mt(data_at_tau, t, tau)
    if (mt < 2 / (alpha * rt)) univ_set <- c(univ_set, t)
    
    m_t_sims <- replicate(B, {
      sim_obs <- c(rnorm(t-1, 0, 1), rnorm(200, 1, 1))
      sim_tau <- run_cusum(sim_obs, A)
      if(is.na(sim_tau)) return(0)
      calc_Mt(sim_obs[1:sim_tau], t, sim_tau)
    })
    if (mt <= quantile(m_t_sims, 1 - alpha * rt)) adapt_set <- c(adapt_set, t)
  }
  
  res <- list(stats = data.frame(T = true_T, Delay = tau - true_T, 
                                 Univ_Size = length(univ_set), Adapt_Size = length(adapt_set),
                                 Univ_Cover = (true_T %in% univ_set), Adapt_Cover = (true_T %in% adapt_set)))
  if(is_viz) {
    res$plot_df <- data.frame(Time = 1:tau, Data = data_at_tau, 
                              In_CI = (1:tau) %in% adapt_set)
  }
  return(res)
}

# --- 4. Parallel Execution for Table 2 ---
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores); registerDoParallel(cl)

table_results <- foreach(i = 1:iterations, .combine = rbind, .packages = c("stats", "dplyr")) %dopar% {
  r1 <- run_single_sim(100)$stats
  r2 <- run_single_sim(500)$stats
  rbind(r1, r2)
}
stopCluster(cl)

# --- 5. Generate Figure 1 Visuals (5 runs) ---
viz_results100 <- lapply(1:5, function(i) {
  out <- run_single_sim(100, is_viz = TRUE)
  out$plot_df$Run <- paste("Run", i)
  out$plot_df
}) %>% bind_rows()

viz_results500 <- lapply(1:5, function(i) {
  out <- run_single_sim(500, is_viz = TRUE)
  out$plot_df$Run <- paste("Run", i)
  out$plot_df
}) %>% bind_rows()

# --- 6. Printing Results ---

# Table 2 Summary 
print("--- REPRODUCED TABLE 2: SETTING I ---")
summary_tab <- table_results %>% group_by(T) %>%
  summarise(Coverage_Univ = mean(Univ_Cover), Coverage_Adapt = mean(Adapt_Cover),
            Size_Univ = mean(Univ_Size), Size_Adapt = mean(Adapt_Size), Delay = mean(Delay))
print(summary_tab)

# Figure 1 Plot 

p=ggplot(viz_results100, aes(x = Time, y = Data)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 100, color = "black", size = 1) + # True T
  geom_point(data = subset(viz_results100, In_CI), aes(x = Time, y = 0), color = "red", shape = 16, size = 1.5) +
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
ggsave("fig1a.png", plot = p, width = 8, height = 5, dpi = 300)






p=ggplot(viz_results500, aes(x = Time, y = Data)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 500, color = "black", size = 1) + # True T
  geom_point(data = subset(viz_results500, In_CI), aes(x = Time, y = 0), color = "red", shape = 16, size = 1.5) +
  facet_wrap(~Run, ncol = 1, scales = "fixed") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), strip.text = element_blank()) +
  labs(title = "",
       x = "Time", y = "Data")+
  theme_minimal() +
  xlim(0, 530) +            # Set x-axis limits
  scale_x_continuous(breaks = seq(0, 530, by = 50), 
                     minor_breaks = seq(0, 530, by = 5)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank()   # Customize facet label size if desired
  )
ggsave("fig1b.png", plot = p, width = 8, height = 5, dpi = 300)






