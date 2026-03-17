library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)

# --- 1. Core Huber-Robust Functions ---



# Huber-robust e-detector (Source: Equation 24 and Page 34)
# Uses the clipped likelihood ratio pi(x) for epsilon-contamination
run_huber_detector <- function(data, threshold, eps) {
  # Constants for Huber LFD (approximate for N(0,1) vs N(1,1))
  # Source: Huber (1965); Saha and Ramdas (2026)
  lexpt=function(l,eps,mu){
    if(mu>0)
      return((l*pnorm(log(l)/mu+mu/2)+1-pnorm(log(l)/mu-mu/2))*(1-eps)-1)
    else
      return((l-l*pnorm(log(l)/mu+mu/2)+pnorm(log(l)/mu-mu/2))*(1-eps)-1)
  }
  
  uexpt=function(u,eps,mu){
    if(mu>0)
      return((1/u-pnorm(log(u)/mu-mu/2)/u+pnorm(log(u)/mu+mu/2))*(1-eps)-1)
    else
      return((pnorm(log(u)/mu-mu/2)/u+1-
                pnorm(log(u)/mu+mu/2))*(1-eps)-1)
  }
  
  c1 <- uniroot(lexpt,eps=eps,mu=1,lower = 0,upper = 1/eps)$root  # lower clipping
  c2 <- uniroot(uexpt,eps=eps,mu=1,lower = 0.1,upper = 1/eps)$root # upper clipping
  
  e_proc <- 0
  for (i in seq_along(data)) {
    # Standard LR for N(1,1) vs N(0,1)
    lr <- exp(data[i] - 0.5)
    # Huber clipping (clipped LR pi(x))
    pi_x <- min(max(lr, c1), c2)
    
    # e-detector recursion: e_i = pi_x * max(e_{i-1}, 1)
    e_proc <- pi_x * max(e_proc, 1)
    
    if (e_proc >= threshold) return(i)
  }
  return(NA)
}

# Mt Statistic for Setting VI (Source: Equation 14 and Page 34)
calc_Mt_robust <- function(data, t, tau, eps) {
  if (t > tau) return(-Inf)
  lexpt=function(l,eps,mu){
    if(mu>0)
      return((l*pnorm(log(l)/mu+mu/2)+1-pnorm(log(l)/mu-mu/2))*(1-eps)-1)
    else
      return((l-l*pnorm(log(l)/mu+mu/2)+pnorm(log(l)/mu-mu/2))*(1-eps)-1)
  }
  
  uexpt=function(u,eps,mu){
    if(mu>0)
      return((1/u-pnorm(log(u)/mu-mu/2)/u+pnorm(log(u)/mu+mu/2))*(1-eps)-1)
    else
      return((pnorm(log(u)/mu-mu/2)/u+1-
                pnorm(log(u)/mu+mu/2))*(1-eps)-1)
  }
  c1 <- uniroot(lexpt,eps=eps,mu=1,lower = 0,upper = 1/eps)$root  # lower clipping
  c2 <- uniroot(uexpt,eps=eps,mu=1,lower = 0.1,upper = 1/eps)$root # upper clipping
  
  # Forward Segment: Product of clipped LRs
  fwd_data <- data[t:tau]
  lr_fwd <- exp(-fwd_data + 0.5)
  pi_fwd <- pmin(pmax(lr_fwd, c1), c2)
  R_t <- max(cumprod(pi_fwd))
  
  # Backward Segment: Product of clipped inverse LRs
  if (t > 1) {
    bwd_data <- data[1:(t-1)]
    lr_bwd <- exp((bwd_data - 0.5))
    pi_bwd <- pmin(pmax(lr_bwd, 1/c2), 1/c1)
    S_t <- max(cumprod(rev(pi_bwd)))
  } else { S_t <- 0 }
  
  return(max(R_t, S_t))
}

# --- 2. Combined Experiment Wrapper ---

run_setting6_combined <- function(true_T, eps, iterations = 200, is_viz = FALSE) {
  alpha <- 0.1
  A <- 1000
  N_sims <- 100
  
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  
  results <- tryCatch({
    foreach(i = 1:iterations, .combine = rbind, 
            .packages = c("stats", "dplyr"),
            .export = c("run_huber_detector", "calc_Mt_robust")) %dopar% {
              
              # Data Generation: epsilon-contamination with Cauchy noise 
              # P0 = (1-eps)N(0,1) + eps*Cauchy(-1,10)
              # P1 = (1-eps)N(1,1) + eps*Cauchy(-1,10)
              n_total <- true_T + 200
              noise <- rcauchy(n_total, -1, 10)
              is_contam <- runif(n_total) < eps
              
              pre_data <- rnorm(true_T - 1, 0, 1)
              post_data <- rnorm(200, 1, 1)
              
              obs <- c(pre_data, post_data)
              obs[is_contam] <- noise[is_contam]
              
              tau <- run_huber_detector(obs, A, eps)
              if (is.na(tau) || tau < true_T) return(NULL)
              
              data_tau <- obs[1:tau]
              
              # rt* Estimation using Monte Carlo from P0 (N(0,1) as proxy for LFD)
              null_stops <- replicate(N_sims, {
                s <- run_huber_detector(rnorm(tau + 1, 0, 1), A, eps)
                if(is.na(s)) tau + 1 else s
              })
              
              in_ci <- logical(tau)
              for (t in 1:tau) {
                rt_star <- (1 + sum(null_stops >= t)) / (N_sims + 1)
                mt <- calc_Mt_robust(data_tau, t, tau, eps)
                if (mt < 2 / (alpha * rt_star)) in_ci[t] <- TRUE
              }
              
              if (is_viz) {
                return(data.frame(Time = 1:tau, Data = data_tau, In_CI = in_ci, Run = i))
              } else {
                return(data.frame(T = true_T, Eps = eps, Delay = tau - true_T, 
                                  Size = sum(in_ci), Covered = (true_T %in% which(in_ci)),
                                  Abs_Dev = abs(median(which(in_ci)) - true_T)))
              }
            }
  }, finally = { stopCluster(cl) })
  
  return(results)
}

# --- 3. Execution & Summary ---

# TABLE 6 Reproduction 
cat("Simulating TABLE 6 (Setting VI)...\n")
res_100_01 <- run_setting6_combined(100, 0.01, iterations = 500)

print(res_100_01 %>% summarise(T=100, Eps=0.01, Coverage=mean(Covered), Size=mean(Size), Delay=mean(Delay)))

res_500_01 <- run_setting6_combined(500, 0.01, iterations = 500)

print(res_500_01 %>% summarise(T=500, Eps=0.01, Coverage=mean(Covered), Size=mean(Size), Delay=mean(Delay)))

res_100_001 <- run_setting6_combined(100, 0.001, iterations = 500)

print(res_100_001 %>% summarise(T=100, Eps=0.001, Coverage=mean(Covered), Size=mean(Size), Delay=mean(Delay)))

res_500_001 <- run_setting6_combined(500, 0.001, iterations = 500)

print(res_500_001 %>% summarise(T=500, Eps=0.001, Coverage=mean(Covered), Size=mean(Size), Delay=mean(Delay)))

cat("\nGenerating Robust Visualization...\n")
viz_runs_3 <- run_setting6_combined(100, 0.01, iterations = 5, is_viz = TRUE)

p=ggplot(viz_runs_3, aes(x = Time, y = Data)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 100, color = "black", size = 1) + # True T
  geom_point(data = subset(viz_runs_3, In_CI), aes(x = Time, y = 0), color = "red", shape = 16, size = 1.5) +
  facet_wrap(~Run, ncol = 1, scales = "free_y") +
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




