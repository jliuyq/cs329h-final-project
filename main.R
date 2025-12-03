## ---------------------------------------------------------------------------
## File: main.R
## Author: Yiqing Liu
## Description: Simulation of Adaptive Pairwise Comparison using BT Models.
##              Compares Random, Entropy, and CAT selection policies.
## ---------------------------------------------------------------------------

# --- Dependencies -----------------------------------------------------------
library(dplyr)

# --- 1. Data Generation & Helper Functions ----------------------------------

#' Generate True Item Strengths
#'
#' Generates latent parameters for items based on IID or Clustered design.
#'
#' @param n_items Integer. Number of items.
#' @param design String. "iid" (standard normal) or "clustered" (mixture).
#' @param sd_iid Numeric. SD for IID design.
#' @param cluster_centers Numeric vector. Means for clusters.
#' @param cluster_sd Numeric. SD within clusters.
#' @return Numeric vector of centered true theta values.
generate_theta <- function(n_items,
                           design = c("iid", "clustered"),
                           sd_iid = 1,
                           cluster_centers = c(-1.5, -0.5, 0.5, 1.5),
                           cluster_sd = 0.2) {
  design <- match.arg(design)
  
  if (design == "iid") {
    theta <- rnorm(n_items, mean = 0, sd = sd_iid)
  } else {
    G <- length(cluster_centers)
    sizes <- rep(floor(n_items / G), G)
    remainder <- n_items - sum(sizes)
    if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1
    
    theta <- numeric(0)
    for (g in seq_len(G)) {
      theta <- c(theta, rnorm(sizes[g], mean = cluster_centers[g], sd = cluster_sd))
    }
  }
  theta - mean(theta) # Center to mean zero
}

#' Create Item Pairs
#'
#' @param n_items Integer. Number of items.
#' @return Matrix with columns "i" and "j" representing all unique pairs.
make_pairs <- function(n_items) {
  pairs <- t(combn(n_items, 2))
  colnames(pairs) <- c("i", "j")
  storage.mode(pairs) <- "integer"
  pairs
}

#' Calculate Bradley-Terry Probabilities
#' 
#' @param theta Numeric vector. Latent strengths.
#' @param pairs Matrix. Item pairs.
#' @return Numeric vector of probabilities that i beats j.
bt_probs <- function(theta, pairs) {
  eta <- theta[pairs[, "i"]] - theta[pairs[, "j"]]
  plogis(eta)
}

# --- 2. Model Fitting (Newton-Raphson MAP) ----------------------------------

#' Log-Posterior for BT Model
logpost_bt <- function(theta, data, tau2) {
  if (nrow(data) == 0) return(-0.5 * sum(theta^2) / tau2) # Prior only
  
  i <- data$i; j <- data$j; y <- data$y
  eta <- theta[i] - theta[j]
  p <- plogis(eta)
  
  loglik <- sum(y * log(p + 1e-10) + (1 - y) * log(1 - p + 1e-10)) # Added epsilon
  logprior <- -0.5 * sum(theta^2) / tau2
  loglik + logprior
}

#' Gradient of Log-Posterior
grad_logpost_bt <- function(theta, data, tau2) {
  n_items <- length(theta)
  grad <- -theta / tau2 # Prior gradient
  
  if (nrow(data) > 0) {
    i <- data$i; j <- data$j; y <- data$y
    p <- plogis(theta[i] - theta[j])
    r <- y - p
    
    # Vectorized accumulation of gradients would be faster, but loop is clear
    for (k in seq_along(r)) {
      grad[i[k]] <- grad[i[k]] + r[k]
      grad[j[k]] <- grad[j[k]] - r[k]
    }
  }
  grad
}

#' Hessian of Log-Posterior (Negative Precision Matrix)
hessian_bt <- function(theta, data, tau2) {
  n_items <- length(theta)
  H <- diag(1 / tau2, n_items) # Start with prior precision
  
  if (nrow(data) > 0) {
    i <- data$i; j <- data$j
    p <- plogis(theta[i] - theta[j])
    w <- p * (1 - p)
    
    for (k in seq_along(w)) {
      ii <- i[k]; jj <- j[k]; wk <- w[k]
      H[ii, ii] <- H[ii, ii] + wk
      H[jj, jj] <- H[jj, jj] + wk
      H[ii, jj] <- H[ii, jj] - wk
      H[jj, ii] <- H[jj, ii] - wk
    }
  }
  H
}

#' Fit BT Model via MAP
#'
#' @param data Dataframe. Columns i, j, y.
#' @param n_items Integer.
#' @param tau2 Numeric. Prior variance.
#' @param init Numeric vector. Initial theta guess.
#' @return List containing theta_hat, Sigma, and convergence info.
bt_map_fit <- function(data, n_items, tau2 = 5, init = NULL) {
  if (is.null(init)) init <- rep(0, n_items)
  
  # Optimization
  f <- function(theta) -logpost_bt(theta, data, tau2)
  g <- function(theta) -grad_logpost_bt(theta, data, tau2)
  
  opt <- optim(par = init, fn = f, gr = g, method = "BFGS", control = list(maxit = 1000))
  
  theta_hat <- opt$par - mean(opt$par) # Center
  H <- hessian_bt(theta_hat, data, tau2)
  
  # Numerical stability check for inversion
  Sigma <- tryCatch(solve(H), error = function(e) diag(1/diag(H)))
  
  list(theta_hat = theta_hat, Sigma = Sigma, convergence = opt$convergence)
}

# --- 3. Selection Policies --------------------------------------------------

select_pair_random <- function(theta_hat, Sigma, pairs, asked) {
  candidates <- which(!asked)
  if (length(candidates) == 0L) stop("No unasked pairs left.")
  if (length(candidates) == 1L) return(candidates)
  sample(candidates, size = 1L)
}

select_pair_entropy <- function(theta_hat, Sigma, pairs, asked) {
  p <- bt_probs(theta_hat, pairs)
  score <- p * (1 - p)
  score[asked] <- -Inf
  
  best <- which(score == max(score, na.rm = TRUE))
  if (length(best) > 1) sample(best, size = 1L) else best
}

contrast_var <- function(Sigma, pairs) {
  d <- diag(Sigma)
  v <- d[pairs[, "i"]] + d[pairs[, "j"]] - 2 * Sigma[cbind(pairs[, "i"], pairs[, "j"])]
  pmax(v, 0)
}

select_pair_cat <- function(theta_hat, Sigma, pairs, asked, top_frac = 0.2) {
  if (all(asked)) stop("No unasked pairs left.")
  
  p <- bt_probs(theta_hat, pairs)
  entropy <- p * (1 - p)
  entropy[asked] <- -Inf
  
  # Two-stage selection: Filter by high entropy, then max contrast variance
  n_unasked <- sum(!asked)
  n_keep <- max(1L, floor(top_frac * n_unasked))
  
  ord <- order(entropy, decreasing = TRUE)
  ord_unasked <- ord[!asked[ord]]
  shortlist <- ord_unasked[seq_len(n_keep)]
  
  v <- contrast_var(Sigma, pairs)
  u <- entropy * v
  u[asked] <- -Inf
  
  shortlist[which.max(u[shortlist])]
}

# --- 4. Metrics & Simulation ------------------------------------------------

#' Compute Evaluation Metrics
#'
#' Calculates correlation and ranking accuracy metrics.
#'
#' @param theta_hat Numeric vector. Estimated parameters.
#' @param theta_true Numeric vector. True parameters.
#' @param k_vec Vector of integers. Top-k levels to evaluate.
#' @return Named vector of metrics.
compute_metrics <- function(theta_hat, theta_true, k_vec = c(4, 8)) {
  # Correlations
  ken <- cor(theta_hat, theta_true, method = "kendall")
  spr <- cor(theta_hat, theta_true, method = "spearman")
  
  # RMSE
  rmse <- sqrt(mean((theta_hat - theta_true)^2))
  
  # Top-K Accuracy
  metrics <- c(kendall = ken, spearman = spr, rmse = rmse)
  
  for (k in k_vec) {
    top_true <- order(theta_true, decreasing = TRUE)[1:k]
    top_hat  <- order(theta_hat, decreasing = TRUE)[1:k]
    acc <- length(intersect(top_true, top_hat)) / k
    names(acc) <- paste0("top", k, "_acc")
    metrics <- c(metrics, acc)
  }
  
  metrics
}

simulate_bt_adaptive <- function(theta_true, pairs_all, T_budget, policy, tau2 = 5) {
  n_items <- length(theta_true)
  n_pairs <- nrow(pairs_all)
  asked <- rep(FALSE, n_pairs)
  
  data_obs <- data.frame(i = integer(0), j = integer(0), y = integer(0))
  kendall <- numeric(T_budget)
  
  # Initial fit (prior only)
  fit <- bt_map_fit(data_obs, n_items = n_items, tau2 = tau2)
  theta_hat <- fit$theta_hat
  Sigma <- fit$Sigma
  
  for (t in seq_len(T_budget)) {
    # 1. Select Pair
    k <- switch(policy,
                random  = select_pair_random(theta_hat, Sigma, pairs_all, asked),
                entropy = select_pair_entropy(theta_hat, Sigma, pairs_all, asked),
                cat     = select_pair_cat(theta_hat, Sigma, pairs_all, asked))
    
    asked[k] <- TRUE
    i_curr <- pairs_all[k, "i"]; j_curr <- pairs_all[k, "j"]
    
    # 2. Simulate Outcome
    p_true <- plogis(theta_true[i_curr] - theta_true[j_curr])
    y_new  <- rbinom(1L, size = 1L, prob = p_true)
    data_obs <- rbind(data_obs, data.frame(i = i_curr, j = j_curr, y = y_new))
    
    # 3. Refit
    fit <- bt_map_fit(data_obs, n_items, tau2, init = theta_hat)
    theta_hat <- fit$theta_hat
    Sigma <- fit$Sigma
    
    # 4. Track Performance
    kendall[t] <- cor(theta_hat, theta_true, method = "kendall")
  }
  
  list(kendall = kendall, theta_hat = theta_hat, data_obs = data_obs)
}

run_replications <- function(n_items, T_budget, design, n_rep, tau2, policies, k_vec) {
  pairs_all <- make_pairs(n_items)
  curves_list <- list()
  metrics_list <- list()
  
  for (rep in seq_len(n_rep)) {
    theta_true <- generate_theta(n_items, design = design)
    
    for (pol in policies) {
      sim <- simulate_bt_adaptive(theta_true, pairs_all, T_budget, pol, tau2)
      
      # Store Curve
      curves_list[[length(curves_list) + 1]] <- data.frame(
        rep = rep, t = seq_len(T_budget), policy = pol, kendall = sim$kendall
      )
      
      # Store Final Metrics
      mets <- compute_metrics(sim$theta_hat, theta_true, k_vec = k_vec)
      metrics_list[[length(metrics_list) + 1]] <- data.frame(
        rep = rep, policy = pol, t_final = T_budget, as.list(mets)
      )
    }
  }
  list(curves = do.call(rbind, curves_list), metrics = do.call(rbind, metrics_list))
}

run_experiment_grid <- function(n_vec, design_vec, mult_vec, n_rep, tau2, policies) {
  all_curves <- list()
  all_metrics <- list()
  
  for (design in design_vec) {
    for (n_items in n_vec) {
      base_T <- ceiling(n_items * log(n_items))
      for (m in mult_vec) {
        T_budget <- m * base_T
        message(sprintf("Running: %s, n=%d, T=%d", design, n_items, T_budget))
        
        res <- run_replications(n_items, T_budget, design, n_rep, tau2, policies, c(4, 8))
        
        res$curves$design <- design; res$curves$n_items <- n_items; res$curves$mult <- m
        res$metrics$design <- design; res$metrics$n_items <- n_items; res$metrics$mult <- m
        
        all_curves[[length(all_curves) + 1]] <- res$curves
        all_metrics[[length(all_metrics) + 1]] <- res$metrics
      }
    }
  }
  list(curves = do.call(rbind, all_curves), metrics = do.call(rbind, all_metrics))
}

# --- 5. Main Execution Block ------------------------------------------------

if (sys.nframe() == 0) {
  set.seed(123)
  
  # Configure Experiment
  # NOTE: Reduced n_rep for quick testing. Increase for final paper.
  exp_config <- list(
    n_vec      = c(20, 40, 80),        # number of items     
    design_vec = c("iid", "clustered"),     # two different types of generation
    mult_vec   = c(1, 2, 3),                # budget
    n_rep      = 50                 # n_replications
  )
  
  message("Starting Experiment...")
  results <- run_experiment_grid(
    n_vec      = exp_config$n_vec,
    design_vec = exp_config$design_vec,
    mult_vec   = exp_config$mult_vec,
    n_rep      = exp_config$n_rep,
    tau2       = 5,
    policies   = c("random", "entropy", "cat")
  )
  
  # Summary
  summary_metrics <- results$metrics %>%
    group_by(design, n_items, mult, T_final=T, policy) %>%
    summarise(
      n_rep = n(),
      across(c(kendall, spearman, rmse, top4_acc, top8_acc),
             list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  print(head(summary_metrics))
  
  # Output
  output_file <- "bt_cat_summary_metrics.csv"
  write.csv(summary_metrics, output_file, row.names = FALSE)
  message(paste("Results saved to:", output_file))
}