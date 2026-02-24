############################################################
# D-optimal exact design for a linear model with correlated errors
# (known covariance structure)
#
# This script demonstrates a practical workflow in R:
# 1) Define a candidate design space
# 2) Build X(xi) for a proposed exact design xi
# 3) Build Sigma(xi) from a known correlation model (AR(1) or exponential)
# 4) Compute log det M(xi), where M = X' Sigma^{-1} X
# 5) Search for a good exact design using a simple exchange algorithm
# 6) Compare to a baseline design and compute relative D-efficiency
#
# This is designed to be easy to modify for your own model/design space.
############################################################

## ---------------------------
## 0) User settings (EDIT ME)
## ---------------------------

set.seed(123)

# Total number of observations allowed (exact design size)
n_obs <- 8

# Candidate design space (e.g., time points you are allowed to sample)
# Example: integer times from 0 to 20
candidate_points <- 0:20

# Regression model: linear in parameters
# Here we use an intercept + slope + quadratic term as an example.
# Replace this with your own regression basis f(x).
regression_basis <- function(x) {
  cbind(1, x, x^2)
}

# Known covariance structure choice:
# "ar1" = Corr(e_i, e_j) = rho^{|t_i - t_j|}
# "exp" = Corr(e_i, e_j) = exp(-phi * |t_i - t_j|)
cov_model <- "exp"

# Known correlation parameter(s)
rho <- 0.8   # used for ar1
phi <- 0.25  # used for exp

# Observation variance (sigma^2). For D-optimal search over designs, this
# is a constant multiplicative factor in Sigma and often does not change the argmax,
# but we include it for completeness.
sigma2 <- 1

# Search settings
n_restarts <- 20
max_iter_per_restart <- 200


## -----------------------------------------
## 1) Helper functions: X(xi) and Sigma(xi)
## -----------------------------------------

# Build the design matrix X for an exact design xi
# xi is a numeric vector of selected design points, length n_obs
build_X <- function(xi) {
  X <- regression_basis(xi)
  # Basic rank check is useful: if rank-deficient, D-criterion is undefined/poor
  if (qr(X)$rank < ncol(X)) {
    return(NULL)
  }
  X
}

# Build covariance matrix Sigma(xi) from known correlation structure
# xi = selected design points (e.g., times)
build_Sigma <- function(xi, model = c("ar1", "exp"), rho = NULL, phi = NULL, sigma2 = 1) {
  model <- match.arg(model)
  
  # Pairwise distances |x_i - x_j|
  D <- as.matrix(dist(xi, method = "manhattan", diag = TRUE, upper = TRUE))
  
  if (model == "ar1") {
    if (is.null(rho) || rho <= 0 || rho >= 1) stop("For AR(1), rho must be in (0,1).")
    R <- rho ^ D
  } else if (model == "exp") {
    if (is.null(phi) || phi <= 0) stop("For exponential correlation, phi must be > 0.")
    R <- exp(-phi * D)
  }
  
  Sigma <- sigma2 * R
  Sigma
}


## ---------------------------------------------------------
## 2) Numerically stable D-criterion: log det M(xi)
## ---------------------------------------------------------
## M(xi) = X' Sigma^{-1} X
## We avoid explicit matrix inversion and use solves / Cholesky where possible.

logdet_D_criterion <- function(xi,
                               cov_model = c("ar1", "exp"),
                               rho = NULL, phi = NULL, sigma2 = 1) {
  cov_model <- match.arg(cov_model)
  
  X <- build_X(xi)
  if (is.null(X)) return(-Inf)  # rank-deficient X -> invalid for D-optimality
  
  Sigma <- build_Sigma(xi, model = cov_model, rho = rho, phi = phi, sigma2 = sigma2)
  
  # Cholesky factorization of Sigma for numerical stability
  # If Sigma is not positive definite (or numerically singular), reject design
  chol_Sigma <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(chol_Sigma)) return(-Inf)
  
  # Compute Sigma^{-1} X without explicitly forming Sigma^{-1}
  # solve(Sigma, X) is more stable than solve(Sigma) %*% X
  SinvX <- tryCatch(solve(Sigma, X), error = function(e) NULL)
  if (is.null(SinvX)) return(-Inf)
  
  M <- crossprod(X, SinvX)  # X' Sigma^{-1} X
  
  # M must be positive definite to compute log det
  chol_M <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(chol_M)) return(-Inf)
  
  # log det(M) = 2 * sum(log(diag(chol(M))))
  logdetM <- 2 * sum(log(diag(chol_M)))
  as.numeric(logdetM)
}


## ---------------------------------------------------
## 3) Exact design representation and initialization
## ---------------------------------------------------
## We represent an exact design as a vector xi of length n_obs.
## Replication is allowed (same point can appear multiple times), but in many
## correlated settings very close/duplicate points may be poor due to correlation.
## This script allows replication but you can restrict it if needed.

# Random initial design (with replacement = TRUE allows replication)
random_design <- function(candidate_points, n_obs) {
  sample(candidate_points, size = n_obs, replace = TRUE)
}

# Optional: baseline design (e.g., roughly equally spaced)
baseline_design <- function(candidate_points, n_obs) {
  idx <- round(seq(1, length(candidate_points), length.out = n_obs))
  candidate_points[idx]
}


## ----------------------------------------------------------
## 4) Simple exchange algorithm for exact D-optimal design
## ----------------------------------------------------------
## This is a practical "swap one point at a time" local search.
## It is easy to understand and aligns with the article's emphasis that
## exact-design search in correlated settings is algorithmic/computational,
## rather than relying on the usual continuous-design equivalence tools. :contentReference[oaicite:1]{index=1}

exchange_search <- function(candidate_points,
                            n_obs,
                            cov_model = c("ar1", "exp"),
                            rho = NULL, phi = NULL, sigma2 = 1,
                            init_design = NULL,
                            max_iter = 200) {
  cov_model <- match.arg(cov_model)
  
  if (is.null(init_design)) {
    xi <- random_design(candidate_points, n_obs)
  } else {
    xi <- init_design
  }
  
  best_val <- logdet_D_criterion(xi, cov_model, rho, phi, sigma2)
  improved <- TRUE
  iter <- 0
  
  # Track optimization path for diagnostics
  history <- data.frame(iter = 0, logdet = best_val)
  
  while (improved && iter < max_iter) {
    improved <- FALSE
    iter <- iter + 1
    
    # Try replacing each current design point with each candidate point
    for (i in seq_len(n_obs)) {
      current_point <- xi[i]
      
      # Randomize candidate order to reduce deterministic path dependence
      for (cand in sample(candidate_points, length(candidate_points), replace = FALSE)) {
        if (cand == current_point) next
        
        xi_try <- xi
        xi_try[i] <- cand
        
        val_try <- logdet_D_criterion(xi_try, cov_model, rho, phi, sigma2)
        
        # Accept strict improvement (small tolerance to avoid floating-point noise)
        if (is.finite(val_try) && (val_try > best_val + 1e-10)) {
          xi <- xi_try
          best_val <- val_try
          improved <- TRUE
        }
      }
    }
    
    history <- rbind(history, data.frame(iter = iter, logdet = best_val))
  }
  
  list(
    design = xi,
    logdet = best_val,
    history = history
  )
}


## ----------------------------------------------------------
## 5) Multi-start strategy (important for local optima)
## ----------------------------------------------------------
## Correlated-design problems can be multimodal; multiple restarts are a simple,
## practical way to improve trust in the final result.

multi_start_search <- function(candidate_points,
                               n_obs,
                               cov_model = c("ar1", "exp"),
                               rho = NULL, phi = NULL, sigma2 = 1,
                               n_restarts = 10,
                               max_iter = 200) {
  cov_model <- match.arg(cov_model)
  
  best_overall <- NULL
  
  for (r in seq_len(n_restarts)) {
    out <- exchange_search(
      candidate_points = candidate_points,
      n_obs = n_obs,
      cov_model = cov_model,
      rho = rho, phi = phi, sigma2 = sigma2,
      init_design = random_design(candidate_points, n_obs),
      max_iter = max_iter
    )
    
    if (is.null(best_overall) || out$logdet > best_overall$logdet) {
      best_overall <- out
      best_overall$restart <- r
    }
  }
  
  best_overall
}


## ----------------------------------------------------------
## 6) Relative D-efficiency (exact-design comparison)
## ----------------------------------------------------------
## For comparing two designs xi_A and xi_B with p parameters:
## D-efficiency(A relative to B) = exp((logdet(M_A)-logdet(M_B))/p)
## This equals 1 if they are equal on the criterion, <1 if A is worse than B.

relative_D_efficiency <- function(xi_A, xi_B,
                                  cov_model = c("ar1", "exp"),
                                  rho = NULL, phi = NULL, sigma2 = 1,
                                  p = NULL) {
  cov_model <- match.arg(cov_model)
  
  if (is.null(p)) {
    Xtmp <- build_X(xi_A)
    if (is.null(Xtmp)) stop("xi_A produces rank-deficient X; cannot determine p.")
    p <- ncol(Xtmp)
  }
  
  ld_A <- logdet_D_criterion(xi_A, cov_model, rho, phi, sigma2)
  ld_B <- logdet_D_criterion(xi_B, cov_model, rho, phi, sigma2)
  
  if (!is.finite(ld_A) || !is.finite(ld_B)) return(NA_real_)
  
  exp((ld_A - ld_B) / p)
}


## ----------------------------------------------------------
## 7) Run the workflow
## ----------------------------------------------------------

# Baseline (e.g., approximately equally spaced)
xi_baseline <- baseline_design(candidate_points, n_obs)
ld_baseline <- logdet_D_criterion(
  xi_baseline, cov_model = cov_model, rho = rho, phi = phi, sigma2 = sigma2
)

# Search for a good design using multiple random restarts
best <- multi_start_search(
  candidate_points = candidate_points,
  n_obs = n_obs,
  cov_model = cov_model,
  rho = rho, phi = phi, sigma2 = sigma2,
  n_restarts = n_restarts,
  max_iter = max_iter_per_restart
)

xi_best <- best$design
ld_best <- best$logdet

# Parameter dimension p (columns of X)
p <- ncol(build_X(xi_best))

# Relative D-efficiency of baseline compared to searched design
# (<1 means baseline is worse than searched design)
eff_baseline_vs_best <- relative_D_efficiency(
  xi_A = xi_baseline,
  xi_B = xi_best,
  cov_model = cov_model, rho = rho, phi = phi, sigma2 = sigma2,
  p = p
)

# Relative D-efficiency of searched design vs baseline
eff_best_vs_baseline <- relative_D_efficiency(
  xi_A = xi_best,
  xi_B = xi_baseline,
  cov_model = cov_model, rho = rho, phi = phi, sigma2 = sigma2,
  p = p
)

## ----------------------------------------------------------
## 8) Summarize results
## ----------------------------------------------------------

cat("=== D-optimal exact design search (correlated errors) ===\n")
cat("Covariance model:", cov_model, "\n")
if (cov_model == "ar1") cat("rho =", rho, "\n")
if (cov_model == "exp") cat("phi =", phi, "\n")
cat("n_obs =", n_obs, "\n")
cat("candidate points =", paste(range(candidate_points), collapse = " to "), "\n\n")

cat("Baseline design (unsorted):", paste(xi_baseline, collapse = ", "), "\n")
cat("Baseline design (sorted):  ", paste(sort(xi_baseline), collapse = ", "), "\n")
cat("log det M (baseline):", round(ld_baseline, 4), "\n\n")

cat("Best design found (restart", best$restart, ", unsorted):",
    paste(xi_best, collapse = ", "), "\n")
cat("Best design found (sorted): ",
    paste(sort(xi_best), collapse = ", "), "\n")
cat("log det M (best):", round(ld_best, 4), "\n\n")

cat("Relative D-efficiency (baseline vs best):", round(eff_baseline_vs_best, 4), "\n")
cat("Relative D-efficiency (best vs baseline):", round(eff_best_vs_baseline, 4), "\n")

# Create a clean table of the final design with replication counts
final_design_table <- as.data.frame(table(xi_best))
names(final_design_table) <- c("design_point", "replicates")
final_design_table$design_point <- as.numeric(as.character(final_design_table$design_point))
final_design_table <- final_design_table[order(final_design_table$design_point), ]

cat("\nFinal exact design table (for implementation):\n")
print(final_design_table, row.names = FALSE)


## ----------------------------------------------------------
## 9) Optional diagnostics / plots
## ----------------------------------------------------------

# Criterion trajectory for best restart
plot(best$history$iter, best$history$logdet, type = "b",
     xlab = "Iteration", ylab = "Best log det M so far",
     main = "Exchange search progress")

# Visualize selected design points (with replication)
plot(candidate_points, rep(0, length(candidate_points)),
     pch = 3, col = "gray",
     xlab = "Candidate point", ylab = "",
     yaxt = "n", main = "Selected design points (replication shown by stack)")
abline(h = 0, col = "gray80")

# Stack repeated points vertically for visibility
tab <- table(xi_best)
for (xv in as.numeric(names(tab))) {
  reps <- as.integer(tab[as.character(xv)])
  points(rep(xv, reps), seq_len(reps), pch = 19)
}