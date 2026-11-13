rm(list = ls())
library(MASS)
library(corrplot)
library(pROC)
library(large)
library(glasso)
library(CVglasso)

# ------------- DGP (fixed Theta, new X each replication) -------------#
# Sigma = blockdiag(dense AR(1), n_units independent 3-node path units)
# Big block: sparse Theta, but non-adjacent pairs have high marginal corr (decoys)
# Small block: in each unit, edge 1-3 has zero marginal corr but strong partial corr (hidden)
R <- 20
p_big <- 41; n_units <- 3; rho <- 0.6          # rho < 1/sqrt(2)
p_small <- 3 * n_units; p <- p_big + p_small; n <- 300

unit <- matrix(c(1, rho, 0,  rho, 1, rho,  0, rho, 1), 3, 3)
Sigma <- matrix(0, p, p)
Sigma[1:p_big, 1:p_big] <- 0.8^abs(outer(1:p_big, 1:p_big, "-"))
idx <- (p_big + 1):p
Sigma[idx, idx] <- kronecker(diag(n_units), unit)

Theta <- solve(Sigma)
Theta[abs(Theta) < 1e-8] <- 0

# ------------- Pair labels -------------#
ut     <- upper.tri(Theta)
d      <- abs(row(Theta) - col(Theta))
blk    <- matrix(FALSE, p, p); blk[idx, idx] <- TRUE
true   <- Theta != 0 & ut
hidden <- true & blk & d == 2                  # zero marginal, strong partial

# ------------- Metrics -------------#
rel.diff <- function(Est, Tru) sum((Est - Tru)^2) / sum(Tru^2)

rmse.off <- function(Theta.hat, Theta) {
  diag(Theta.hat) <- 0; diag(Theta) <- 0
  rel.diff(Est = Theta.hat, Tru = Theta)
}

# AUROC on |Theta.hat|; pos = which true edges count as positives
# 'small' = small-block edges vs all non-edges (incl. big-block decoys)
auc_pairs <- function(Theta.hat, pos) {
  keep  <- (pos | !true) & ut
  score <- round(abs(Theta.hat[keep]), 5)
  as.numeric(auc(roc(pos[keep] + 0, score, direction = "<",
                     levels = c(0, 1), quiet = TRUE)))
}

evaluate <- function(Theta.hat) {
  c(AUC      = auc_pairs(Theta.hat, true),
    RMSE_off = rmse.off(Theta.hat, Theta))
}

# ------------- Fit methods -------------#
fit_all <- function(X) {
  fits <- list()
  fits$LARGE <- fit_large(X = X, alpha = 0.02, thr = 0.005, maxit = 50,
                          penalize_diag = FALSE, verbose = FALSE)$Theta
  cv <- CVglasso(X = X, S = cov(X), crit.cv = "loglik", K = 5,
                 diagonal = FALSE, trace = "none")
  fits$CVglasso <- glasso(s = cov(X), rho = cv$Tuning[2],
                          penalize.diagonal = FALSE)$wi
  fits
}

# ------------- Replications -------------#
res <- list()
for (r in 1:R) {
  set.seed(r)
  X <- mvrnorm(n, rep(0, p), Sigma)
  fits <- fit_all(X)
  if (r == 1) { X1 <- X; fits1 <- fits }       # keep rep 1 for plots
  res[[r]] <- data.frame(r = r, method = names(fits),
                         t(sapply(fits, evaluate)), row.names = NULL)
  cat("rep", r, "done\n")
}
res <- do.call(rbind, res)

# ------------- Summary: mean (SD) over replications -------------#
metrics <- c("AUC", "RMSE_off")
summ <- aggregate(res[, metrics], by = list(method = res$method),
                  FUN = function(x) sprintf("%.3f (%.3f)", mean(x), sd(x)))
cat("\nMean (SD) over", R, "replications\n")
print(summ, row.names = FALSE)

# ------------- Plots for replication 1 -------------#
par(mfrow = c(2, 2))
corrplot(cov2cor(Sigma), method = "color", tl.pos = "n", title = "Sigma", mar = c(0, 0, 1, 0))
corrplot(cov2cor(Theta), method = "color", tl.pos = "n", title = "Theta", mar = c(0, 0, 1, 0))
for (m in names(fits1))
  corrplot(cov2cor(fits1[[m]]), method = "color", tl.pos = "n",
           title = paste0(m, " (rep 1)"), mar = c(0, 0, 1, 0))