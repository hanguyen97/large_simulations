# Load libraries
library(large)
library(stringr)
library(corrplot)
library(CVglasso)
library(glasso)
library(huge)
library(parallel)
library(reshape2)
library(ggplot2)
library(PCGLASSO)

#-------------------------------------------#
#---------------- 1. Paths  ----------------#
#-------------------------------------------#

setwd("~/Documents/Cornell/Projects/FC_SC/Code/Data_analysis/fs86_hpf_tsclean_993subj/")
in_path <- "~/Documents/Cornell/Projects/FC_SC/Code/Data_analysis/fs86_hpf_tsclean_993subj/"

out_path <- "~/Documents/Cornell/Projects/FC_SC/Code/Results/Plots"
r <- 86
files <- list.files(path = in_path, pattern = "\\.RData$")



#-------------------------------------------#
#------------ 3. Estimate FC ---------------#
#-------------------------------------------#
set.seed(1)
(sj <- sample(1:length(files), 1))
regions <- 1:86

# Main loop body wrapped in a function
analyze_subject <- function(i) {
  
  load(files[i])
  fname <- basename(files[i])
  id <- sub("_TS.*", "", fname)
    # 
  idx <- 1:1200
  
  X1 <- timeseries[idx, regions]
  X2 <- timeseries[idx + 1200, regions]
  X3 <- timeseries[idx + 2400, regions]
  X4 <- timeseries[idx + 3600, regions]
  pen <- TRUE
  
  # ATTglasso
  Theta.hat.att.glasso1 <- fit_large(X=X1, alpha=0.02, thr=0.05, maxit=20, penalize_diag=pen, verbose=TRUE)
  Theta.hat.att.glasso2 <- fit_large(X=X2, alpha=0.02, thr=0.05, maxit=20, penalize_diag=pen, verbose=TRUE)
  Theta.hat.att.glasso3 <- fit_large(X=X3, alpha=0.02, thr=0.05, maxit=20, penalize_diag=pen, verbose=TRUE)
  Theta.hat.att.glasso4 <- fit_large(X=X4, alpha=0.02, thr=0.05, maxit=20, penalize_diag=pen, verbose=TRUE)
  out.att.glasso <- (Theta.hat.att.glasso1$Theta + Theta.hat.att.glasso2$Theta + 
                       Theta.hat.att.glasso3$Theta + Theta.hat.att.glasso4$Theta) / 4
  
  # huge EBIC
  out.huge <- huge(X1, method = "glasso", verbose=FALSE)
  out.ebic <- huge.select(out.huge, criterion = "ebic", verbose=FALSE)
  Theta.hat.ebic1 <- glasso(s=cov(X1), rho=out.ebic$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X2, method = "glasso", verbose=FALSE)
  out.ebic <- huge.select(out.huge, criterion = "ebic", verbose=FALSE)
  Theta.hat.ebic2 <- glasso(s=cov(X2), rho=out.ebic$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X3, method = "glasso", verbose=FALSE)
  out.ebic <- huge.select(out.huge, criterion = "ebic", verbose=FALSE)
  Theta.hat.ebic3 <- glasso(s=cov(X3), rho=out.ebic$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X4, method = "glasso", verbose=FALSE)
  out.ebic <- huge.select(out.huge, criterion = "ebic", verbose=FALSE)
  Theta.hat.ebic4 <- glasso(s=cov(X4), rho=out.ebic$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge.ebic <- (Theta.hat.ebic1+ Theta.hat.ebic2 + 
                      Theta.hat.ebic3 + Theta.hat.ebic4) / 4
  
  
  # huge RIC
  out.huge <- huge(X1, method = "glasso", verbose=FALSE)
  out.ric <- huge.select(out.huge, criterion = "ric", verbose=FALSE)
  Theta.hat.ric1 <- glasso(s=cov(X1), rho=out.ric$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X2, method = "glasso", verbose=FALSE)
  out.ric <- huge.select(out.huge, criterion = "ric", verbose=FALSE)
  Theta.hat.ric2 <- glasso(s=cov(X2), rho=out.ric$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X3, method = "glasso", verbose=FALSE)
  out.ric <- huge.select(out.huge, criterion = "ric", verbose=FALSE)
  Theta.hat.ric3 <- glasso(s=cov(X3), rho=out.ric$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X4, method = "glasso", verbose=FALSE)
  out.ric <- huge.select(out.huge, criterion = "ric", verbose=FALSE)
  Theta.hat.ric4 <- glasso(s=cov(X4), rho=out.ric$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge.ric <- (Theta.hat.ric1+ Theta.hat.ric2 + 
                     Theta.hat.ric3 + Theta.hat.ric4) / 4
  
  
  # huge STARS
  out.huge <- huge(X1, method = "glasso", verbose=FALSE)
  out.stars <- huge.select(out.huge, criterion = "stars", verbose=FALSE)
  Theta.hat.stars1 <- glasso(s=cov(X1), rho=out.stars$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X2, method = "glasso", verbose=FALSE)
  out.stars <- huge.select(out.huge, criterion = "stars", verbose=FALSE)
  Theta.hat.stars2 <- glasso(s=cov(X2), rho=out.stars$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X3, method = "glasso", verbose=FALSE)
  out.stars <- huge.select(out.huge, criterion = "stars", verbose=FALSE)
  Theta.hat.stars3 <- glasso(s=cov(X3), rho=out.stars$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge <- huge(X4, method = "glasso", verbose=FALSE)
  out.stars <- huge.select(out.huge, criterion = "stars", verbose=FALSE)
  Theta.hat.stars4 <- glasso(s=cov(X4), rho=out.stars$opt.lambda, penalize.diagonal=pen)$wi
  
  out.huge.stars <- (Theta.hat.stars1+ Theta.hat.stars2 + 
                       Theta.hat.stars3 + Theta.hat.stars4) / 4
  
  
  # CVglasso
  out.cv.loglik <- CVglasso(X=X1, S=cov(X1), crit.cv="loglik", K=5, diagonal=pen, trace="none")
  Theta.hat.cv.loglik1 <- glasso(s=cov(X1), rho=out.cv.loglik$Tuning[2], penalize.diagonal=pen)$wi
  
  out.cv.loglik <- CVglasso(X=X2, S=cov(X2), crit.cv="loglik", K=5, diagonal=pen, trace="none")
  Theta.hat.cv.loglik2 <- glasso(s=cov(X2), rho=out.cv.loglik$Tuning[2], penalize.diagonal=pen)$wi
  
  out.cv.loglik <- CVglasso(X=X3, S=cov(X3), crit.cv="loglik", K=5, diagonal=pen, trace="none")
  Theta.hat.cv.loglik3 <- glasso(s=cov(X3), rho=out.cv.loglik$Tuning[2], penalize.diagonal=pen)$wi
  
  out.cv.loglik <- CVglasso(X=X4, S=cov(X4), crit.cv="loglik", K=5, diagonal=pen, trace="none")
  Theta.hat.cv.loglik4 <- glasso(s=cov(X4), rho=out.cv.loglik$Tuning[2], penalize.diagonal=pen)$wi
  out.cv.loglik <- (Theta.hat.cv.loglik1+ Theta.hat.cv.loglik2 + 
                      Theta.hat.cv.loglik3 + Theta.hat.cv.loglik4) / 4
  
  
  # huge TIGER
  out.huge.tiger <- huge.tiger(X1, lambda=sqrt(log(dim(X1)[2])/dim(X1)[1]), 
                               verbose=FALSE)
  Theta.hat.tiger1 <- out.huge.tiger$icov[[1]]
  
  out.huge.tiger <- huge.tiger(X2, lambda=sqrt(log(dim(X2)[2])/dim(X2)[1]), 
                               verbose=FALSE)
  Theta.hat.tiger2 <- out.huge.tiger$icov[[1]]
  
  out.huge.tiger <- huge.tiger(X3, lambda=sqrt(log(dim(X3)[2])/dim(X3)[1]), 
                               verbose=FALSE)
  Theta.hat.tiger3 <- out.huge.tiger$icov[[1]]
  
  out.huge.tiger <- huge.tiger(X4, lambda=sqrt(log(dim(X4)[2])/dim(X4)[1]), 
                               verbose=FALSE)
  Theta.hat.tiger4 <- out.huge.tiger$icov[[1]]
  
  out.huge.tiger <- (Theta.hat.tiger1+ Theta.hat.tiger2 + 
                       Theta.hat.tiger3 + Theta.hat.tiger4) / 4
  
  compute_bic <- function(Theta_hat, S, n) {
    num_edges <- sum(Theta_hat[lower.tri(Theta_hat)] != 0)
    log_likelihood <- (n/2) * (determinant(Theta_hat, logarithm=TRUE)$modulus -
                                 sum(diag(S %*% Theta_hat)))
    bic <- log(n) * num_edges - 2 * log_likelihood
    return(as.numeric(bic))
  }
  
  best_rhos <- numeric(4)
  for (rho_i in 1:4) {
    lambda.grid <- seq(from=0, to=1, by = 0.1); l <- length(lambda.grid)
    bic.arr = numeric(l)
    
    for(i in 1:l){
      # Estimate
      lambda <- lambda.grid[i]
      Theta.hat.pcglasso <- pcglasso(S = cov(get(paste0("X", rho_i))), rho = lambda)
      if(any(is.na(Theta.hat.pcglasso))){
        bic.arr[i] = NA
      } else {
        bic.arr[i]<- compute_bic(Theta_hat=Theta.hat.pcglasso, S=cov(get(paste0("X", rho_i))), 
                                 n=nrow(get(paste0("X", rho_i))))
      } 
    }
    best_rhos[rho_i] <- lambda.grid[which.min(bic.arr)]
  }
  
  Theta.hat.pc1 <- pcglasso(S = cov(X1), rho = best_rhos[1])
  Theta.hat.pc2 <- pcglasso(S = cov(X2), rho = best_rhos[2])
  Theta.hat.pc3 <- pcglasso(S = cov(X3), rho = best_rhos[3])
  Theta.hat.pc4 <- pcglasso(S = cov(X4), rho = best_rhos[4])
  
  out.huge.pc <- (Theta.hat.pc1+ Theta.hat.pc2 + 
                       Theta.hat.pc3 + Theta.hat.pc4) / 4
  
  return(list(ATT = out.att.glasso, 
              EBIC = out.huge.ebic,
              CV = out.cv.loglik, 
              STARS = out.huge.stars,
              RIC = out.huge.ric, 
              TIGER = out.huge.tiger,
              PC = out.huge.pc))
}

# Run in parallel across 5 subjects (adjust `1:5` as needed)
# results <- mclapply(1, analyze_subject, mc.cores = detectCores() - 1)
results <- analyze_subject(sj)

#-------------------------------------------#
#--------- 4. Plot FC networks -------------#
#-------------------------------------------#


prec2corr <- function(Theta) {
  D <- 1 / sqrt(diag(Theta))
  P <- diag(D) %*% Theta %*% diag(D)  # correlation matrix
  return(P)
}

mat2df <- function(M, method) {
  df <- melt(M)
  colnames(df) <- c("Var1", "Var2", "value")
  df$Method <- method
  df
}

plot_corr <- function(df, method, sparsity) {
  # set diagonal values to NA
  df$value[df$Var1 == df$Var2] <- NA
  
  ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
    # # base heatmap
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red",
                        mid = "white", midpoint = 0,
                        na.value = "white",
                        limits = c(-1, 1)) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    guides(fill = "none") +
    labs(title = paste0(method, " (FC)  p = ", r, "  n = ", 1200,
                        ", sparsity = ", round(sparsity), "%")) 
}

calc_sparsity <- function(Theta.hat) {
  (1 - sum(abs(Theta.hat[upper.tri(Theta.hat)]) > 1e-6) / choose(86,2))*100
}

df <- mat2df(prec2corr(results$ATT), "Autotune")
p1 <- plot_corr(df, "LARGE", calc_sparsity(results$ATT))

df <- mat2df(prec2corr(results$EBIC), "EBIC")
p2 <- plot_corr(df, "EBIC", calc_sparsity(results$EBIC))

df <- mat2df(prec2corr(results$RIC), "RIC")
p3 <- plot_corr(df, "RIC", calc_sparsity(results$RIC))

df <- mat2df(prec2corr(results$STARS), "StARS")
p4 <- plot_corr(df, "StARS", calc_sparsity(results$STARS))

df <- mat2df(prec2corr(results$CV), "CV")
p5 <- plot_corr(df, "CV", calc_sparsity(results$CV))

df <- mat2df(prec2corr(results$TIGER), "TIGER")
p6 <- plot_corr(df, "TIGER", calc_sparsity(results$TIGER))

df <- mat2df(prec2corr(results$PC), "PCGLASSO")
p7 <- plot_corr(df, "PCGLASSO", calc_sparsity(results$PC))

library(gridExtra)
grid.arrange(p1, p2, # p3,
             p4, p5, p6,
             p7, ncol = 3)

