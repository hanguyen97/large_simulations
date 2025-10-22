rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
i_input <- as.numeric(args[1])

setwd("~/ATTglasso/")
library(foreach)
library(doMC)
registerDoMC(40)

options(cores = 40)

# Load data
load("out/DGP.RData")
R <- 50
for (i in i_input) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  
  print(paste(i, graph.type, graph.param, p, n, same.diag, R))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R)
  
  load(paste0("out/", data.name, "_data.RData"))
  
  # Run autotune glasso and benchmark
  compare.arr = foreach(r = 1:R, .combine = rbind) %dopar% {

    pen <- FALSE
    
    set.seed(r)
    source("library_test.R")
    
    library(huge)
    library(glasso)
    library(CVglasso)
    library(large)
    library(clime)
    library(PCGLASSO)
    
    # Get ground truth
    Theta <- Theta_list[[r]]
    X <- X_list[[r]]
        
    # Run autotune glasso 2%
    if (p==100) {
      att_thr=0.005
    } else if (p==300) {
      att_thr=0.05
    }
    startT <- Sys.time()
    out.att.glasso <- fit_large(X=X, alpha=0.02, thr=att_thr, maxit=20, penalize_diag=pen, verbose=TRUE)
    endT <- Sys.time()-startT
    Theta.hat.att.glasso <- out.att.glasso$Theta

    temp.att.glasso_a2 <- get_perf_row(Theta.hat=Theta.hat.att.glasso, 
                                       Theta=Theta, 
                                       method="large", 
                                       r=r, 
                                       converged=out.att.glasso$converged, 
                                       endT=endT) 

    # Run huge ebic
    startT <- Sys.time()
    out.huge <- huge(X, method = "glasso", verbose=FALSE)
    out.huge.ebic <- huge.select(out.huge, criterion = "ebic", verbose=FALSE)
    endT <- Sys.time()-startT

    Theta.hat.huge.ebic <- glasso(s=cov(X),rho=out.huge.ebic$opt.lambda, penalize.diagonal=pen)$wi
    temp.huge.ebic <- get_perf_row(Theta.hat=Theta.hat.huge.ebic, 
                                   Theta=Theta, 
                                   method="huge.ebic", 
                                   r=r, 
                                   converged=TRUE, 
                                   endT=endT) 
    # Run huge ric
    startT <- Sys.time()
    out.huge <- huge(X, method = "glasso", verbose=FALSE)
    out.huge.ric <- huge.select(out.huge, criterion = "ric", verbose=FALSE)

    endT <- Sys.time()-startT
    if (is.null(out.huge.ric$opt.icov)) {
      Theta.hat.huge.ric <- NULL
    } else {
      Theta.hat.huge.ric <- glasso(s=cov(X),rho=out.huge.ric$opt.lambda, penalize.diagonal=pen)$wi
    }
    temp.huge.ric <- get_perf_row(Theta.hat=Theta.hat.huge.ric, 
                                   Theta=Theta, 
                                   method="huge.ric", 
                                   r=r, 
                                   converged=TRUE, 
                                   endT=endT) 
    
    # Run huge Stars
    startT <- Sys.time()
    out.huge <- huge(X, method = "glasso", verbose=FALSE)
    out.huge.stars <- huge.select(out.huge, criterion = "stars", verbose=FALSE)

    endT <- Sys.time()-startT
    if (is.null(out.huge.stars$opt.icov)) {
      Theta.hat.huge.stars <- NULL
    } else {
      Theta.hat.huge.stars <- glasso(s=cov(X),rho=out.huge.stars$opt.lambda, penalize.diagonal=pen)$wi
    }
    temp.huge.stars <- get_perf_row(Theta.hat=Theta.hat.huge.stars, 
                                  Theta=Theta, 
                                  method="huge.stars", 
                                  r=r, 
                                  converged=TRUE, 
                                  endT=endT) 
    
    # Run CVglasso loglik - use original glasso function
    startT <- Sys.time()
    out.cv.loglik = CVglasso(X=X, S=cov(X), crit.cv="loglik", K=5, diagonal=pen, trace="none")
    endT <- Sys.time()-startT
    Theta.hat.cv.loglik <- glasso(s=cov(X),rho=out.cv.loglik$Tuning[2], penalize.diagonal=pen)$wi
    temp.cv.loglik <- get_perf_row(Theta.hat=Theta.hat.cv.loglik, 
                                    Theta=Theta, 
                                    method="cv.loglik", 
                                    r=r, 
                                    converged=TRUE, 
                                    endT=endT) 
    
    # Run huge TIGER
    startT <- Sys.time()
    out.huge.tiger <- huge.tiger(X, lambda=sqrt(log(dim(X)[2])/dim(X)[1]), 
                                 verbose=FALSE)
    endT <- Sys.time()-startT
    Theta.hat.huge.tiger <- out.huge.tiger$icov[[1]]
    temp.huge.tiger <- get_perf_row(Theta.hat=Theta.hat.huge.tiger, 
                                    Theta=Theta, 
                                    method="huge.tiger", 
                                    r=r, 
                                    converged=TRUE, 
                                    endT=endT) 

    # Run PC glasso
    compute_bic <- function(Theta_hat, S, n) {
      num_edges <- sum(Theta_hat[lower.tri(Theta_hat)] != 0)
      log_likelihood <- (n/2) * (determinant(Theta_hat, logarithm=TRUE)$modulus -
      sum(diag(S %*% Theta_hat)))
      bic <- log(n) * num_edges - 2 * log_likelihood
      return(as.numeric(bic))
    }

    startT <- Sys.time()
    lambda.grid <- seq(from=0, to=1, by = 0.1); l <- length(lambda.grid)
    bic.arr = numeric(l)
    
    for(i in 1:l){
      # Estimate
      lambda <- lambda.grid[i]
      Theta.hat.pcglasso <- pcglasso(S = cov(X), rho = lambda)
      if(is.null(Theta.hat.pcglasso)){
          bic.arr[i] = NA
        } else {
          bic.arr[i]<- compute_bic(Theta_hat=Theta.hat.pcglasso, S=cov(X), n=nrow(X))
        } 
    }
    best_lambda <- lambda.grid[which.min(bic.arr)]
    Theta.hat.pcglasso <- pcglasso(S = cov(X), rho = best_lambda)

    endT <- Sys.time()-startT
    temp.pcglasso <- get_perf_row(Theta.hat= Theta.hat.pcglasso, 
                                    Theta=Theta, 
                                    method="pcglasso", 
                                    r=r, 
                                    converged=TRUE, 
                                    endT=endT) 
    
    # # Run CLIME
    # startT <- Sys.time()
    # re.clime <- clime(X, standardize=FALSE)
    # re.cv <- cv.clime(re.clime)
    # re.clime.opt <- clime(X, re.cv$lambdaopt)
    # endT <- Sys.time()-startT
    # Theta.hat.clime <- re.clime.opt$Omegalist[[1]]
    # temp.clime<- get_perf_row(Theta.hat=Theta.hat.clime, 
    #                                 Theta=Theta, 
    #                                 method="clime", 
    #                                 r=r, 
    #                                 converged=TRUE, 
    #                                 endT=endT) 
    
    # Output
    temp <- rbind(temp.att.glasso_a2, 
                  # temp.att.glasso_a5, 
                  temp.huge.ebic, 
                  temp.huge.ric, 
                  temp.huge.stars,
                  temp.cv.loglik,
                  temp.huge.tiger,
                  temp.pcglasso
                  # ,temp.clime
                  )
    
    temp
  }
  
  save(compare.arr, file=paste0("out/", data.name, "_compare.arr.RData"))
  
}
