library(pROC)
library(plot3D)
library(gdata)
library(clime)
library(Matrix)

col_palette <- gray((100:0/100)^1.5)

rel.diff <- function(Est, Tru){
  return(sum((Est - Tru)^2)/sum(Tru^2))
}

ROC.Theta <- function(Theta.hat, Theta){
  truth <- (!!upperTriangle(round(Theta, 5), diag = FALSE)) + 0
  pred <- upperTriangle(round((Theta.hat), 5), diag = FALSE)
  roc.obj <- pROC::roc(truth, pred, direction = "<", levels = c(0, 1))
  return(roc.obj)
}

KL.div <- function(Est, Tru) {
  
  Est <- (Est + t(Est)) / 2
  Tru <- (Tru + t(Tru)) / 2
  
  Sigma.tru <- solve(Tru)
  trace_term <- sum(Sigma.tru * Est)  
  logdet_term <- determinant(Sigma.tru %*% Est, logarithm = TRUE)$modulus
  p <- ncol(Tru)
  
  # Final KL divergence
  D_kl <- 0.5 * (trace_term - as.numeric(logdet_term) - p)
  
  return(D_kl)
}

performance <- function(Theta, Theta.hat){
  truth <- as.factor((!!upperTriangle(round(Theta, 5), diag = FALSE)) + 0)
  preds <- as.factor((!!upperTriangle(round(Theta.hat, 5), diag = FALSE)) + 0)
  t <- table(truth, preds)
  if(ncol(t) == 1)
    t <- cbind(t, c(0, 0))
  a <- t[1, 1]
  b <- t[1, 2]
  c <- t[2, 1]
  d <- t[2, 2]

  pre <- ifelse(b+d, d/(b+d), NA)
  recall <- ifelse(c+d, d/(c+d), NA)
  acc <- ifelse(a+b+c+d, (a+d)/(a+b+c+d), NA)
  f1 <- ifelse(b+c+2*d, d/(d + (b+c)/2), NA)
  fdr <- ifelse(b+d, b/(b+d), NA) 
  sen <- ifelse(d+c, d/(d+c), NA)
  spe <- ifelse(a+b, a/(a+b), NA)

  return(c(pre, recall, acc, f1, fdr, 
           sen, spe, 
           tp=d, fp=b, fn=c, tn=a))
}

get_band1_graph <- function(p, rho1=0.3, same.diag = 0,
                            diag_sc = c(10, 1, 0.5)) {
  if (same.diag == 0) {
    p1 <- round(p/3)
    p2 <- round(p/3)
    p3 <- p - (p1+p2)
    p_blocks <- c(p1, p2, p3)
    
    Theta_blocks <- list()
    for (i in 1:3) {
      p <- p_blocks[i]
      A <- matrix(data=0, nrow=p, ncol=p)
      diag(A) <- 1
      diag(A[1:(p-1), 2:p]) <- rho1
      diag(A[2:p, 1:(p-1)]) <- rho1
      Theta_blocks[[i]] <- diag_sc[i] * A
    }
    
    Theta <- as.matrix(bdiag(Theta_blocks[[1]], Theta_blocks[[2]], Theta_blocks[[3]]))
    return(Theta)
    
  }
}

get_band2_graph <- function(p, rho1=0.3, rho2=0.2, same.diag = 0,
                            diag_sc = c(10, 1, 0.5)) {
  if (p < 12) {
    stop("p must be at least 12")
  }
  if (same.diag == 0) {
    p1 <- round(p/3)
    p2 <- round(p/3)
    p3 <- p - (p1+p2)
    p_blocks <- c(p1, p2, p3)
    
    Theta_blocks <- list()
    for (i in 1:3) {
      p <- p_blocks[i]
      A <- matrix(data=0, nrow=p, ncol=p)
      diag(A) <- 1
      diag(A[1:(p-1), 2:p]) <- rho1
      diag(A[2:p, 1:(p-1)]) <- rho1
      diag(A[1:(p-2), 3:p]) <- rho2
      diag(A[3:p, 1:(p-2)]) <- rho2
      Theta_blocks[[i]] <- diag_sc[i] * A
    }
    
    Theta <- as.matrix(bdiag(Theta_blocks[[1]], Theta_blocks[[2]], Theta_blocks[[3]]))
    return(Theta)
    
  }
}

get_scalefree_graph <- function(p, same.diag = 0,
                                diag_sc = c(10, 1, 0.5)) {
  Omega_0 <- huge.generator(d=p, graph="scale-free", verbose=FALSE)$omega
  
  # Standardize to have unit diagonal
  if (same.diag > 0) {
    D <- sqrt(diag(Omega_0))
    Omega_0 <- diag(sqrt(same.diag)/D) %*% Omega_0 %*% diag(sqrt(same.diag)/D)
  } else {
    blocks <- diag_sc
    block_sizes <- sort(rep(blocks, length.out = p))  
    D <- sqrt(diag(Omega_0))
    Omega_0 <- diag(sqrt(block_sizes)/D) %*% Omega_0 %*% diag(sqrt(block_sizes)/D)
  }
  
  return(Omega_0)
}


get_block_graph <- function(p, prob=0.02, rho=0.5, same.diag = 0,
                             diag_sc = c(10, 1, 0.5)) {
  if (same.diag == 0) {
    p1 <- round(p/3)
    p2 <- round(p/3)
    p3 <- p - (p1+p2)
    p_blocks <- c(p1, p2, p3)
    
    Theta_blocks <- list()
    for (i in 1:3) {
      p <- p_blocks[i]
      B <- matrix(0, p, p)
      off_diag_indices <- which(upper.tri(B), arr.ind = TRUE)
      random_vals <- runif(nrow(off_diag_indices)) < prob
      B[off_diag_indices] <- rho * random_vals
      B <- (B + t(B))/2  # Make symmetric
      diag(B) <- 1
      Theta_blocks[[i]] <- diag_sc[i] * B
    }
    
    Theta <- as.matrix(bdiag(Theta_blocks[[1]], Theta_blocks[[2]], Theta_blocks[[3]]))
    return(Theta)
    
  }
}

get_hub_graph <- function(p, rho=0.3, same.diag = 0, 
                            diag_sc = c(10, 1, 0.5)) {
  if (same.diag == 0) {
    p1 <- round(p/3)
    p2 <- round(p/3)
    p3 <- p - (p1+p2)
    p_blocks <- c(p1, p2, p3)
    
    Theta_blocks <- list()
    for (i in 1:3) {
      p <- p_blocks[i]
      B <- matrix(0, p, p)
      G <- huge.generator(graph = "hub", d=p, g=p/5, verbose=FALSE)$omega > 1e-16
      B <- rho * G
      B <- (B + t(B))/2  # Make symmetric
      diag(B) <- 1
      Theta_blocks[[i]] <- diag_sc[i] * B
    }
    
    Theta <- as.matrix(bdiag(Theta_blocks[[1]], Theta_blocks[[2]], Theta_blocks[[3]]))
    return(Theta)
    
  }
}

get_random_graph <- function(p, rho = 0.5, same.diag = 0, prob = 0.01,
                                diag_sc = c(10, 1, 0.5)) {
  B <- matrix(0, p, p)
  off_diag_indices <- which(upper.tri(B), arr.ind = TRUE)
  random_vals <- runif(nrow(off_diag_indices)) < prob
  B[off_diag_indices] <- rho * random_vals
  B <- (B + t(B))/2  # Make symmetric
  diag(B) <- 1
  
  # Standardize to have unit diagonal
  if (same.diag==0) {
    block_sizes <- sort(rep(diag_sc, length.out = p))  
    D <- sqrt(diag(B))
    out <- diag(sqrt(block_sizes)/D) %*% B %*% diag(sqrt(block_sizes)/D)
  }
  
  alpha <- min(eigen(out)$values)
  if (alpha < 0) {
    # print(alpha)
    diag(out) <- diag(out) + (abs(alpha) + 0.005)
  }
  return(out)
}

get_perf_row <- function(Theta.hat, Theta, method, r, converged, endT) {

  if (is.null(Theta.hat) || mean(is.na(Theta.hat)) > 0.5) {
    return(data.frame(method=method, 
                      rep=r, 
                      rmse=NA,
                      rmse.off=NA,
                      rel.rmse=NA,
                      kl.div=NA,
                      auroc=NA,
                      sen=NA,
                      spe=NA,
                      pre=NA,
                      recall=NA,
                      acc=NA,
                      f1=NA,
                      fdr=NA,
                      converged=FALSE,
                      run.time=NA))
  }
  units(endT) <- "secs"

  rmse <- rel.diff(Est=Theta.hat, Tru=Theta)

  Theta.hat.off <- Theta.hat
  diag(Theta.hat.off) <- 0
  Theta.off <- Theta
  diag(Theta.off) <- 0
  rmse.off <- rel.diff(Est=Theta.hat.off, Tru=Theta.off)

  perf <- performance(Theta=Theta, Theta.hat=Theta.hat)

  ROC.Theta(Theta.hat=Theta.hat, Theta=Theta)$auc

  return(data.frame(method=method, 
             rep=r, 
             rmse=rmse,
             rmse.off=rmse.off,
             kl.div = KL.div(Est=Theta.hat, Tru=Theta),
             auroc=ROC.Theta(Theta.hat=Theta.hat, Theta=Theta)$auc,
             sen=perf[6],
             spe=perf[7],
             pre=perf[1],
             recall=perf[2],
             acc=perf[3],
             f1=perf[4],
             fdr=perf[5],
             converged=converged,
             run.time=endT))       
}

