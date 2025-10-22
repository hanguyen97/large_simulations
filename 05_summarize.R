rm(list = ls())

setwd("~/ATTglasso/")

library(data.table)

load("out/DGP.RData")
dgp.df <- data.table(dgp.df)

R <- 50

perf <- data.table()

for (i in c(1:nrow(dgp.df))) {
  graph.type <- dgp.df[i, ]$graph.type
  graph.param <- dgp.df[i, ]$graph.param
  p <- dgp.df[i, ]$p
  n <- dgp.df[i, ]$n
  same.diag <- dgp.df[i, ]$same.diag
  
  print(paste(i, graph.type, graph.param, p, n, same.diag, R))
  data.name <- paste0(graph.type, graph.param, 
                      "_p", p, "_n", n, "_diag", same.diag, "_R", R)

  load(paste0("out/", data.name, "_compare.arr.RData"))
  compare.arr$dgp <- data.name
  perf <- rbind(perf, compare.arr, fill=TRUE)
  
}

# Output all metrics
write.csv(perf, "out/all_perf.csv")

# Summarize 
metrics <- c("auroc", "rmse", "rmse.off", "kl.div", "sen", "spe", "pre", "recall", "fdr", "acc", "f1", "run.time", "converged")
sum_tbl <- perf[, lapply(.SD, function(x) paste0(round( mean(x),2), " (", round(sd( x ),2), ")")), 
          by=c("dgp", "method"), .SDcols = metrics]
write.csv(sum_tbl, "out/sum_perf.csv")

