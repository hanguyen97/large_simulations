#-------------------------------------------#
#---------------- 1. Paths  ----------------#
#-------------------------------------------#

setwd("")
in_path <- ""

r <- 86
files <- list.files(path = in_path, pattern = "\\.RData$")


#-------------------------------------------#
#------------ 2. Plot diagonals ------------#
#-------------------------------------------#
set.seed(1)
(sj <- sample(1:length(files), 1))
regions <- 1:86

load(files[sj])
fname <- basename(files[sj])  
print(fname)

id <- sub("_TS.*", "", fname) 
X = timeseries[1:1200,]

hist(diag(solve(cov(X))), 
     main = "", 
     xlab = "Sample Partial Variances", 
     col = "black", border = "white", 
     breaks=86)


