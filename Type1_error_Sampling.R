library(MASS)
library(sfsmisc)
library(parallel)

datagen_v <- function(n, mean.gA = 0, var.gA = 1, mean.gB = 0, 
                      var.gB = 1, var.gC = 1, var.gD = 1,
                      perc.out, val.out = 10) {
  x <- factor(rep(c("A","B", "C", "D"), each = n))  
  
  y.b <- rnorm(n , mean.gB, var.gB) 
  y.c <- rnorm(n , 0, var.gC) 
  y.d <- rnorm(n , 0, var.gD) 
  
  y.a <- rnorm(n, mean.gA, var.gA) 
  if(perc.out != 0) {
    y.a[1 : ceiling(perc.out * n) ] <- rnorm(ceiling(perc.out * n), val.out, sd = 10) }
  
  y <- c(y.a, y.b, y.c, y.d)
  
  dat = data.frame(x, y)  # dataframe of data in use for simulation
  
  ## p- value of F-test from lm(), rlm-m, and rlm-mm
  lm <- summary(lm(y ~ x, data = dat))
  lm.fp <- pf(lm$fstatistic[1], lm$fstatistic[2], lm$fstatistic[3], lower.tail = F)  # LSE
  
  
  rlm.m <- rlm(y ~ x, data = dat, psi = psi.huber, init = "lts", method = "M", contrasts = list(x = "contr.sum")) # M
  rlm.fp <- f.robftest(rlm.m)$p.value
  
  rlm.mm <- rlm(y ~ x, data = dat, method = "MM", contrasts = list(x = "contr.sum")) #MM
  rlmm.fp <- f.robftest(rlm.mm)$p.value
  
  
  c(lm.fp, rlm.fp, rlmm.fp) #return vector of p-value from 3 regression methods
}

datagen_cv <- Vectorize(datagen_v, vectorize.args = c("n", "perc.out"))   


cl <- makeCluster(detectCores()-1) # set num of cores in use
clusterEvalQ(cl, library(MASS)) # pass library 
clusterEvalQ(cl, library(sfsmisc))
clusterExport(cl,c("datagen_cv"))
clusterSetRNGStream(cl) # random seeds simulation

Type1 <- parSapply(cl, 1:100, function(i,...) { 
  set.seed(NULL)    
  P <- replicate(500, {
    set.seed(NULL)
    dat = datagen_cv(n = c(20, 50, 100), mean.gA = 0, var.gA = 2, 
                     var.gB = 0.5, var.gC = 1.5, var.gD = 1,
                     perc.out = 0, val.out = 10)
    as.vector(dat)
  })   # 500 runs each of p-values from 3 different group sample size 20, 50, 100
  rowMeans(unlist(P) <= 0.05) # type1 error 
})
stopCluster(cl)