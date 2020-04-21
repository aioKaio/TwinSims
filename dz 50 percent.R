
#### An illustration of the idea that dizygotes share 50% of their segregating genes ####

ng <- 100  # number of genes
nt <- 1000 # number of twin-pairs

musicalityt1 <- vector() # vector to hold musicality scores for twin 1
musicalityt2 <- vector() # vector to hold musicality scores for twin 2

for(t in 1:nt){
  f1 <- rnorm(ng) # generate genes on first chromosome for father
  f2 <- rnorm(ng) # generate genes on second chromosome for father
  m1 <- rnorm(ng) # generate genes on first chromosome for mother
  m2 <- rnorm(ng) # generate genes on second chromosome for mother
  
  t1f <- vector() # vector to hold genes of twin 1 from father
  t1m <- vector() # vector to hold genes of twin 1 from mother
  t2f <- vector() # vector to hold genes of twin 2 from father
  t2m <- vector() # vector to hold genes of twin 2 from mother
  
  for(i in 1:ng){
    t1f[i] <- c(f1,f2)[sample(c(i,i+ng),1)] # select from genes of father for twin 1
    t1m[i] <- c(m1,m2)[sample(c(i,i+ng),1)] # select from genes of mother for twin 1
    t2f[i] <- c(f1,f2)[sample(c(i,i+ng),1)] # select from genes of father for twin 2
    t2m[i] <- c(m1,m2)[sample(c(i,i+ng),1)] # select from genes of mother for twin 2
  }
  
  musicalityt1[t] <- sum(t1f + t1m) # sum genes to get musicality-score for twin 1
  musicalityt2[t] <- sum(t2f + t2m) # sum genes to get musicality-score for twin 2
}

cor(musicalityt1, musicalityt2)  # correlation of musicality-scores from twin 1 and twin 2
plot(musicalityt1, musicalityt2) # plot of musicality-scores from twin 1 and twin 2
