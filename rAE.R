


#### Script to investigate correlation between A and E ####

library(MASS)
n <- 1000 # pairs of twins in each zygosity

# input values for a, c, e, and cor(A,C)
a <- .3 # A
c <- .2 # C
e <- .5 # E
rAE <- .3 # correlation between A and E

# vectors to hold values from the simulation
Av <- vector()
Cv <- vector()
Ev <- vector()
rmzv <- vector()
rdzv <- vector()

# loop that simulates samples of 2000 twin pairs 1000 times
for(i in 1:1000){
  ## MZ
  #matrix with correlations between variance components for twin 1 and twin 2
  sigma <- matrix(c(1  ,0  ,rAE,1  ,0  ,rAE,
                    0  ,1  ,0  ,0  ,1  ,0  ,
                    rAE,0  ,1  ,rAE,0  ,0  ,
                    1  ,0  ,rAE,1  ,0  ,rAE,
                    0  ,1  ,0  ,0  ,1  ,0  ,
                    rAE,0  ,0  ,rAE,0  ,1  ), 6, 6, byrow = T)
  
  dtmz <- mvrnorm(n, mu = rep(0,6), Sigma = sigma)
  
  dtmz <- data.frame(dtmz)
  names(dtmz) <- c("A1", "C1", "E1", "A2", "C2", "E2")
  
  mz1 <- sqrt(a)*dtmz$A1 + sqrt(c)*dtmz$C1 + sqrt(e)*dtmz$E1
  mz2 <- sqrt(a)*dtmz$A2 + sqrt(c)*dtmz$C2 + sqrt(e)*dtmz$E2
  
  ## DZ
  #matrix with correlations between variance components for twin 1 and twin 2
  sigma <- matrix(c(1    ,0  ,rAE  ,.5   ,0  ,rAE/2,
                    0    ,1  ,0    ,0    ,1  ,0    ,
                    rAE  ,0  ,1    ,rAE/2,0  ,0    ,
                    .5   ,0  ,rAE/2,1    ,0  ,rAE  ,
                    0    ,1  ,0    ,0    ,1  ,0    ,
                    rAE/2,0  ,0    ,rAE  ,0  ,1    ), 6, 6, byrow = T)
  
  dtdz <- mvrnorm(n, mu = rep(0,6), Sigma = sigma)
  
  dtdz <- data.frame(dtdz)
  names(dtdz) <- c("A1", "C1", "E1", "A2", "C2", "E2")
  
  dz1 <- sqrt(a)*dtdz$A1 + sqrt(c)*dtdz$C1 + sqrt(e)*dtdz$E1
  dz2 <- sqrt(a)*dtdz$A2 + sqrt(c)*dtdz$C2 + sqrt(e)*dtdz$E2
  
  (rmz <- cor(mz1, mz2))
  (rdz <- cor(dz1, dz2))
  
  (A <- 2*(rmz - rdz))
  (C <- rmz - A)
  (E <- 1 - rmz)
  Av[i] <- A
  Cv[i] <- C
  Ev[i] <- E
  rmzv[i] <- rmz
  rdzv[i] <- rdz
}

# mean estimates for A, C, E, rmz, and rdz from the simulation
mean(Av)
mean(Cv)
mean(Ev)
mean(rmzv)
mean(rdzv)


### theoretical prediction

# additional variance
ae2r = 2*sqrt(a)*sqrt(e) * rAE

# Total variance
V = a + c + e + ae2r
# %
c(a, c, e, ae2r) / V

# Predicted rMZ and rDZ
rmz = (a + c + ae2r) / V
rdz = (a/2 + c + ae2r/2) / V

# predicted values for a, c, and e, calculated in two ways each
a/V + ae2r/V
c/V 
e/V 
2*(rmz - rdz) 
rmz - 2*(rmz - rdz) 
1 - rmz 
