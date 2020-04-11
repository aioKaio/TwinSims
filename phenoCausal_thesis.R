


#### simulations of phenotypic causality ####

rm(list = ls())
require(OpenMx)

### functions ###
# these can be minimized

# function to generate correlated variables
kamran  <- function(x, ym, ys, co){
  y <- rnorm(length(x),0,ys)*(sqrt(1-co^2)) + (x-ym)*co + ym
  return(as.numeric(y))
}

# function to generate univariate twin data
twinVar <- function(n,v,a,c,e){
  # n = nr of twin-pairs (50% mz, 50% dz)
  # v = name for variable to simulate 
  # a,c,e = variance components
  n <- n/2
  Amz1 <- rnorm(n); Cmz1 <- rnorm(n); Emz1 <- rnorm(n)
  Amz2 <- Amz1; Cmz2 <- Cmz1; Emz2 <- rnorm(n)
  
  mz1 <- sqrt(a)*Amz1 + sqrt(c)*Cmz1 + sqrt(e)*Emz1
  mz2 <- sqrt(a)*Amz2 + sqrt(c)*Cmz2 + sqrt(e)*Emz2
  
  Adz1 <- rnorm(n); Cdz1 <- rnorm(n); Edz1 <- rnorm(n)
  Adz2 <- kamran(Adz1, 0, 1, .5); Cdz2 <- Cdz1; Edz2 <- rnorm(n)
  
  dz1 <- sqrt(a)*Adz1 + sqrt(c)*Cdz1 + sqrt(e)*Edz1
  dz2 <- sqrt(a)*Adz2 + sqrt(c)*Cdz2 + sqrt(e)*Edz2
  
  dtmz <- data.frame(cbind(c(mz1, mz2), c(1:n, 1:n), rep(1, 2*n), rep(1, 2*n)))
  names(dtmz) <- c(paste(v,"t1", sep =""), "pairnum", "zygo", "z")
  dtmz <- dtmz[order(dtmz$pairnum),]
  dtmz[,paste(v,"t2", sep ="")] <- dtmz[, paste(v,"t1", sep ="")][1:nrow(dtmz) + c(1,-1)]
  dtmz$TwinID <- rep(c(1,2), n)
  
  dtdz <- data.frame(cbind(c(dz1, dz2), c((n+1):(2*n), (n+1):(2*n)), rep(2, 2*n), rep(.5, 2*n)))
  names(dtdz) <- c(paste(v,"t1", sep =""), "pairnum", "zygo", "z")
  dtdz <- dtdz[order(dtdz$pairnum),]
  dtdz[,paste(v,"t2", sep ="")] <- dtdz[, paste(v,"t1", sep ="")][1:nrow(dtdz) + c(1,-1)]
  dtdz$TwinID <- rep(c(1,2), n)
  
  dt <- rbind(dtmz,dtdz)[,c(paste(v,"t1", sep =""), paste(v,"t2", sep =""), "pairnum", "zygo", "z", "TwinID")]
  
  return(dt)
}

# univariate heritability analysis
uniACE  <- function(dt, var){
  
  selVars <- c(paste(var,"t1", sep =""), paste(var,"t2", sep =""))
  mzData <- subset(dt, zygo==1, selVars)
  dzData <- subset(dt, zygo==2, selVars)
  
  pars <- list(mxMatrix( type="Full", 1, 1, free=TRUE, values=.5, label="a11", name="a" ),
               mxMatrix( type="Full", 1, 1, free=TRUE, values=.5, label="c11", name="c" ),
               mxMatrix( type="Full", 1, 1, free=TRUE, values=.5, label="e11", name="e" ),
               mxAlgebra( expression=a %*% t(a), name="A" ),
               mxAlgebra( expression=c %*% t(c), name="C" ),
               mxAlgebra( expression=e %*% t(e), name="E" ),
               mxAlgebra( expression=A+C+E, name="V" ))
  
  twinACEModel  <- mxModel( "ACE", pars,
                            mxModel( pars,
                                     mxMatrix( type="Full", 1, 2, free=TRUE, mean(colMeans(cbind(mzData, dzData),na.rm=TRUE)),
                                               label="mean", name="expMean" ),
                                     mxAlgebra( expression=rbind( cbind(V, A+C),
                                                                  cbind(A+C, V)), name="expCovMZ" ),
                                     mxData( observed=mzData, type="raw" ),
                                     mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars ),
                                     mxFitFunctionML(),
                                     name="MZ" ),
                            mxModel( pars,
                                     mxMatrix( type="Full", 1, 2, free=TRUE, mean(colMeans(cbind(mzData, dzData),na.rm=TRUE)),
                                               label="mean", name="expMean" ),
                                     mxAlgebra( expression=rbind( cbind(V, 0.5%x%A+C),
                                                                  cbind(0.5%x%A+C , V)), name="expCovDZ" ),
                                     mxData( observed=dzData, type="raw" ),
                                     mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars ),
                                     mxFitFunctionML(),
                                     name="DZ" ),
                            mxFitFunctionMultigroup(c("MZ", "DZ") ) )
  
  # Run Model
  twinACEFit   <- mxRun(twinACEModel, intervals=T)
  
  estVA     <- mxEval(a*a, twinACEFit)
  estVC     <- mxEval(c*c, twinACEFit)
  estVE     <- mxEval(e*e, twinACEFit)
  estVP     <- (estVA+estVC+estVE)
  estPropVA <- estVA/estVP
  estPropVC <- estVC/estVP
  estPropVE <- estVE/estVP
  
  return(cbind(estPropVA,estPropVC,estPropVE))
}

### simulation ###

np <- 1000 # nr of twin-pairs (50% mz, 50% dz)

## varians i X
errVx <- .4                    # variance in X from measurement error
aceVx <- .6*(1-errVx)          # variance in X from unique A, C, E
cauV  <- .2*(1-errVx)          # variance in X from the causal effect of Y
AconV <- .1*(1-errVx)          # variance in X from A-confounds with Y
CconV <- .1*(1-errVx)          # variance in X from C-confounds with Y
EconV <- .0*(1-errVx)          # variance in X from E-confounds with Y
conV  <- AconV + CconV + EconV # sum of variance in X from all confounds with Y

## variance in Y
errVy <- .2           # variance in Y from measurement error
aceVy <- 1-conV-errVy # variance in Y from unique A, C, E

## simulation of values for X, Y, and confounds before causal effects, confounding, and error are added
dtX    <- twinVar(np,"X",  a = .35*aceVx, c = .35*aceVx, e = .3*aceVx)
dtY    <- twinVar(np,"Y",  a = .35*aceVy, c = .35*aceVy, e = .3*aceVy) 
dtAcon <- twinVar(np,"AZ", a = 1,         c = 0,         e = .0) 
dtCcon <- twinVar(np,"CZ", a = 0,         c = 1,         e = .0) 
dtEcon <- twinVar(np,"EZ", a = 0,         c = 0,         e = .1) 

dtY <- merge(dtY, dtAcon, by = c("pairnum", "zygo", "TwinID", "z"), sort = F)
dtY <- merge(dtY, dtCcon, by = c("pairnum", "zygo", "TwinID", "z"), sort = F)
dtY <- merge(dtY, dtEcon, by = c("pairnum", "zygo", "TwinID", "z"), sort = F)
dt  <- merge(dtX, dtY,    by = c("pairnum", "zygo", "TwinID", "z"), sort = F)

## confounding
dt$Xt1 <- dt$Xt1 + dt$AZt1*sqrt(AconV) 
dt$Xt1 <- dt$Xt1 + dt$CZt1*sqrt(CconV) 
dt$Xt1 <- dt$Xt1 + dt$EZt1*sqrt(EconV) 
dt$Xt2 <- dt$Xt2 + dt$AZt2*sqrt(AconV)
dt$Xt2 <- dt$Xt2 + dt$CZt2*sqrt(CconV)
dt$Xt2 <- dt$Xt2 + dt$EZt2*sqrt(EconV)

dt$Yt1 <- dt$Yt1 + dt$AZt1*sqrt(AconV)
dt$Yt1 <- dt$Yt1 + dt$CZt1*sqrt(CconV)
dt$Yt1 <- dt$Yt1 + dt$EZt1*sqrt(EconV)
dt$Yt2 <- dt$Yt2 + dt$AZt2*sqrt(AconV)
dt$Yt2 <- dt$Yt2 + dt$CZt2*sqrt(CconV)
dt$Yt2 <- dt$Yt2 + dt$EZt2*sqrt(EconV)

## causal effect
K <- sqrt(cauV + conV^2) - conV # value that makes var(X + Y*K) equal to 1-errVx

dt$Xt1 <- dt$Xt1 + dt$Yt1*K
dt$Xt2 <- dt$Xt2 + dt$Yt2*K
# the causal effect might end up weaker than what's specified, if there are also confounds
# the variance of X rises by the right amount, but much of this variance is from the confounds
# due to the variance sum law, and so on

## measurement error
dt$Xt1 <- dt$Xt1 + rnorm(length(dt$Xt1), 0, sqrt(errVx))
dt$Xt2 <- dt$Xt2 + rnorm(length(dt$Xt2), 0, sqrt(errVx))

dt$Yt1 <- dt$Yt1 + rnorm(length(dt$Yt1), 0, sqrt(errVy))
dt$Yt2 <- dt$Yt2 + rnorm(length(dt$Yt2), 0, sqrt(errVy))

# sanity check: these should equal ~1
var(dt$Xt1)
var(dt$Yt1)

# adjustment to how OpenMx apparently dislikes 0 as average?
dt$Xt1 <- dt$Xt1 + 5
dt$Xt2 <- dt$Xt2 + 5
dt$Yt1 <- dt$Yt1 + 5
dt$Yt2 <- dt$Yt2 + 5

### Biometric ###
Pmz <- dt[dt$TwinID == 1 & dt$z == 1,   c("Yt1", "Xt1", "Yt2", "Xt2")]
Pdz <- dt[dt$TwinID == 1 & dt$z == 0.5, c("Yt1", "Xt1", "Yt2", "Xt2")]

P = rbind(Pmz, Pdz)
colnames(P) = c("t1_y", "t1_x", "t2_y", "t2_x")
datW = data.frame(P,
                  z = rep(c(1, 0.5), times = c(nrow(Pmz), nrow(Pdz))),
                  pair = 1:(nrow(Pmz) + nrow(Pdz)))

mod_cause = mxModel("mod_cause",
                    # Parameters
                    mxMatrix("Symm", 2, 2, c(T, T, T), c(1.0, 0.5, 1.0), c("va11", "va21", "va22"), name = "Va"),
                    mxMatrix("Symm", 2, 2, c(T, T, T), c(1.0, 0.5, 1.0), c("vc11", "vc21", "vc22"), name = "Vc"),
                    mxMatrix("Symm", 2, 2, c(T, F, T), c(1.0, 0.0, 1.0), c("ve11", "ve21", "ve22"), name = "Ve"),
                    mxMatrix("Full", 2, 2, c(F, T, F, F), c(0.0, 0.5, 0.0, 0.0), c("p11", "p21", "p12", "p22"), name = "p"),
                    # Covariance
                    mxMatrix("Iden", 2, 2, name = "I"),
                    mxAlgebra(solve(I - p) %&% (Va + Vc + Ve), name = "Cw"),
                    mxAlgebra(solve(I - p) %&% (data.z*Va + Vc), name = "Cb"),
                    mxAlgebra(rbind(cbind(Cw, Cb),
                                    cbind(Cb, Cw)), name = "C"),
                    # Means
                    mxMatrix("Full", 2, 1, T, 0, c("mx", "my", "mx", "my"), name = "Intr"),
                    mxAlgebra(solve(I - p) %*% Intr, name = "M1"),
                    mxAlgebra(t(rbind(M1, M1)), name = "M"),
                    # Expectations
                    mxExpectationNormal("C", "M", colnames(datW)[1:4]),
                    mxFitFunctionML(),
                    # Data
                    mxData(datW,"raw"))

# full ACE-ß
fit_cause = mxTryHard(mod_cause)
summary(fit_cause)
output <- data.frame(cbind(round(fit_cause$output$estimate, 4), round(fit_cause$output$standardErrors, 4)))
names(output) <- c("estimate", "std.error"); output

# without the causal path
mod_common = omxSetParameters(mod_cause, "p21", F , 0, name = "mod_common")
fit_common = mxRun(mod_common)

mxCompare(fit_cause, list(fit_common))

uniACE(dt, "X")
uniACE(dt, "Y")


