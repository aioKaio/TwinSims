

#### simulation of twin data with all assumptions fulfilled ####

# function to create variables that correlate
kamran <- function(x, co){
  y <- rnorm(length(x))*(sqrt(1-co^2))+x*co
  return(as.numeric(y))
}


# Specify variance components
a <- .4
c <- .1
e <- .5

# Draw from standard normal to get A, C, and E vectors for each person
AV <- rnorm(2000)
CV <- rnorm(2000)
EV <- rnorm(2000)

# multiply by squareroot of a, c, and e and add up
data <- sqrt(a)*AV + sqrt(c)*CV + sqrt(e)*EV

Amz1 <- AV[1:1000]       # A-vector for MZ-twin 1
Cmz1 <- CV[1:1000]       # C-vector for MZ-twin 1
Emz1 <- EV[1:1000]       # E-vector for MZ-twin 1

Adz1 <- AV[1001:2000]    # A-vector for DZ-twin 1
Cdz1 <- CV[1001:2000]    # C-vector for DZ-twin 1
Edz1 <- EV[1001:2000]    # E-vector for DZ-twin 1

mz1 <- data[1:1000]      # trait-scores for MZ-twin 1
dz1 <- data[1001:2000]   # trait-scores for DZ-twin 1

Amz2 <- Amz1             # A-vector for MZ-twin 2
Cmz2 <- Cmz1             # C-vector for MZ-twin 2
Emz2 <- rnorm(1000)      # E-vector for MZ-twin 2

Adz2 <- kamran(Adz1, .5) # A-vector for DZ-twin 2
Cdz2 <- Cdz1             # C-vector for DZ-twin 2
Edz2 <- rnorm(1000)      # E-vector for DZ-twin 2

mz2 <- sqrt(a)*Amz2 + sqrt(c)*Cmz2 + sqrt(e)*Emz2 # trait-scores for MZ-twin 2
dz2 <- sqrt(a)*Adz2 + sqrt(c)*Cdz2 + sqrt(e)*Edz2 # trait-scores for DZ-twin 2

rmz <- cor(mz1, mz2)
rdz <- cor(dz1, dz2)

A <- 2*(rmz - rdz)
C <- 2*rdz - rmz
E <- 1 - rmz

A
C
E
