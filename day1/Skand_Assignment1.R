library(xtable)

setwd("/Users/skand/Dropbox/Galichon_mec")

#=========================== Create Data ==========================================================
z1.loc <- rep((1:10)*0.1, 10)
z2.loc <- unlist(lapply((1:10)*0.1, function(x) rep(x, 10)))
My <- z1.loc*z2.loc

z.loc <- cbind(z1.loc, z2.loc)
sqdist <- function(x,y) (x[1] - y[1])^2 + (x[2] - y[2])^2

Cyz <- matrix(0 ,100, 100)
for (i in 1:100) {
  for (j in 1:100) {
    Cyz[i,j] <- sqdist(z.loc[i,], z.loc[j,])
  }
}

Uxz <- matrix(-10, 100, 100)
diag(Uxz) <- 0

#=========================== Define demand and supply functions ===================================
Dz <- function(p) {
  p.matrix <- matrix(rep(p, each = 100), 100, 100)
  D.exp <- exp(Uxz - p.matrix)
  prob.xz <- D.exp/(1 + rowSums(D.exp))
  sum.z <- colSums(prob.xz)
  return(sum.z)
}

Sz <- function(p) {
  p.matrix <- matrix(rep(p, each = 100), 100, 100)
  D.exp <- exp(p.matrix - Cyz)
  prob.yz <- D.exp/(1 + rowSums(D.exp))
  sum.y <- colSums(prob.yz*My)
  return(sum.y)
}

Ez <- function(p) {
  E <- Sz(p) - Dz(p)
  return(E)
}

#=============================== Gradient descent =================================================
time0 <- Sys.time()
p <- rep(0,100)
iter <- 0
precision <- 10
while (precision > 10e-5) {
  S <- Sz(p) 
  D <- Dz(p)
  precision <- max(abs(S - D)/(S + D))
  p <- p - 0.2*(S - D)
  iter <- iter + 1
  print(precision)
}
time1 <- Sys.time()

results_gd <- c(precision, sum(Dz(p)), sum(Dz(p)*p), iter, time1 - time0)

#======================================== Newton descent ==========================================

# Gradient function
E.gradient <- function(p) {
  p.matrix <- matrix(rep(p, each = 100), 100, 100)
  D.exp <- exp(Uxz - p.matrix)
  D.prob.xz <- D.exp/(1 + rowSums(D.exp))
  S.exp <- exp(p.matrix - Cyz)
  S.prob.yz <- S.exp/(1 + rowSums(S.exp))
  
  Gradient <- matrix(0, 100, 100)
  # Diagonal elements
  diag(Gradient) <- colSums(S.prob.yz*(1-S.prob.yz)*My) + colSums(D.prob.xz*(1-D.prob.xz))

  # Off diagonal elements
  for (i in 1:100) {
    j = i + 1
    while (j <= 100) {
      Gradient[i,j] <- - sum(S.prob.yz[,i]*S.prob.yz[,j]*My) - sum(D.prob.xz[,i]*D.prob.xz[,j])
      Gradient[j,i] <- Gradient[i,j]
      j <- j + 1
    }
  }
  

  return(Gradient)
}

# Algorithm
time0 <- Sys.time()
p <- rep(0,100)
iter <- 0
precision <- 10
while (precision > 10e-5) {
  S <- Sz(p) 
  D <- Dz(p)
  precision <- max(abs(S - D)/(S + D))
  p <- p - 0.2* solve(E.gradient(p)) %*% (S - D)
  iter <- iter + 1
  print(precision)
}
time1 <- Sys.time()

results_nd <- c(precision, sum(Dz(p)), sum(Dz(p)*p), iter, time1 - time0)

#================================== Coordinate descent ============================================

# define E function with separate argument for p
E.coordinate <- function(p_i) {
  P0[i] <- p_i
  return(Ez(P0)[i]^2)
}

# Algorithm
time0 <- Sys.time()
P0 <- rep(0, 100)
P1 <- rep(0, 100)
iter <- 0
precision <- 10
while (precision > 10e-5) {
  S <- Sz(P0) 
  D <- Dz(P0)
  precision <- max(abs(S - D)/(S + D))
  for (i in 1:100) {
  P1[i] <- optim(P0[i], E.coordinate, 
                 method = "Brent", lower = P0[i] - 100, upper = P0[i] + 100)$par
  }
  P0 <- P1
  iter <- iter + 1
  print(precision)
}
time1 <- Sys.time()

results_cd <- c(precision, sum(Dz(P0)), sum(Dz(P0)*P0), iter, time1 - time0)

#============================== Combine results ===================================================
results_all <- cbind(results_gd, results_nd, results_cd)
row.names(results_all) <- c("Precision", "Total rides", "Total price", "niter", "time taken in sec")
colnames(results_all) <- c("Gradient descent", "Newton descent", "Coordinate descent")

digitmat <- matrix(rep(c(7,3,3,0,3), 4), 5, 4)
print(xtable(results_all, digits = digitmat, align = "l|c|c|c"), 
      file = "Assignment1/Results1.tex",
      floating = FALSE,
      type = "latex", sanitize.text.function = function(x) {x})

#============================ Supply with subsidies ==============================================

# Supply function with tax
Sz_tax <- function(p) {
  p.matrix <- matrix(rep(p, each = 100), 100, 100)
  p.matrix <- p.matrix + log(1 + exp(-p.matrix))
  D.exp <- exp(p.matrix - Cyz)
  prob.yz <- D.exp/(1 + rowSums(D.exp))
  sum.y <- colSums(prob.yz*My)
  return(sum.y)
}

# Excess supply funtion with tax
Ez_tax <- function(p) {
  E <- Sz_tax(p) - Dz(p)
  return(E)
}

# define E function with separate argument for p
E.coordinate_tax <- function(p_i) {
  P0[i] <- p_i
  return(Ez_tax(P0)[i]^2)
}

# Algorithm
time0 <- Sys.time()
P0 <- rep(0, 100)
P1 <- rep(0, 100)
iter <- 0
precision <- 10
while (precision > 10e-5) {
  S <- Sz_tax(P0) 
  D <- Dz(P0)
  precision <- max(abs(S - D)/(S + D))
  for (i in 1:100) {
    P1[i] <- optim(P0[i], E.coordinate, 
                   method = "Brent", lower = P0[i] - 10, upper = P0[i] + 20)$par
  }
  P0 <- P1
  iter <- iter + 1
  print(precision)
}
time1 <- Sys.time()

results_cd <- c(precision, sum(Dz(P0)), sum(Dz(P0)*P0), iter, time1 - time0)

