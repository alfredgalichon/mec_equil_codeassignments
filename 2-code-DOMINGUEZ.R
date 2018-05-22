rm(list=ls())
library(tictoc)
library(doParallel)

############### 1. Parallelizing the Coordinate Updates Algorithm - Jacobi Version

#Set up grids

coordinates = seq(0.1,1,0.1)

X <- merge(coordinates, coordinates)
Y <- merge(coordinates, coordinates)
Z <- merge(coordinates, coordinates)

m_y  <- matrix(Y$x*Y$y)
m_yM <- matrix(m_y,100,100)

#Distance and utility function

d <- matrix(0,100,100)

for (i in c(1:100)){
  for (j in c(1:100)){
    d[i,j] <- sqrt( (X$x[i] - Z$x[j])^2 + (X$y[i] - Z$y[j])^2 )}
}

#u <- -10*(d >= 0.1) #for some reason this assigns =0 to some elements of d=0.1, so I'll use a value just below 0.1
u <- -10*(d >= 0.0999999999)
c <- d*d

#Demand and supply

D <- function(p){
  
  p <- matrix(t(p),100,100)
  num_xz <- exp(u - p)
  den_x  <- matrix(matrix(1 + rowSums(num_xz)),100,100)
  quot_x <- num_xz/den_x
  return( matrix(colSums(quot_x)) )
}

S <- function(p){
  
  p <- matrix(t(p),100,100)
  num_yz <- exp(p - c)
  den_y  <- matrix(matrix(1 + rowSums(num_yz)),100,100)
  quot_y <- m_yM*num_yz/den_y
  
  return( matrix(colSums(quot_y)) )
}

Ez <- function(p_z){
  
  if (i==1){
    p_v <- matrix( c( p_z, p_nz[1:99] ) )
  } else {
    if(i==100){
      p_v <- matrix( c( p_nz[1:99], p_z ) )
    } else {
      p_v <- matrix( c( p_nz[1:(i-1)], p_z, p_nz[i:99] ) )
    }
  }
  
  return( ( S(p_v)[i] - D(p_v)[i] )^2 )
}

##### Non - parallelized version

p_0 <- matrix( rep(0,100) )
p <- p_0
tol  <- 10^-5
precision <- 1
iter <- 0

tic()
while (precision > tol){
  
  precision <- max( abs(S(p_0)-D(p_0))/(S(p_0)+D(p_0)) )
  print(precision)
  
  
  for(i in c(1:100)){
    
    if (i==1){
      p_nz <- matrix( c( p_0[2:100] ) )
    } else {
      if(i==100){
        p_nz <- matrix( c( p_0[1:99] ) )
      } else {
        p_nz <- matrix( c( p_0[1:(i-1)], p_0[(i+1):100] ) )
      }
    }
    
    res = optim(p_0[i], Ez,  method="BFGS")
    p[i] = res$par
  } 
  
  p_0 <- p
  
  iter <- iter + 1
}
toc()

iter_1 = iter
precision_1 = precision
Dstar_1 = D(p)
Sstar_1 = S(p)
I_1 = sum(D(p))
Ip_1 =  Dstar_1*p

cbind("Precision is", precision_1)
cbind("Total number of rides is", I_1)
cbind("Average price of a ride is", sum(Ip_1))
cbind("Time taken on this machine, in number of iterations, is", iter_1)


##### Parallelized version

p_0 <- matrix( rep(0,100) )
p <- p_0
tol  <- 10^-5
precision <- 1
iter <- 0

n_cores <- 6

tic()
while (precision > tol){
  
  precision <- max( abs(S(p_0)-D(p_0))/(S(p_0)+D(p_0)) )
  print(precision)
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  pp <- foreach(i=1:100, .combine=c) %dopar% {
    
    if (i==1){
      p_nz <- matrix( c( p_0[2:100] ) )
    } else {
      if(i==100){
        p_nz <- matrix( c( p_0[1:99] ) )
      } else {
        p_nz <- matrix( c( p_0[1:(i-1)], p_0[(i+1):100] ) )
      }
    }
    res = optim(p_0[i], Ez,  method="BFGS")
    c(res$par)
  } 

  stopCluster(cl)
  
  iter <- iter + 1
  
  p_0 <- pp
}
toc()

iter_2 = iter
precision_2 = precision
Dstar_2 = D(p_0)
Sstar_2 = S(p_0)
I_2 = sum(D(p_0))
Ip_2 =  Dstar_2*p_0

cbind("Precision is", precision_2)
cbind("Total number of rides is", I_2)
cbind("Average price of a ride is", sum(Ip_2))
cbind("Time taken on this machine, in number of iterations, is", iter_2)



############### 2. Gale-Shapley Algorithm
set.seed(0)

alpha <- -(d >= 0.5)
gamma <- -d
oo <- -2

populateP <-function(Parg){
  
  M <- matrix(0,100,100)
  
  for(i in c(1:100)){
    M[i,Parg[i]]=1
  }
  return(M)
}

populateA <- function(M){
  
  A <- matrix(1,100,100)
  
  for(i in c(1:100)){
    for(j in c(1:100)){
      if (M[i,j]==1){
        A[i,]=0
        A[,j]=0
      } else {}
    }}
  return(A)}



A_0 <-  matrix( rep(1, 100,100), 100, 100 )
P_0 <-  matrix( rep(0, 100,100), 100, 100 )
M_0 <-  matrix( rep(0, 100,100), 100, 100 )

gap <-1
iter <-0

tic()
while (gap > 0){
  
  Pval <- alpha*(A_0==1) + (A_0==0)*oo
  Parg <- matrix( max.col(Pval) )
  P <- populateP(Parg)*(M_0==0)*(A_0==1)
  P_ <- t(P)
  
  Eval <- gamma*(P_==1) + (P_==0)*oo
  Earg <- matrix( max.col(Eval) )
  E <- populateP(Earg)*(P_==1)*(t(M_0)==0)
  
  M <- M_0 + P*t(E)
  A <- A_0 - P*t(E)
  A <- populateA(M)
  
  gap <- max(abs(A - A_0))
 
  iter <- iter + 1
  
  A_0 <- A
  P_0 <- P
  M_0 <- M
  
}
toc()

cbind("Passenger welfare is",  c(sum(M_0*alpha + (M_0==0)*oo)) )
cbind("Driver welfare is",  c(sum(M_0*gamma + (M_0==0)*oo)) )



############### 3. Gale-Shapley Algorithm with tie breaking
set.seed(0)

tiebreakersalpha <- 0.01*matrix(runif(1000),100,100)
tiebreakersgamma <- 0.01*matrix(runif(1000),100,100)
alpha <- -(d >= 0.5) + tiebreakersalpha
gamma <- -d + tiebreakersgamma

A_0 <-  matrix( rep(1, 100,100), 100, 100 )
P_0 <-  matrix( rep(0, 100,100), 100, 100 )
M_0 <-  matrix( rep(0, 100,100), 100, 100 )

gap <-1
iter <-0

tic()
while (gap > 0){
  
  Pval <- alpha*(A_0==1) + (A_0==0)*oo
  Parg <- matrix( max.col(Pval) )
  P <- populateP(Parg)*(M_0==0)*(A_0==1)
  P_ <- t(P)
  
  Eval <- gamma*(P_==1) + (P_==0)*oo
  Earg <- matrix( max.col(Eval) )
  E <- populateP(Earg)*(P_==1)*(t(M_0)==0)
  
  M <- M_0 + P*t(E)
  A <- A_0 - P*t(E)
  A <- populateA(M)
  
  gap <- max(abs(A - A_0))
  
  iter <- iter + 1
  
  A_0 <- A
  P_0 <- P
  M_0 <- M
  
}
toc()

cbind("Passenger welfare is",  c(sum(M_0*alpha + (M_0==0)*oo)) )
cbind("Driver welfare is",  c(sum(M_0*gamma + (M_0==0)*oo)) )
