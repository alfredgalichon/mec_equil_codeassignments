rm(list=ls())
library(tictoc)

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

#u <- -10*(d >= 0.1)
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

############### 1. Gradient Descent Method

p_0 <- rep(0,100)
epsilon <- 0.2
tol  <- 10^-5
precision <- 1
iter <- 0

tic()
while (precision > tol){
   
  iter <- iter + 1
  
  p <- p_0 - epsilon*( S(p_0)-D(p_0) )
  
  precision <- max( abs(S(p_0)-D(p_0))/(S(p_0)+D(p_0)) )

  p_0 <- p
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


############### 2. Newton method

DE <- function(p){
  
  DE <- matrix(0,100,100)
  
  for (i in c(1:100)){
    for (j in c(1:100)){
      DE[i,j] <- sum(  -m_y*exp(p[i]-c[,i])*exp(p[j]-c[,j])/((1+   sum(exp(p[j]-c[,j])) )^2)  ) -  sum( ( exp(u[,i]-p[i])*exp(u[,j]-p[j]) )/( (1+ sum( exp(u[,j]-p[j]) )  )^2  ) )
      }
  }
  
  Ss <- S(p)
  Ds <- D(p)
  
  p <- matrix(t(p),100,100)
  num_xz <- exp(u - p)
  den_x  <- matrix(matrix(1 + rowSums(num_xz)),100,100)
  quot_x <- (num_xz/den_x)*(num_xz/den_x)
  sum_x <- colSums(quot_x)
  
  num_yz <- exp(p - c)
  den_y  <- matrix(matrix(1 + rowSums(num_yz)),100,100)
  quot_y <- m_yM*( (num_yz/den_y)*(num_yz/den_y) )
  sum_y <- colSums(quot_y)
  
  for(i in c(1:100)){
    DE[i,i] <- Ss[i] + Ds[i] - sum_y[i] - sum_x[i]
  }
        
  return(DE)
}

p_0 <- rep(0,100)
epsilon <- 0.2
tol  <- 10^-5
precision <- 1
iter <- 0

tic()
while (precision > tol){
  
  iter <- iter + 1
  
  p <- p_0-epsilon*solve(DE(p_0))%*%( S(p_0)-D(p_0) )
  
  precision <- max( abs(S(p_0)-D(p_0))/(S(p_0)+D(p_0)) )
  
  p_0 <- p
  
  print(precision)
  
}
toc()

iter_2 = iter
precision_2 = precision
Dstar_2 = D(p)
Sstar_2 = S(p)
I_2 = sum(D(p))
Ip_2 =  Dstar_2*p

cbind("Precision is", precision_2)
cbind("Total number of rides is", I_2)
cbind("Average price of a ride is", sum(Ip_2))
cbind("Time taken on this machine, in number of iterations, is", iter_2)
  

############### 3a. Coordinate Updates - Jacobi Version

p_0 <- matrix( rep(0,100) )
tol  <- 10^-5
precision <- 1
iter <- 0

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

iter_3 = iter
precision_3 = precision
Dstar_3 = D(p)
Sstar_3 = S(p)
I_3 = sum(D(p))
Ip_3 =  Dstar_3*p

cbind("Precision is", precision_3)
cbind("Total number of rides is", I_3)
cbind("Average price of a ride is", sum(Ip_3))
cbind("Time taken on this machine, in number of iterations, is", iter_3)




############### 4. Coordinate Updates - Jacobi Version - with subsidies

S <- function(p){
  
  p <- matrix(t(p),100,100)
  Fp <- p + log(1+exp(-p)) 
  num_yz <- exp(Fp - c)
  den_y  <- matrix(matrix(1 + rowSums(num_yz)),100,100)
  quot_y <- m_yM*num_yz/den_y
  
  return( matrix(colSums(quot_y)) )
}


p_0 <- matrix( rep(0,100) )
tol  <- 10^-5
precision <- 1
iter <- 0

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

iter_4 = iter
precision_4 = precision
Dstar_4 = D(p)
Sstar_4 = S(p)
I_4 = sum(D(p))
Ip_4 =  Dstar_4*p

cbind("Precision is", precision_4)
cbind("Total number of rides is", I_4)
cbind("Average price of a ride is", sum(Ip_4))
cbind("Time taken on this machine, in number of iterations, is", iter_4)


