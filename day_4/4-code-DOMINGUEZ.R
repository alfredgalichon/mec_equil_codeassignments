rm(list=ls())
library(tictoc)

######### 1 and 2. Taxless benchmark and Linear taxes

#Set up grids

coordinates = seq(0.1,1,0.1)

X <- merge(coordinates, coordinates)
Y <- merge(coordinates, coordinates)

m_y  <- matrix(Y$x*Y$y)
m_yM <- matrix(m_y,100,100)

#Distance and utility function

d <- matrix(0,100,100)

for (i in c(1:100)){
  for (j in c(1:100)){
    d[i,j] <- sqrt( (X$x[i] - Y$x[j])^2 + (X$y[i] - Y$y[j])^2 )}
}

d2 <- d*d
alpha <- -0.1*d2
#alpha <- -0.1*d2/0.8 #Change here for linear taxes

gamma <- d
K <- exp( 0.5*(alpha+gamma) )

obj_x <- function(mu_xi){
   return( sum( sqrt(mu_xi*mu_0y)*K[i,] ) + mu_xi -1 )
}

obj_y <- function(mu_yi){
  return( sum( sqrt(mu_x*mu_yi)*K[,j] ) + mu_yi - m_y[j] )
}


# Solve system of equations for mu_x0 and mu_0y

mu_0y <- m_y
mu_y <- m_y
mu_0x <- matrix( rep(1,100) )
mu_x <- mu_0x 
obj_x_val <- rep(1,100)
obj_y_val <- rep(1,100)
tol <- 10^(-5)
gap <- 1
iter <-0

while( gap > tol ){
  
  for( i in c(1:100) ){
    mu_x[i] <- uniroot(obj_x, c(0,1), tol=10^(-18))$root
    obj_x_val[i] <- obj_x( mu_x[i] )
  }
  
  for( j in c(1:100) ){
    mu_y[j] <- uniroot(obj_y, c(0,1), tol=10^(-18))$root
    obj_y_val[j] <- obj_y( mu_y[j] )
  }
  
  gap <- max( max( abs( mu_0x - mu_x ) ) , max( abs( mu_0y - mu_y ) ) )
  
  mu_0x <- mu_x
  mu_0y <- mu_y
  
  iter <- iter + 1
  
  cat("Iteration", iter, "Gap", gap, "\n")

}

mu_xy <- matrix(0,100,100)
w_xy <- matrix(0,100,100)

for( i in c(1:100) ){
  
  for(j in c(1:100) ){
    mu_xy[i,j] <- sqrt( mu_0x[i]*mu_0y[j] )*K[i,j]
    w_xy[i,j] <- gamma[i,j] - log( mu_xy[i,j]/mu_0x[i] )
  }
  
}
  
cat('The total number of matched pairs is', sum(mu_xy))
cat('The average wage is', sum(mu_xy*w_xy)/sum(mu_xy) )


