### HW for Math+Econ+Code Day 1 
## Danyan Zha
## May 21, 2018

rm(list=ls())

n <- 10
X1 = rep(0.1*c(1:10), n)
X2 = rep(0.1*c(1:10),each=n)
X = cbind(X1,X2)
Y = X
Z = X

nX = rep(1,dim(X)[1])
nY = Y[,1]*Y[,2] 

dist = matrix(0,n^2, n^2)
for (i in 1:n^2){
  for (j in 1:n^2){
    dist[i,j]=sqrt(sum((X[i,]-Z[j,])^2))
  }
}

U = -10 *(dist-0.099999999>0)
C = dist^2

# 
demand <- function(p,U,nX){
  A_xz = exp(t(t(U)-p))
  B_x = 1 + rowSums(A_xz)
  C_xz = nX*A_xz/B_x
  return(colSums(C_xz))
}

supply <- function(p,C,nY){
  A_yz = exp(t(p-t(C)))
  B_y = 1 + rowSums(A_yz)
  C_yz = nY*A_yz/B_y
  return(colSums(C_yz))
}

indutils <- function(p,U,C,nX,nY){
  A_xz = exp(t(t(U)-p))
  B_x = log(1 + rowSums(A_xz))
  
  A_yz = exp(t(p-t(C)))
  B_y = log(1 + rowSums(A_yz))
  
  return(sum(nX*B_x)+sum(nY*B_y))
}

excess_supply <- function(p,U,C,nX,nY){
  return(supply(p,C,nY)-demand(p,U,nX))
}

Hessian <- function(p,U,C,nX,nY){
  DE <- matrix(0,n^2, n^2)
  A_xz = exp(t(t(U)-p))
  B_x = 1 + rowSums(A_xz)
  mu_xz = A_xz/B_x
  
  A_yz = exp(t(p-t(C)))
  B_y = 1 + rowSums(A_yz)
  mu_yz <- A_yz/B_y
  
  for (i in 1:n^2){
    for (j in 1:n^2){
      DE[i,j]= -sum(nX*mu_xz[,i]*mu_xz[j])-sum(nY*mu_yz[,i]*mu_yz[,j])
    }
  }
  for (i in 1:n^2){
    DE[i,i]=sum(nX*mu_xz[,i]*(1-mu_xz[,i]))+sum(nY*mu_yz[,i]*(1-mu_yz[,i]))
  }
  return(DE)
}

## Direct Gradient Method
  tol = 1e-6
  maxiter <- 1000000
  p0 = rep(0, n^2)
  eps = 0.2
  
  ptm=proc.time()
  cont = TRUE
  iter = 0
  p = p0
  
  while (cont){
    iter = iter + 1 
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    p = p-eps*(S-D)
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
  }
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  avgprice = sum(demand(p,U,nX)*p)/rides_demand
  val = indutils(p,U,C,nX,nY)
  p_Gradient=p
  if (iter >= maxiter ) 
    {print('Maximum number of iterations reached in Direct Gradient Method.')
  } else {
    print(paste0("Direct gradient converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("The average price of a ride is = ", avgprice))
    print(paste0("Value of the minimization problem = ", val))
  }


### Gradient Method in optimization formulation
  ptm=proc.time()
  res = optim(p0,indutils,excess_supply,U=U, C=C, nX=nX, nY=nY, method="BFGS",control=list(maxit=maxiter))
  p=res$par

  time = proc.time()-ptm
  S = supply(p,C,nY)
  D = demand(p,U,nX)
  rides_demand = sum(D)
  rides_supply = sum(S)
  avgprice = sum(demand(p,U,nX)*p)/rides_demand
  error = max(abs(S-D)/(S+D))
  p_Gradient_opt=p
  if (is.null(res$message)) {
    print(paste0("Gradient converged in ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("The average price of a ride is = ", avgprice))
    print(paste0("Value of the minimization problem = ", res$value))
  } else {
    print(res$message)
    }


## Newton descent
  ptm=proc.time()
  cont = TRUE
  iter = 0
  p = p0
  
  while (cont){
    iter = iter + 1
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    H = Hessian(p,U,C,nX,nY)
    p = p- eps*as.vector(solve(H) %*% (S-D))
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
  }
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  avgprice = sum(demand(p,U,nX)*p)/rides_demand
  val = indutils(p,U,C,nX,nY)
  p_Newton=p
  
  if (iter >= maxiter ) 
  {print('Maximum number of iterations reached in Newton Method.')
  } else {
    print(paste0("Newton method converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("The average price of a ride is = ", avgprice))
    print(paste0("Value of the minimization problem = ", val))
  }


## Coordinate Descent
# Jacobi Descent
  update <- function(p,U,C,nX,nY){
    A_xz = exp(t(t(U)-p))
    B_x = 1 + rowSums(A_xz)
    C_xz = nX*exp(U)/B_x
    
    A_yz = exp(t(p-t(C)))
    B_y = 1 + rowSums(A_yz)
    C_yz = nY*exp(-C)/B_y
    
    p_new = 1/2*(log(colSums(C_xz))-log(colSums(C_yz)))
    return(p_new)
  }

  ptm=proc.time()
  cont = TRUE
  iter = 0
  p = p0

  while (cont){
    iter = iter + 1 
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    p = update(p,U,C,nX,nY)
    S = supply(p,C,nY)
    D = demand(p,U,nX) 
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
  }

  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  avgprice = sum(demand(p,U,nX)*p)/rides_demand
  val = indutils(p,U,C,nX,nY)
  p_Jacobi = p
  
  if (iter >= maxiter ) 
  {print('Maximum number of iterations reached in Jacobi Coordinate Descent Method.')
  } else {
    print(paste0("Jacobi Coordinate Descent Method converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("The average price of a ride is = ", avgprice))
    print(paste0("Value of the minimization problem = ", val))
  }

# Gauss Seidel Descent
  ptm=proc.time()
  cont = TRUE
  iter = 0
  p = p0
  
  while (cont){
    iter = iter + 1
    for (j in 1:dim(Z)[1]){
      p_new = update(p,U,C,nX,nY)
      p[j]=p_new[j]
    }
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
  }
  
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  avgprice = sum(demand(p,U,nX)*p)/rides_demand
  val = indutils(p,U,C,nX,nY)
  p_Gauss=p
  if (iter >= maxiter ) 
  {print('Maximum number of iterations reached in Gauss-Seidel Coordinate Descent Method.')
  } else {
    print(paste0("Gauss-Seidel Coordinate Descent Method converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("The average price of a ride is = ", avgprice))
    print(paste0("Value of the minimization problem = ", val))
  }

# Summary of Previous findings:
  # [1] "Direct gradient converged in 296 steps and 0.325000000000045s."
  # [1] "Precision error is = 9.79475702593096e-07"
  # [1] "The total number of rides demand = 30.0742984777543"
  # [1] "The total number of rides supply = 30.0742553396218"
  # [1] "The average price of a ride is = 0.844430416585444"
  # [1] "Value of the minimization problem = 191.668911118198"
  # 
  # [1] "Gradient converged in 0.0199999999999818s."
  # [1] "Precision error is = 7.16272517309403e-05"
  # [1] "The total number of rides demand = 30.0759458634938"
  # [1] "The total number of rides supply = 30.074241636321"
  # [1] "The average price of a ride is = 0.844352148277277"
  # [1] "Value of the minimization problem = 191.668911202612"
  # 
  # [1] "Newton method converged in 60 steps and 4.12900000000002s."
  # [1] "Precision error is = 9.74517195854596e-07"
  # [1] "The total number of rides demand = 30.0742149201625"
  # [1] "The total number of rides supply = 30.0742560356446"
  # [1] "The average price of a ride is = 0.844434401869166"
  # [1] "Value of the minimization problem = 191.668911118195"
  # 
  # [1] "Jacobi Coordinate Descent Method converged in 30 steps and 0.0710000000000264s."
  # [1] "Precision error is = 8.38923706543084e-07"
  # [1] "The total number of rides demand = 30.0743045856514"
  # [1] "The total number of rides supply = 30.0742552882611"
  # [1] "The average price of a ride is = 0.844430117608689"
  # [1] "Value of the minimization problem = 191.668911118212"
  # 
  # [1] "Gauss-Seidel Coordinate Descent Method converged in 21 steps and 1.40899999999999s."
  # [1] "Precision error is = 6.49126557689387e-07"
  # [1] "The total number of rides demand = 30.0742777456828"
  # [1] "The total number of rides supply = 30.0742555124653"
  # [1] "The average price of a ride is = 0.844431406422646"
  # [1] "Value of the minimization problem = 191.668911118167"
  
## Optional: Change Supply Function
  



  
  
  
