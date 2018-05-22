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

Hessian <- function(p,U,C, nX, nY){
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
  tol = 1e-5
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
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
    p = p-eps*(S-D)
  }
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  val = indutils(p,U,C,nX,nY)
  p_Gradient=p
  if (iter >= maxiter ) 
    {print('Maximum number of iterations reached in Direct Gradient Method.')
  } else {
    print(paste0("Direct gradient converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("Value of the minimization problem = ", val))
  }


### Gradient Method in optimization formulation
  ptm=proc.time()
  res = optim(p0,indutils,excess_supply,U=U, C=C, nX=nX, nY=nY, method="BFGS",control=list(maxit=maxiter))
  p=res$par
  
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  error= max(abs(rides_supply-rides_demand)/(rides_supply+rides_demand))
  p_Gradient_opt=p
  if (is.null(res$message)) {
    print(paste0("Gradient converged in ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
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
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
    H = Hessian(p,U,C, nX, nY)
    p = p- eps*as.vector(solve(H) %*% (S-D))
  }
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  val = indutils(p,U,C,nX,nY)
  p_Newton=p
  
  if (iter >= maxiter ) 
  {print('Maximum number of iterations reached in Newton Method.')
  } else {
    print(paste0("Newton method converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
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
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
    p = update(p,U,C,nX,nY)
  }

  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  val = indutils(p,U,C,nX,nY)
  p_Jacobi = p
  
  if (iter >= maxiter ) 
  {print('Maximum number of iterations reached in Jacobi Coordinate Descent Method.')
  } else {
    print(paste0("Jacobi Coordinate Descent Method converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("Value of the minimization problem = ", val))
  }

# Gauss Seidel Descent
  ptm=proc.time()
  cont = TRUE
  iter = 0
  p = p0
  
  while (cont){
    iter = iter + 1
    S = supply(p,C,nY)
    D = demand(p,U,nX)
    error = max(abs(S-D)/(S+D))
    if ((error < tol ) | (iter>=maxiter)){cont=FALSE}
      for (j in 1:dim(Z)[1]){
        p_new = update(p,U,C,nX,nY)
        p[j]=p_new[j]
      }
  }
  
  time = proc.time()-ptm
  rides_demand = sum(demand(p,U,nX))
  rides_supply = sum(supply(p,C,nY))
  val = indutils(p,U,C,nX,nY)
  p_Gauss=p
  if (iter >= maxiter ) 
  {print('Maximum number of iterations reached in Gauss-Seidel Coordinate Descent Method.')
  } else {
    print(paste0("Gauss-Seidel Coordinate Descent Method converged in ",iter, " steps and ", time[1], "s."))
    print(paste0("Precision error is = ", error))
    print(paste0("The total number of rides demand = ", rides_demand))
    print(paste0("The total number of rides supply = ", rides_supply))
    print(paste0("Value of the minimization problem = ", val))
  }

# Summary of Previous findings:
  # [1] "Direct gradient converged in 297 steps and 0.298000000000002s."
  # [1] "Precision error is = 9.79475702593096e-07"
  # [1] "The total number of rides demand = 30.0742966681671"
  # [1] "The total number of rides supply = 30.0742553546667"
  # [1] "Value of the minimization problem = 191.668911118195"
  # 
  # [1] "Gradient converged in 0.0190000000000055s."
  # [1] "Precision error is = 2.83328655085832e-05"
  # [1] "The total number of rides demand = 30.0759458634938"
  # [1] "The total number of rides supply = 30.074241636321"
  # [1] "Value of the minimization problem = 191.668911202612"
  # 
  # [1] "Newton method converged in 61 steps and 4.0569999999999s."
  # [1] "Precision error is = 9.74517195854596e-07"
  # [1] "The total number of rides demand = 30.0742230325944"
  # [1] "The total number of rides supply = 30.0742559679333"
  # [1] "Value of the minimization problem = 191.668911118181"
  # 
  # [1] "Jacobi Coordinate Descent Method converged in 31 steps and 0.0630000000001019s."
  # [1] "Precision error is = 8.38923706543084e-07"
  # [1] "The total number of rides demand = 30.0742873849622"
  # [1] "The total number of rides supply = 30.0742554314715"
  # [1] "Value of the minimization problem = 191.668911118178"
  
  # [1] "Gauss-Seidel Coordinate Descent Method converged in 22 steps and 1.58600000000001s."
  # [1] "Precision error is = 6.49126557689387e-07"
  # [1] "The total number of rides demand = 30.0742670640525"
  # [1] "The total number of rides supply = 30.0742556010412"
  # [1] "Value of the minimization problem = 191.668911118158"

## Optional: Change Supply Function



  
  
  
