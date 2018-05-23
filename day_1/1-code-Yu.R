# Homework 1 for Math and Econ Code - Equilibrium
# Yue Yu, yy2558@columbia.edu

rm(list=ls())
setwd("C:/Users/Yue/Dropbox/temp/MathEcon/equilibrium/assignment 1")
con <- file("results.log")
sink(con, append=FALSE, type=c("output", "message"))

#1. Model setup
n = 10
X1 = rep(0.1*(1:10),n)
X2 = c(matrix(rep(0.1*(1:10),n), nrow = n, ncol = n, byrow = TRUE)) 
nodes = cbind(X1,X2)
nX = rep(1,n^2)
mY = nodes[,1]*nodes[,2]
C = matrix(0,n^2, n^2)
for (i in (1:n^2)) {
  for (j in (i:n^2)) {
    C[i,j] = (nodes[i,1] - nodes[j,1])^2 + (nodes[i,2] - nodes[j,2])^2
    C[j,i] = C[i,j]
  }
}
dist = sqrt(C)
U = -10*(dist > 0.09999999999999)
P0 = rep(0,n^2)

#2. Define functions
## Dz(P)
Demand <- function(P,U,nX) {
  Uxz = exp(t(t(U)-P))
  D = colSums(Uxz*nX/(rowSums(Uxz)+1)) 
  return(D)
}

## Sz(P)
Supply <- function(P,C,mY) {
  Cyz = exp(t(-t(C)+P))
  S = colSums(Cyz*mY/(rowSums(Cyz)+1)) 
  return(S)
}


#3. Gradient descent method
tol = 1e-5
epsilon = 0.2
criteria = 10
iter = 0
iter_max = 5000
P_old = P0
Sp_old = Supply(P0,C,mY)
Dp_old = Demand(P0,U,nX)
start_time <- Sys.time()
while (criteria> tol) {
  iter = iter + 1
  P_new = P_old - epsilon*(Sp_old-Dp_old)
  Sp_new = Supply(P_new,C,mY)
  Dp_new = Demand(P_new,U,nX)
  criteria = max(abs(Sp_new-Dp_new)/(Sp_new+Dp_new))
  Sp_old = Sp_new
  Dp_old = Dp_new 
  P_old = P_new
  #print(iter)
  if (iter>=iter_max) break
}
end_time = Sys.time()
time_length = end_time - start_time
sumIz = sum(Dp_new)
avgP = mean(Dp_new*P_new)  

print("Gradient descent method results")
print(sprintf("(i) Precision number: %g", criteria)) 
print(sprintf("(ii) Total number of rides: %f", sumIz)) 
print(sprintf("(iii) Average price of ride: %f", avgP)) 
print(sprintf("(iv) Number of iteration: %g", iter)) 
print(sprintf("(v) Time length: %f", time_length))
print("-----------------------------------------")

#4. Newton descent method
Hessian <- function(P,U,C,nX,mY,n) {
  JP = matrix(0,n^2,n^2)
  Uxz = exp(t(t(U)-P))
  MUxz = Uxz/(rowSums(Uxz)+1)
  Cyz = exp(t(P-t(C)))
  MUyz = Cyz/(rowSums(Cyz)+1)
  
  for (i in 1:n^2){
    for (j in 1:n^2){
      JP[i,j]= -sum(nX*MUxz[,i]*MUxz[j])-sum(mY*MUyz[,i]*MUyz[,j])
    }
  }
  for (i in 1:n^2){
    JP[i,i]=sum(nX*MUxz[,i]*(1-MUxz[,i]))+sum(mY*MUyz[,i]*(1-MUyz[,i]))
  }
  return(JP)
}


tol = 1e-5
epsilon = 0.2
criteria = 10
iter = 0
iter_max = 5000
P_old = P0
Sp_old = Supply(P0,C,mY)
Dp_old = Demand(P0,U,nX)
start_time <- Sys.time()
while (criteria> tol) {
  iter = iter + 1
  JP = Hessian(P_old,U,C,nX,mY,n)
  P_new = c(as.matrix(P_old - epsilon*solve(JP)%*%(Sp_old-Dp_old)))
  Sp_new = Supply(P_new,C,mY)
  Dp_new = Demand(P_new,U,nX)
  criteria = max(abs(Sp_new-Dp_new)/(Sp_new+Dp_new))
  Sp_old = Sp_new
  Dp_old = Dp_new 
  P_old = P_new
  #print(iter)
  if (iter>=iter_max) break
}
end_time = Sys.time()
time_length = end_time - start_time
sumIz = sum(Dp_new)
avgP = mean(Dp_new*P_new)  

print("Newton descent method results")
print(sprintf("(i) Precision number: %g", criteria)) 
print(sprintf("(ii) Total number of rides: %f", sumIz)) 
print(sprintf("(iii) Average price of ride: %f", avgP)) 
print(sprintf("(iv) Number of iteration: %g", iter)) 
print(sprintf("(v) Time length: %f", time_length)) 
print("-----------------------------------------")

#5. Jacobi coordinate updates
Jacobi <- function(P,U,C,nX,mY){
  Uxz = nX*exp(U)/(rowSums(exp(t(t(U)-P)))+1)
  Cyz = mY*exp(-C)/(rowSums(exp(t(P-t(C))))+1)
  P1 = 1/2*(log(colSums(Uxz))-log(colSums(Cyz)))
  return(P1)
}

tol = 1e-5
epsilon = 0.2
criteria = 10
iter_max = 5000
iter = 0
P_old = P0
Sp_old = Supply(P0,C,mY)
Dp_old = Demand(P0,U,nX)
start_time <- Sys.time()
while (criteria> tol) {
  iter = iter + 1
  P_new = Jacobi(P_old,U,C,nX,mY)
  Sp_new = Supply(P_new,C,mY)
  Dp_new = Demand(P_new,U,nX)
  criteria = max(abs(Sp_new-Dp_new)/(Sp_new+Dp_new))
  Sp_old = Sp_new
  Dp_old = Dp_new 
  P_old = P_new
  #print(iter)
  if (iter>=iter_max) break
}
end_time = Sys.time()
time_length = end_time - start_time
sumIz = sum(Dp_new)
avgP = mean(Dp_new*P_new)  

print("Jacobi method results")
print(sprintf("(i) Precision number: %g", criteria)) 
print(sprintf("(ii) Total number of rides: %f", sumIz)) 
print(sprintf("(iii) Average price of ride: %f", avgP)) 
print(sprintf("(iv) Number of iteration: %g", iter)) 
print(sprintf("(v) Time length: %f", time_length)) 
print("-----------------------------------------")
#6. Gauss Seidel Coordinate updates
tol = 1e-5
epsilon = 0.2
criteria = 10
iter_max = 5000
iter = 0
P_old = P0
Sp_old = Supply(P0,C,mY)
Dp_old = Demand(P0,U,nX)
start_time <- Sys.time()
while (criteria> tol) {
  iter = iter + 1
  for (i in 1:n^2){
    P_new = Jacobi(P_old,U,C,nX,mY)
    P_old[i]=P_new[i]
  }
  Sp_new = Supply(P_new,C,mY)
  Dp_new = Demand(P_new,U,nX)
  criteria = max(abs(Sp_new-Dp_new)/(Sp_new+Dp_new))
  Sp_old = Sp_new
  Dp_old = Dp_new 
  P_old = P_new
  #print(iter)
  if (iter>=iter_max) break
}
end_time = Sys.time()
time_length = end_time - start_time
sumIz = sum(Dp_new)
avgP = mean(Dp_new*P_new)  

print("Gauss Seidel results")
print(sprintf("(i) Precision number: %g", criteria)) 
print(sprintf("(ii) Total number of rides: %f", sumIz)) 
print(sprintf("(iii) Average price of ride: %f", avgP)) 
print(sprintf("(iv) Number of iteration: %g", iter)) 
print(sprintf("(v) Time length: %f", time_length)) 
print("-----------------------------------------")
sink(file = NULL)
sink()