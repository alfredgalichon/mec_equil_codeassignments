rm(list=ls())

X <- seq(0.1, 1, 0.1)

grid <-merge(X,X)

nx <- matrix(1,10,10)
my <- apply(grid,1,prod)

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

dist_matrix <- matrix(0,100,100)

for (i in 1:100){
  for (j in 1:100){dist_matrix[i,j]= euc.dist(grid[i,], grid[j,])}
    }

uzx <-function(x){
  if(x>=0.099999){-10*x}
  else{0}
}

uzx <- apply(dist_matrix,1:2,uzx)

prob_zx <-function(p) {
  p <- matrix(t(p),100,100)
  exp_uzx= exp(uzx -  p)
  denom <- matrix(matrix(1 + rowSums(exp_uzx)),100,100)
  exp_uzx/denom
}

demand <- function(p) {
  d=matrix(rowSums(prob_zx(p)),100,1)
  d
}

czy = dist_matrix^2
prob_zy <- function(p){
  p <- matrix(t(p),100,100)
  exp_uzy= exp(p-czy)
  denom <- matrix(matrix(1 + rowSums(exp_uzy)),100,100)
  exp_uzy/denom
}

supply <- function(p){
  s=t(my)%*%prob_zy(p)
  matrix(s,100,1)
}

excess_supply <-function(p){supply(p)-demand(p)}

precision <- function(p){
  num = abs(supply(p)-demand(p))
  denom = supply(p)+demand(p)
  max(num/denom)
}

######## Gradient Descent ########

GradDescentStep <-function(p,epsilon){
  p-epsilon*excess_supply(p)
}

p0 = matrix(2,100,1)

GradDescentStep(p0,0.01)

GradDescent <-function(p,epsilon){
  initial_time = Sys.time()
  iter = 0
  while(precision(p)>= 10^(-5) & iter <= 5000){
    print(precision(p))
    aux_p = GradDescentStep(p,epsilon)
    p = aux_p
    iter = iter + 1
  }
  iter = iter -1
  final_time =  Sys.time()
  rbind(p, iter, precision(p), final_time-initial_time)
}

res = GradDescent(p0,1)
p_star=matrix(res[1:100],10,10)
p_star_long=matrix(res[1:100],100,1)
num_iter=res[101]
prec=res[102]
lapsed_time = res[103]

number_rides_demand = matrix(demand(p_star),10,10)
total_number_rides <- sum(number_rides_demand)
average_price_ride <- sum(number_rides_demand*p_star)/total_number_rides


######## Newton Descent ########

gradient_supply <-function(p){
  prob = prob_zy(p)
  product = my*t(prob)
  one_minus_prob =1-prob
  matrix = one_minus_prob%*%product
  diag_terms <- diag(diag(matrix))
  non_diag_terms =  -prob%*%product
  non_diag_terms = non_diag_terms-diag(diag(non_diag_terms))
  non_diag_terms+diag_terms
}

gradient_demand <-function(p){
  prob = prob_zx(p)
  product = nx*t(prob)
  one_minus_prob =-(1-prob)
  matrix = one_minus_prob%*%product
  diag_terms <- diag(diag(matrix))
  non_diag_terms =  prob%*%product
  non_diag_terms = non_diag_terms-diag(diag(non_diag_terms))
  non_diag_terms+diag_terms
}

gradient_excess_supply<-function(p){gradient_supply(p)-gradient_demand(p)}


NewtonDescentStep <-function(p,epsilon){
  grad = gradient_excess_supply(p)
  p-epsilon*solve(grad)%*%excess_supply(p)
}

p0 = rep(2,100)
NewtonDescentStep(p0,0.02)

NewtonDescent <-function(p,epsilon){
  initial_time = Sys.time()
  iter = 0
  while(precision(p)> 10^(-5) & iter <= 1000){
    print(precision(p))
    aux_p = NewtonDescentStep(p,epsilon)
    p <- aux_p
    iter = iter + 1
  }
  iter = iter -1
  final_time =  Sys.time()
  rbind(p, iter, precision(p), final_time-initial_time)
}

res2 = NewtonDescent(p0,0.2)
p_star_long2=matrix(res[1:100],100,1)

p_star2=matrix(res[1:100],10,10)
num_iter2=res[101]
prec2=res[102]
lapsed_time2 = res[103]

number_rides_demand2 = matrix(demand(p_star2),10,10)
total_number_rides2 <- sum(number_rides_demand2)
average_price_ride2 <- sum(number_rides_demand2*p_star2)/total_number_rides2




