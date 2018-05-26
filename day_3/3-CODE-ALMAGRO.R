#### NYC Subway ####
rm(list=ls())

# Reading data
arcs  <- read.csv("~/Desktop/CEcourse/Material/NYCSubway/arcs.csv", header = TRUE)
nodes <- read.csv("~/Desktop/CEcourse/Material/NYCSubway/nodes.csv", header = TRUE)


origin <- arcs$from_stop_nb
destination <- arcs$to_stop_nb
distance <-arcs$dis_line
N<-length(origin)

node_destination <- 471
node_origin <- 452

id <- (nodes$stop_id)
n<-length(id)
sz <- matrix(0,n,1)
sz[node_origin]=1
sz[node_destination]=-1


#### Bellman-Ford algorithm ####

price <- matrix(exp(15),n,1)
aux_price = price
price[node_origin]=0
predecesor<-rep(0,n)

tol = 10^(-5)
max_iter = 1000

err = 1
iter = 0
while(err>tol && iter< max_iter){
  for (i in 1:n){
    connected_to<- subset(origin,destination == i)
    distance_connected<-distance[destination==i]
    price_connected<-price[connected_to]
    c_p_connected <- distance_connected+price_connected
    min_c_p <- min(c_p_connected)
    aux_price[i]<- min(price[i],min_c_p)
    if (aux_price[i]!=price[i]){
      print(i)
      predecesor[i]<-connected_to[which(c_p_connected==min_c_p)]
    }
  }
  err = sum(abs(aux_price-price))
  print(which.min(price_connected))
  price<-aux_price
  iter = iter +1
  cat('Iteration = ', iter,', error = ', err,'\n',sep='')
}


node = node_destination
shortest_path <- node
iter = 1
while (node!=node_origin && iter <max_iter){
  shortest_path<-cbind(shortest_path,predecesor[node])
  node = predecesor[node]
  iter = iter +1
  print(node)
}

cat('Shortest path = ', rev(shortest_path))


#### Regularized Problem ####
library("nloptr")

nabla = matrix(0,N,n)

for (i in 1:n){
  connected_to<- which(destination == i)
  connected_from<-which(origin == i)
  nabla[connected_to,i]=1
  nabla[connected_from,i]=-1
}

regular_flow <- function(mu,a){
  return(t(distance)%*%mu+t(log(mu))%*%mu)
}

grad_regular_flow <- function(mu,a){
  return(distance+log(mu)+1)
}

constraint <- function(mu,a){
  return( t(nabla)%*%mu-sz )
}