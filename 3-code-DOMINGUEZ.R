rm(list=ls())


######### 3. Bellman - Ford algorithm

arcs <- read.csv("arcs.csv")
nodes <- read.csv("nodes.csv")

node_o <- 452
node_d <- 471

x <- matrix(arcs$from_stop_nb)
y <- matrix(arcs$to_stop_nb)
c <- arcs$dis_line

z <- nodes$stop_nb
s <- matrix(0, length(z),1)
s[node_o]=1
s[node_d]=-1


p_0 <- rep(exp(500),length(z))
p_0[node_o] <- 0
p <- rep(0,length(z))
x_p <- rep(0,length(z))

#Get prices 

gap <- 1
iter <- 0 
while ( gap> 0 ){
  
  for(i in c(1:length(z)) ){
    
    c_r <- c[y==i]
    x_r <- subset(x,y == i)
    p_r <- p_0[x_r]
    obj_r <- c_r+p_r
    p[i]<- min(p_0[i], min(obj_r))
    
    if (p[i]!=p_0[i]){
      x_p[i]<-x_r[which(obj_r==min(obj_r))] }
  }
  
  gap <- max(abs(p-p_0))
  p_0 <- p
  iter <- iter+1
  
  cat("Iteration = ", iter, ", Price gap = ", gap, '\n', sep='')
  
}

#Back out optimal path

current_z <- node_d
path <- current_z

while ( current_z!=node_o ){
  path<-cbind(x_p[current_z], path)
  current_z = x_p[current_z]
}

cat('The optimal path is = ', path)




  
  

