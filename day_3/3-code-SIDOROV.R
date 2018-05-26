library('Matrix')
library('slam')
library('matrixcalc')
library('gurobi')
gtime0 = proc.time()
nodes <- read.csv("nodes.csv")
arcs <- read.csv("arcs.csv")
goal<- matrix(0,length(nodes$stop_name),1)
start<- 452
finish <- 471
goal[finish] <-1 
goal[start]<- (-1)
nabt<- matrix(0,length(nodes$stop_name),length(arcs$dis_line))
for (i in 1:length(arcs$dis_line)){
  nabt[arcs$from_stop_nb[i],i]<- (-1)
  nabt[arcs$to_stop_nb[i],i]<- 1
}

params <- list(OutputFlag=0 )
result <- gurobi (list(A=nabt,obj=arcs$dis_line,modelsense="min",rhs=goal,sense="="), params ) 
t<- (1:length(result$x)) * result$x
mixedpath<-arcs[t[t>0],1:2]
current<-start
path<-vector(mode =  "character",length =length(t[t>0])+1 )
path[1]<- as.character(nodes$stop_name[start])
len<-length(t[t>0])
for(i in  1:len){
  for (j in 1:len){
    if(mixedpath[j,1]==current){
      path[i+1]<-as.character(nodes$stop_name[mixedpath[j,2]])
      new<-mixedpath[j,2]
      break
    }
  }
  current<-new
  if(new==finish){break}
}
gurobipath<-path[1:(i+1)]
print(gurobipath)
gtime1 = proc.time()
#####################################
############################Bellman- Ford
bftime0 = proc.time()
dist<-arcs$dis_line
numnods<- length(nodes$stop_name) #number of nodes
Large<-100*sum(dist) #100*sum(dist) is infinity ;-)
p0<-matrix(Large,numnods,1) 
start<- 452
finish <- 471
p0[start]<- 0
p<-p0
########## martix M: for two nodes tells the arc number
M<- matrix(0,numnods,numnods)
for (i in 1:length(arcs$dis_line)){
  if(M[arcs$from_stop_nb[i],arcs$to_stop_nb[i]]==0){
    M[arcs$from_stop_nb[i],arcs$to_stop_nb[i]]<- i}
  else{ if (dist[M[arcs$from_stop_nb[i],arcs$to_stop_nb[i]]]<dist[i]){
    M[arcs$from_stop_nb[i],arcs$to_stop_nb[i]]<- i
    print("Hey, it's a better ark!")
  }}
}

for (iter in 1:1000){
  psearch<-p
  for(z in 1:numnods){
    for (zpr in 1:numnods){
      #psearch[z]<-p[z]
      if(M[zpr,z]>0){
        if(dist[M[zpr,z]]+p[zpr]<psearch[z]){
          psearch[z]<-dist[M[zpr,z]]+p[zpr]
        }
      }
    }
  }
  precision<- max(abs(psearch-p))
  if(precision<10^(-6)){
    break
  }
  p<-psearch
}

#########recovering the path
res<-matrix(0,length(dist),1) 
for(i in 1:length(dist)){
  res[i]<-abs(p[arcs$to_stop_nb[i]]-p[arcs$from_stop_nb[i]]-dist[i])<10^(-5)}
t<- (1:length(dist)) * res
mixedpath<-arcs[t[t>0],1:2]


repath<-vector(mode =  "character",length =length(t[t>0])+1 ) #we will go in reverse: from the last to the first point
current<-finish
repath[1]<- as.character(nodes$stop_name[finish])
len<-length(t[t>0])
for(i in  1:len){
  for (j in 1:len){
    if(mixedpath[j,2]==current){
      repath[i+1]<-as.character(nodes$stop_name[mixedpath[j,1]])
      new<-mixedpath[j,1]
      break
    }
  }
  current<-new
  print(new)
  
  if(new==start){break}
}
repath<-repath[1:(i+1)]
BellmanFordpath<-rev(repath)
print(BellmanFordpath)
bftime1 = proc.time()
print(gtime1-gtime0)
print(bftime1-bftime0)
