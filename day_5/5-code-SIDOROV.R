
############################
###########Question 1#######
############################

library(gurobi)
Nf<-2 #Number of firms
Nw<-5 #Number of workers

fi<-matrix(runif(Nw*Nf)+1,Nw,Nf)

#Nf<-dim(fi)[2]
#Nw<-dim(fi)[1]
pset<-matrix(c(1,0),2,1)
for(i in 2:Nw){
 pset<- rbind(cbind(matrix(1,dim(pset)[1],1),pset),cbind(matrix(0,dim(pset)[1],1),pset))
}
Fi<-(pset %*% (fi^2))^(0.5) #productivity of each bundle of workers
Nbun<-dim(Fi)[1]

###########################
#Algorithm
###########################
minwage<-matrix(0,Nw,Nf)
for(iter in 1:2000){
Profits<-Fi-pset%*% minwage
choice<-t(matrix(rep(1,Nbun),Nbun,1)%*%apply(Profits, 2, max)==Profits)
if(sum(apply(choice, 1, sum))!=Nf){  stop("More than one choice for a firm") }
offerees<-(choice %*% pset) #workers that get offers
offers<-offerees*t(minwage) #wage offers
prereply<-(matrix(rep(1,Nf),Nf,1)%*%apply(offers, 2, max)==offers)*matrix(rep(1:Nf,Nw),Nf,Nw)
reply<-matrix(rep(1,Nf),Nf,1)%*%apply(prereply, 2, max)==prereply
minwage<-minwage+t(offerees*(1-reply))*0.001
Ndecl<-sum(offerees*(1-reply))
#print(Ndecl)
if(Ndecl==0){
  print(iter)
  #print(minwage)
  #print("Individual worker productivities(rows are firms)")
  #print(t(fi))
  #print("Wages of the workers (rows are firms)")
  #print((offerees)*t(minwage))
  break
}
}

#gurobi
###########################
A<- cbind(kronecker(matrix(1,Nf,1),pset),kronecker(diag(Nf),matrix(1,Nbun,1)))
rhs<-matrix(Fi,Nbun*Nf,1)
obj<-matrix(1,1,Nw+Nf)
result		<- gurobi (list(A=A,obj=obj,modelsense="min",rhs=rhs,sense=">="), params=NULL ) 
if (result$status=="OPTIMAL") { } else {stop("optimization problem with Gurobi")}
pivec <- result$x; 
Lvec <- matrix((result$pi),Nbun,Nf)
###printing

if(Ndecl==0){
  #print(iter)
  #print(minwage)
  print("Individual worker productivities(rows are firms)")
  print(t(fi))
  print("Wages of the workers (rows are firms)")
  print((offerees)*t(minwage))
  #break
}

print("Wages of the workers in gurobi case")
print((t(Lvec)%*%pset)*(matrix(rep(1,Nf),Nf,1)%*%pivec[1:Nw]))

