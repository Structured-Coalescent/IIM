

CI.prof.likelihood<-function(params1,i){
ci<-vector(length=2)
params<-params1[-i]
params[1]<-params[1]+0.25 		#to avoid "false convergence"

	 

##lower bound search


points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10
#the "round" function avoids, 
#for instance, truncating 2.3*100
#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1-j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni.bi.2,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni.bi.2,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni.bi.2,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[1]<-points[j-1]

##upper bound search:    

points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10
 									 									#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1+j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni.bi.2,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni.bi.2,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni.bi.2,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[2]<-points[j-1]

ci
} 

CI.prof.likelihood(c(3.2731545,3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,9.7776228,
0.0000000,0.06416968),9)


params1<-c(3.2731545,3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,9.7776228,
0.0000000,0.06416968)



negll.IIM.uni.bi.2(c(3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,9.7776228,
0.0000000,0.06416968),3.2731545,8)

negll.IIM.uni.bi.2<-function(params,target,i){

M<-vector(length=2,mode="numeric")

if(i==1){
a<-target/params[1]
theta<-params[1]
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==2){
a<-params[1]/target
theta<-target
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==3){
a<-params[1]/params[2]
theta<-params[2]
b<-target/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==4){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(target/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==5){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,target/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==6){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-target/theta
tau0<-params[6]/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==7){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-target/theta+tau1
M[1]<-params[7]*2
M[2]<-params[8]*2/b
}

if(i==8){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-params[7]/theta+tau1
M[1]<-target*2
M[2]<-params[8]*2/b
}

if(i==9){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-params[7]/theta+tau1
M[1]<-params[8]*2
M[2]<-target*2/b
}




if (M[1]>0 & M[2]>0){

R1<-matrix(rep(r.a,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1a,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,3),nrow=3,byrow=TRUE)


Qt<-matrix(ncol=4,nrow=4)
Qt[,1]<-c(-(1+M[1]),0,M[1],1)
Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
Qt[,4]<-c(0,0,0,0)
X<-eigen(Qt)$vectors
nu<-eigen(Qt)$values[-4]
P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
C<-solve(X) #(%*%P.0)
A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;



loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)

loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)

loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
}

if (M[1]==0 & M[2]>0){

R1<-matrix(rep(r.a,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1a,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,2),nrow=2,byrow=TRUE)


A1<-1
nu1<--1

A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
nu21<--1
A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
nu22<--M[2]/2
A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
nu23<--(1/b+M[2])

A2<-c(A21,A22,A23)
nu2<-c(nu21,nu22,nu23)

A31<-M[2]/(M[2]-2)
nu31<--1
A32<-2/(2-M[2])
nu32<--M[2]/2

A3<-c(A31,A32)
nu3<-c(nu31,nu32)

loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)

loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)

loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)

}

if (M[1]>0 & M[2]==0){

R1<-matrix(rep(r.a,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
R3<-matrix(rep(r.bc,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1a,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
X3<-matrix(rep(x3bc,2),nrow=2,byrow=TRUE)

A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
nu11<--1/b
A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
nu12<--M[1]/2
A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
nu13<--(1+M[1])

A1<-c(A11,A12,A13)
nu1<-c(nu11,nu12,nu13)

A2<-1
nu2<--1/b

A31<-b*M[1]/(b*M[1]-2)
nu31<--1/b
A32<-2/(2-b*M[1])
nu32<--M[1]/2

A3<-c(A31,A32)
nu3<-c(nu31,nu32)

loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)

loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)

loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)

}

if (M[1]==0 & M[2]==0){

R1<-matrix(rep(r.a,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1a,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,3),nrow=3,byrow=TRUE)

loglike1<-sum(log(
(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
+
exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
   ))


loglike2<-sum(log(
(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
+
exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
   ))


loglike3<-sum(log(
exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
))

}


-(loglike1+loglike2+loglike3)
}

############# with a unidirectional likelihood function

CI.prof.likelihood.2<-function(params1,i){
ci<-vector(length=2)
params<-params1[-i]
params[1]<-params[1]+0.25 		#to avoid "false convergence"

	 

##lower bound search


points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10	#the "round" function avoids, 
     									#for instance, truncating 2.3*100 									#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1-j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[1]<-points[j-1]

##upper bound search:    

points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10
 									 									#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1+j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[2]<-points[j-1]

ci
} 

CI.prof.likelihood.2(c(3.2731545,3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,9.7776228,
0.06416968),7)

negll.IIM.uni8.analytic.1<-function(params,target,i){


if(i==1){
a<-target/params[1]
theta<-params[1]
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*2/b
}

if(i==2){
a<-params[1]/target
theta<-target
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*2/b
}

if(i==3){
a<-params[1]/params[2]
theta<-params[2]
b<-target/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*2/b
}

if(i==4){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(target/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*2/b
}

if(i==5){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,target/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*2/b
}

if(i==6){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-target/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*2/b
}

if(i==7){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-target/theta+tau1
M2<-params[7]*2/b
}


if(i==8){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-params[7]/theta+tau1
M2<-target*2/b
}



if (M2>0){

R1<-matrix(rep(r.a,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1a,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,2),nrow=2,byrow=TRUE)


A1<-1
nu1<--1

A21<-(b*M2^2)/((-2 + M2)*(1-b+b*M2))
nu21<--1
A22<-4*b*M2/((2-M2)*(2+b*M2))
nu22<--M2/2
A23<-(1/b)/(1/b+M2)+b^2*M2^2/((2+b*M2)*(1-b+b*M2)*(1/b+M2))
nu23<--(1/b+M2)

A2<-c(A21,A22,A23)
nu2<-c(nu21,nu22,nu23)

A31<-M2/(M2-2)
nu31<--1
A32<-2/(2-M2)
nu32<--M2/2

A3<-c(A31,A32)
nu3<-c(nu31,nu32)

loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)

loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)

loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)

}


if (M2==0){

R1<-matrix(rep(r.a,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1a,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,3),nrow=3,byrow=TRUE)

loglike1<-sum(log(
(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
+
exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
   ))


loglike2<-sum(log(
(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
+
exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
   ))


loglike3<-sum(log(
exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
))

}


-(loglike1+loglike2+loglike3)
}

##### For T0 instead of V:

negll.IIM.uni8.analytic.1.T0<-function(params,target,i){


if(i==1){
a<-target/params[1]
theta<-params[1]
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta
M2<-params[7]*2/b
}

if(i==2){
a<-params[1]/target
theta<-target
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta
M2<-params[7]*2/b
}

if(i==3){
a<-params[1]/params[2]
theta<-params[2]
b<-target/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta
M2<-params[7]*2/b
}

if(i==4){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(target/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta
M2<-params[7]*2/b
}

if(i==5){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,target/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta
M2<-params[7]*2/b
}

if(i==6){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-target/theta
tau0<-params[6]/theta
M2<-params[7]*2/b
}

if(i==7){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-target/theta
M2<-params[7]*2/b
}


if(i==8){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-params[7]/theta
M2<-target*2/b
}



if (M2>0){

R1<-matrix(rep(r.a,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1a,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,2),nrow=2,byrow=TRUE)


A1<-1
nu1<--1

A21<-(b*M2^2)/((-2 + M2)*(1-b+b*M2))
nu21<--1
A22<-4*b*M2/((2-M2)*(2+b*M2))
nu22<--M2/2
A23<-(1/b)/(1/b+M2)+b^2*M2^2/((2+b*M2)*(1-b+b*M2)*(1/b+M2))
nu23<--(1/b+M2)

A2<-c(A21,A22,A23)
nu2<-c(nu21,nu22,nu23)

A31<-M2/(M2-2)
nu31<--1
A32<-2/(2-M2)
nu32<--M2/2

A3<-c(A31,A32)
nu3<-c(nu31,nu32)

loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)

loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)

loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)

}


if (M2==0){

R1<-matrix(rep(r.a,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1a,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,3),nrow=3,byrow=TRUE)

loglike1<-sum(log(
(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
+
exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
   ))


loglike2<-sum(log(
(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
+
exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
   ))


loglike3<-sum(log(
exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
))

}


-(loglike1+loglike2+loglike3)
}




CI.prof.likelihood.3<-function(params1,i){
ci<-vector(length=2)
params<-params1[-i]
params[1]<-params[1]+0.25 		#to avoid "false convergence"

	 

##lower bound search


points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10	#the "round" function avoids, 
     									#for instance, truncating 2.3*100 									#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1-j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.T0,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.T0,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.T0,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[1]<-points[j-1]

##upper bound search:    

points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10
 									 									#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1+j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.T0,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.T0,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.T0,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[2]<-points[j-1]

ci
} 


CI.prof.likelihood.3(c(3.2731545,3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,16.7076085,
0.06416968),7)


###### For M2b/theta instead of M2

negll.IIM.uni8.analytic.1.M2bovertheta<-function(params,target,i){


if(i==1){
a<-target/params[1]
theta<-params[1]
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*theta/b
}

if(i==2){
a<-params[1]/target
theta<-target
b<-params[2]/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*theta/b
}

if(i==3){
a<-params[1]/params[2]
theta<-params[2]
b<-target/theta
c<-c(params[3]/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*theta/b
}

if(i==4){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(target/theta,params[4]/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*theta/b
}

if(i==5){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,target/theta)
tau1<-params[5]/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*theta/b
}

if(i==6){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-target/theta
tau0<-params[6]/theta+tau1
M2<-params[7]*theta/b
}

if(i==7){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-target/theta+tau1
M2<-params[7]*theta/b
}


if(i==8){
a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-params[7]/theta+tau1
M2<-target*theta/b
}



if (M2>0){

R1<-matrix(rep(r.a,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1a,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,2),nrow=2,byrow=TRUE)


A1<-1
nu1<--1

A21<-(b*M2^2)/((-2 + M2)*(1-b+b*M2))
nu21<--1
A22<-4*b*M2/((2-M2)*(2+b*M2))
nu22<--M2/2
A23<-(1/b)/(1/b+M2)+b^2*M2^2/((2+b*M2)*(1-b+b*M2)*(1/b+M2))
nu23<--(1/b+M2)

A2<-c(A21,A22,A23)
nu2<-c(nu21,nu22,nu23)

A31<-M2/(M2-2)
nu31<--1
A32<-2/(2-M2)
nu32<--M2/2

A3<-c(A31,A32)
nu3<-c(nu31,nu32)

loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)

loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)

loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)

}


if (M2==0){

R1<-matrix(rep(r.a,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r.bc,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1a,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3bc,3),nrow=3,byrow=TRUE)

loglike1<-sum(log(
(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
+
exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
   ))


loglike2<-sum(log(
(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
+
exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
   ))


loglike3<-sum(log(
exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
))

}


-(loglike1+loglike2+loglike3)
}


params1<-c(3.2731545,3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,16.7076085,
0.0382242)
i<-8
CI.prof.likelihood.4<-function(params1,i){
ci<-vector(length=2)
params<-params1[-i]
params[1]<-params[1]+0.25 		#to avoid "false convergence"

##lower bound search


points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10	
#the "round" function avoids, 
#for instance, truncating 2.3*100
#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1-j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.M2bovertheta,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.M2bovertheta,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1-j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
if(points[j]<0) break
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.M2bovertheta,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[1]<-points[j-1]

##upper bound search:    

points<-params1[i]
points1<-trunc(round(params1[i]*10,digits=1))/10
 									 									#to 229, instead of 230 
for (j in 1:100){
points<-c(points,points1+j*0.1)
}

w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.M2bovertheta,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*100,digits=2))/100
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.01)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.M2bovertheta,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}

points1<-trunc(round(points[j-1]*1000,digits=3))/1000
points<-points1
for (j in 1:100){
points<-c(points,points1+j*0.001)
}
w<-vector(length=0)
for (j in 1:length(points)){
x<-nlminb(params,
target=points[j],i=i,negll.IIM.uni8.analytic.1.M2bovertheta,lower=0)
if (x$objective==0|x$convergence!=0){break}
w<-c(w,2*(x$objective-89899.22))
if (w[j]>3.84){break}
}
ci[2]<-points[j-1]

ci
} 


CI.prof.likelihood.4(c(3.2731545,3.3575070,1.9288036,6.6231554,2.6468602,6.9300011,9.7776228,
0.0382242),8)

