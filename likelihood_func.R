### LOG-LIKELIHOOD FUNCTIONS FOR MODELS "Isolation", "IM1", "IM2", 
### "IIM1" and "IIM2" -- (Figure 6 has the diagrams for these models)



#### "Isolation"

negll.iso<-function(params){

a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
tau<-params[4]/theta


R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r2,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

loglike1<-log(((theta*R1[1,])^X1[1,])/				
((1+theta*R1[1,])^(X1[1,]+1))*							
(1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
exp(-tau*(1+theta*R1[1,]))*
(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))

loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
((1+b*theta*R2[1,])^(X2[1,]+1))*							
(1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
exp(-tau*(1/b+theta*R2[1,]))*
(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))

loglike3<-log(exp(-tau*theta*R3[1,])*			
(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))

-sum(c(loglike1,loglike2,loglike3))
}



#####  "IM1"

negll.IM.1<-function(params){

a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
tau<-params[4]/theta
M<-vector(length=2,mode="numeric")
M[1]<-params[5]
M[2]<-params[5]


R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

if(params[5]>0){


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

loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
((theta*R1-nu)^(X1+1))*
(1-ppois(X1,tau*(-nu+theta*R1)))+
exp(-tau*(theta*R1-nu))*
(a*theta*R1)^X1/
(1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))

loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
((theta*R2-nu)^(X2+1))*
(1-ppois(X2,tau*(-nu+theta*R2)))+
exp(-tau*(theta*R2-nu))*
(a*theta*R2)^X2/
(1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))

loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
((theta*R3-nu)^(X3+1))*
(1-ppois(X3,tau*(-nu+theta*R3)))+
exp(-tau*(theta*R3-nu))*
(a*theta*R3)^X3/
(1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))

}

if(params[5]==0){

loglike1<-log(((theta*R1[1,])^X1[1,])/				
((1+theta*R1[1,])^(X1[1,]+1))*							
(1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
exp(-tau*(1+theta*R1[1,]))*
(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))

loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
((1+b*theta*R2[1,])^(X2[1,]+1))*							
(1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
exp(-tau*(1/b+theta*R2[1,]))*
(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))

loglike3<-log(exp(-tau*theta*R3[1,])*			
(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))

}


-sum(c(loglike1,loglike2,loglike3))
}




#### "IM2"

negll.IM.2<-function(params){

a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
tau<-params[4]/theta
M<-vector(length=2,mode="numeric")
M[1]<-params[5]
M[2]<-params[6]




if (M[1]>0 & M[2]>0){

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

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


loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
((theta*R1-nu)^(X1+1))*
(1-ppois(X1,tau*(-nu+theta*R1)))+
exp(-tau*(theta*R1-nu))*
(a*theta*R1)^X1/
(1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))

loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
((theta*R2-nu)^(X2+1))*
(1-ppois(X2,tau*(-nu+theta*R2)))+
exp(-tau*(theta*R2-nu))*
(a*theta*R2)^X2/
(1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))

loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
((theta*R3-nu)^(X3+1))*
(1-ppois(X3,tau*(-nu+theta*R3)))+
exp(-tau*(theta*R3-nu))*
(a*theta*R3)^X3/
(1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))

}

if (M[1]==0 & M[2]>0){

R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)


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


loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
(-nu1+theta*R1)^(X1+1)*
(1-ppois(X1,tau*(-nu1+theta*R1)))+
exp(-tau*(theta*R1-nu1))*
(a*theta*R1)^X1/
(1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))

loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
(-nu2+theta*R2)^(X2+1)*
(1-ppois(X2,tau*(-nu2+theta*R2)))+
exp(-tau*(theta*R2-nu2))*
(a*theta*R2)^X2/
(1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))

loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
(-nu3+theta*R3)^(X3+1)*
(1-ppois(X3,tau*(-nu3+theta*R3)))+
exp(-tau*(theta*R3-nu3))*
(a*theta*R3)^X3/
(1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))

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

loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
(-nu1+theta*R1)^(X1+1)*
(1-ppois(X1,tau*(-nu1+theta*R1)))+
exp(-tau*(theta*R1-nu1))*
(a*theta*R1)^X1/
(1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))

loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
(-nu2+theta*R2)^(X2+1)*
(1-ppois(X2,tau*(-nu2+theta*R2)))+
exp(-tau*(theta*R2-nu2))*
(a*theta*R2)^X2/
(1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))

loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
(-nu3+theta*R3)^(X3+1)*
(1-ppois(X3,tau*(-nu3+theta*R3)))+
exp(-tau*(theta*R3-nu3))*
(a*theta*R3)^X3/
(1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))


}

if (M[1]==0 & M[2]==0){

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

loglike1<-log(((theta*R1[1,])^X1[1,])/				
((1+theta*R1[1,])^(X1[1,]+1))*							
(1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
exp(-tau*(1+theta*R1[1,]))*
(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))

loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
((1+b*theta*R2[1,])^(X2[1,]+1))*							
(1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
exp(-tau*(1/b+theta*R2[1,]))*
(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))

loglike3<-log(exp(-tau*theta*R3[1,])*			
(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))

}

-sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}




##### "IIM1"

negll.IIM.1<-function(params){

a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(1,params[3]/theta)
tau1<-params[4]/theta
tau0<-params[5]/theta+tau1
M<-vector(length=2,mode="numeric")
M[1]<-params[6]
M[2]<-params[7]


if (M[1]>0 & M[2]>0){

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)


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

R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)


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

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)

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

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

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



####### "IIM2"

negll.IIM.2<-function(params){

a<-params[1]/params[2]
theta<-params[2]
b<-params[3]/theta
c<-c(params[4]/theta,params[5]/theta)
tau1<-params[6]/theta
tau0<-params[7]/theta+tau1
M<-vector(length=2,mode="numeric")
M[1]<-params[8]
M[2]<-params[9]

if (M[1]>0 & M[2]>0){

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)


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

R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)


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

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)

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

R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)

X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

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



