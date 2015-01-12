### MAXIMUM-LIKELIHOOD POINT ESTIMATION AND WALD CONFIDENCE INTERVALS 
### FOR MODELS "Isolation", "IM1", "IM2", "IIM1" and "IIM2" 
### (Figure 6 has the diagrams for these models)

### Running the code below will yield, for each model:
### 1) The list "par.[model]" of maximum-likelihood estimates;
### 2) The list "se.[model]" of Wald standard errors.
### NOTE: the estimates are meaningful if (and only if) 'nlminb'
### returns "relative convergence" or "absolute convergence". If not,
### other initial values should be used.

### Before you run the code below, you need:
### -- to install and load the packages 'numDeriv' (contains the function 
### 'hessian') and 'stats' (contains the function 'nlminb').
### -- to run the file 'likelihood_func', containing the likelihood functions;
### -- to create vectors 'x1', 'x2' and 'x3' containing the numbers of
### different nucleotides between pairs of sequences from different loci ('x1'
### if both sequences come from population 1, 'x2' if both sequences come
### from population 2 and 'x3' if there is one sequence from each population);
### -- to create vectors 'r1', 'r2' and 'r3' with the relative mutation rates 
### of the loci in 'x1', 'x2' and 'x3' respectively.


##### "Isolation"

par.iso<-nlminb(rep(1.5,4),negll.iso,lower=0.000001)$par
par.iso<-list(A=par.iso[1],theta=par.iso[2],B=par.iso[3],V=par.iso[4])

hessian.iso<-hessian(negll.iso,par.iso)
se.iso<-sqrt(diag(solve(hessian.iso)))
se.iso<-list(A=se.iso[1],theta=se.iso[2],B=se.iso[3],V=se.iso[4])

##### "IM1"

par.IM.1<-nlminb(rep(1.5,5),negll.IM.1,
lower=c(0.000001,0.000001,0.000001,0.000001,0))$par
par.IM.1<-list(A=par.IM.1[1],theta=par.IM.1[2],B=par.IM.1[3],V=par.IM.1[4],
M=par.IM.1[5])

hessian.IM.1<-hessian(negll.IM.1,par.IM.1)
se.IM.1<-sqrt(diag(solve(hessian.IM.1)))
se.IM.1<-list(A=se.IM.1[1],theta=se.IM.1[2],B=se.IM.1[3],
V=se.IM.1[4],M=par.IM.1[5])

##### "IM2"

par.IM.2<-nlminb(rep(1.4,6),negll.IM.2,
lower=c(0.000001,0.000001,0.000001,0.000001,0,0))$par
par.IM.2<-list(A=par.IM.2[1],theta=par.IM.2[2],B=par.IM.2[3],
V=par.IM.2[4],M1=par.IM.2[5],M2=par.IM.2[6])

hessian.IM.2<-hessian(negll.IM.2,par.IM.2)
se.IM.2<-sqrt(diag(solve(hessian.IM.2)))
se.IM.2<-list(A=se.IM.2[1],theta=se.IM.2[2],B=se.IM.2[3],
V=se.IM.2[4],M1=se.IM.2[5],M2=se.IM.2[6])


##### "IIM1"

par.IIM.1<-nlminb(rep(1.4,7),negll.IIM.1,
lower=c(0.000001,0.000001,0.000001,0.000001,0.000001,0,0))$par
par.IIM.1<-list(A=par.IIM.1[1],theta=par.IIM.1[2],B=par.IIM.1[3],
T1=par.IIM.1[4],V=par.IIM.1[5],M1=par.IIM.1[6],M2=par.IIM.1[7])

hessian.IIM.1<-hessian(negll.IIM.1,par.IIM.1)
se.IIM.1<-sqrt(diag(solve(hessian.IIM.1)))
se.IIM.1<-list(A=se.IIM.1[1],theta=se.IIM.1[2],B=se.IIM.1[3],
T1=se.IIM.1[4],V=se.IIM.1[5],M1=se.IIM.1[6],M2=se.IIM.1[7])



##### "IIM2"

par.IIM.2<-nlminb(rep(1.4,9,negll.IIM.2,lower=c(0.000001,0.000001,
0.000001,0.000001,0.000001,0.000001,0.000001,0,0))$par
par.IIM.2<-list(A=par.IIM.2[1],theta=par.IIM.2[2],B=par.IIM.2[3],
C1=par.IIM.2[4],C2=par.IIM.2[5],T1=par.IIM.2[6],V=par.IIM.2[7],
M1=par.IIM.2[8],M2=par.IIM.2[9])

hessian.IIM.2<-hessian(negll.IIM.2,par.IIM.2)
se.IIM.2<-sqrt(diag(solve(hessian.IIM.2)))
se.IIM.2<-list(A=se.IIM.2[1],theta=se.IIM.2[2],B=se.IIM.2[3],
C1=se.IIM.2[4],C2=se.IIM.2[5],T1=se.IIM.2[6],V=se.IIM.2[7],
M1=se.IIM.2[8],M2=se.IIM.2[9])
