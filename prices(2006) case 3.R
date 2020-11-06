#
# case 3
# 500 indivdiuals 
# 250 admix pop1, 250 admix pop2
# ancestry risk r=2
# to simulate case 4, set r=3.

# nsnps=50 
# nind=500
# n1=250 		## total number of population 1
# n2=nind-n1 	## total number of population 2
# R=1.5 		## relative risk of disease in cases vs. controls
# r=2 		## ancestry risk

case_3<-function(nsnps, nind, n1, R, r)
{

 nsnps=nsnps 
 nind=nind
 n1=n1
 R=R
 r=r 		
n2=nind-n1
maffilter=.90
fst=.01

###################################
## category 1 (Random SNPs)
###################################
p<-runif(nsnps,0.01,maffilter)
freq<-rbeta(nsnps,p*(1-fst)/fst,(1-p)*(1-fst)/fst)
snps1<-sapply(1:nsnps, function(x) sample(c(0,1,2),nind,replace=TRUE,prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))

###################################
## category 2 (Differentiated SNPs)
###################################
## ancestry allele frequency Pop1
p1<-rep(.80, nsnps)  
## ancestry allele frequency Pop2
p2<-rep(.20, nsnps)  
## Admix population 1
###################################
a=runif(1,0.1, 0.9)
a1=a
p=a*p1+(1-a)*p2
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))

pop1<-sapply(1:nsnps, function(x) sample(c(0,1,2),n1,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))
## Admix population2
###################################
a=runif(1,0.1, 0.9)
a2=a
p=a*p1+(1-a)*p2
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))

pop2<-sapply(1:nsnps, function(x) sample(c(0,1,2),n1,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))
snps2=rbind(pop1,pop2)
popID=c(rep(1,n1), rep(2,n2))

snps2=rbind(pop1,pop2)
popID=c(rep(1,n1), rep(2,n2))


###################################
## category 3 (Causal SNPs) 
################################### 
prop1=r^a1*0.5*log(r)/(r-1) ## disease probability of pop1
prop2=r^a2*0.5*log(r)/(r-1) ## disease probability of pop2
snps3 = matrix(NA, nind, nsnps)

###################################
## Sample cases/controls from POP1
###################################
nd1=round(prop1*n1) 		## number of cases from pop1
nc1=n1-nd1 				## number of controls from pop1

vec1=rep(0, n1)
vec1[sample(1:n1, size=nd1, replace = FALSE)]=1  ## cases/control status of pop1

a=a1									## allele frequencies
p=a*p1+(1-a)*p2
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))

pop1_cases<-sapply(1:nsnps, function(x) sample(c(0,1,2),nd1,replace=TRUE, prob=c(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(freq[x]^2*R^2)/sum((1-freq[x])^2,2*R*freq[x]*(1-freq[x]),R^2*freq[x]^2))))  
## sample cases
pop1_controls<-sapply(1:nsnps, function(x) sample(c(0,1,2),nc1,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))		 
## sample controls

snps3[which(vec1==1),]=pop1_cases
snps3[which(vec1==0),]=pop1_controls
###################################
## Sample cases/controls from POP2
###################################
nd2=round(prop2*n2) 				## number of cases from pop2
nc2=n2-nd2 				## number of controls from pop2

vec2=rep(0, n2)
vec2[sample(1:n2, size=nd2, replace = FALSE)]=1  
## cases/control status of pop1

a=a2
p=a*p1+(1-a)*p2
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))
pop2_cases<-sapply(1:nsnps, function(x) sample(c(0,1,2),nd2,replace=TRUE, prob=c(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(freq[x]^2*R^2)/sum(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(R^2*freq[x]^2))))) 	
## sample cases
pop2_controls<-sapply(1:nsnps, function(x) sample(c(0,1,2),nc2,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))			
## sample controls

snps3[n1+which(vec2==1),]=pop2_cases
snps3[n1+which(vec2==0),]=pop2_controls
disease=c(vec1, vec2)   ## disease status for pop1 and pop2

#################
colnames(snps1)=c(1:nsnps)
colnames(snps2)=c((nsnps+1):(2*nsnps))
colnames(snps3)=c((2*nsnps+1):(3*nsnps))
dt=cbind(snps1,snps2,snps3)

return(list(dt=dt, disease=disease, popID=popID));
}


##p=a*p1+(1-a)*p2					## allele frequencies
#freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))


