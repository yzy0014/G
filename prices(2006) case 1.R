#######################################################################
# case 1
######################################################################
# 250 cases, 250 controls
# 60% cases from pop1; rest cases and all controls from pop2
# 250 pop1, 250 pop2
# nsnps=50 
# nind=500
# n1=250 		## total number of population 1
# n2=nind-n1 	## total number of population 2
# nd=250  		## total number of cases
# nc=nind-nd  ## total number of controls
# R=1.5 		## relative risk of disease in cases vs. controls
# a=0.6 		## proportion of cases from pop1


case_1<-function(nsnps, nind, n1, nd, R, a){
 nsnps=nsnps 
 nind=nind
 n1=n1
 R=R
 a=a 
n2=nind-n1
nc=nind-nd
maffilter=.90
fst=.01
###################################
## category 1 (Random SNPs)
###################################
p<-runif(nsnps,0.01,maffilter)
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))
snps1<-sapply(1:nsnps, function(x) sample(c(0,1,2),nind,replace=TRUE,prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))

###################################
## category 2 (Differentiated SNPs)
###################################
## Population 1, allel frequency=.80
p<-rep(.80, nsnps)
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))
pop1<-sapply(1:nsnps, function(x) sample(c(0,1,2),n1,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))

## Population 2, allel frequency=.20
p<-rep(.20, nsnps)
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))
pop2<-sapply(1:nsnps, function(x) sample(c(0,1,2),n2,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))
snps2=rbind(pop1,pop2)
popID=c(rep(1,n1), rep(2,n2))

###################################
## category 3 (Causal SNPs)  
## Sample cases/controls from POP1
###################################
snps3 = matrix(NA, nind, nsnps)
a=.60
nd1=round(a*nd) 		## number of cases from pop1
nc1=n1-nd1 				## number of controls from pop1

vec1=rep(0, n1)
vec1[sample(1:n1, size=nd1, replace = FALSE)]=1  ## cases/control status of pop1

p<-rep(.80, nsnps)					## allele frequencies
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))

pop1_cases<-sapply(1:nsnps, function(x) sample(c(0,1,2),nd1,replace=TRUE, prob=c(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(freq[x]^2*R^2)/sum(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(R^2*freq[x]^2)))))  ## sample cases
pop1_controls<-sapply(1:nsnps, function(x) sample(c(0,1,2),nc1,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))		 ## sample controls

snps3[which(vec1==1),]=pop1_cases
snps3[which(vec1==0),]=pop1_controls
###################################
## Sample cases/controls from POP2
###################################
nd2=nd-nd1 				## number of cases from pop2
nc2=n2-nd2 				## number of controls from pop2

vec2=rep(0, n2)
vec2[sample(1:n2, size=nd2, replace = FALSE)]=1  ## cases/control status of pop2

p<-rep(.20, nsnps)					## allele frequencies
freq<-sapply(1:length(p),function(x) rbeta(1,p[x]*(1-fst)/fst,(1-p[x])*(1-fst)/fst))

pop2_cases<-sapply(1:nsnps, function(x) sample(c(0,1,2),nd2,replace=TRUE, prob=c(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(freq[x]^2*R^2)/sum(((1-freq[x])^2),(2*R*freq[x]*(1-freq[x])),(R^2*freq[x]^2))))) 	## sample cases
pop2_controls<-sapply(1:nsnps, function(x) sample(c(0,1,2),nc2,replace=TRUE, prob=c(((1-freq[x])^2),(2*freq[x]*(1-freq[x])),(freq[x]^2))))				## sample controls


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








