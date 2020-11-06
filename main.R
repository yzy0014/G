
source("~/MS code/binary_pop.R")
source("~/MS code/plot function.R")
source("~/MS code/power_cal.R")
source("~/MS code/test_case1.R")
source("~/MS code/test_case2.R")
source("~/MS code/test_case3.R")
source("~/MS code/prices(2006) case 1.R")
source("~/MS code/prices(2006) case 2.R")
source("~/MS code/prices(2006) case 3.R")


##################################################################################


## case 1 & 2
##################################################################################
#Discrete subpopulations with moderate ancestry differences bewteen Cases and Controls
####################################################################################
## snps: #of SNPs in each category(Random, Differentiated, Causal)
## nind:# of individuals in cases and controls
## n1: # of individuals in pop 1
## nd: # of diseased individuals in pop 1
## R: relative risk: the ratio= prob(diseased in case group)/Prob(diseased in control)
## a and 1-a:admix allele frequencies from pop 1 and pop2, respectively


m1=test_case1(num_loop<-400, relative_risk<-1.5, alpha<-.005/150, nsnps<-50, nind<-1000, n1<-500, nd<-500, a<-0.6)
m2=test_case2(num_loop<-200, relative_risk<-1.5, alpha<-.005/150, nsnps<-50, nind<-1000, n1<-250, nd<-500, a<-0.5)

save(m1, file="case 2(200).RData")
save(m2, file="case 2(200).RData")

#######################
## Table & plots 
## case 1
#######################
## generate Prices et al(2006) table 1
alp<-c(seq(0, .05, by=.01))
ls=prices_table(alp, m1)
ls
## generate ROC curve
alp<-c(seq(0,.00001, by=.00000001),seq(0.00001,.001, by=.00001),seq(0.001, 1, by=.001))
dt=power_typeI(alp, m1) ## obtain power and typeI-error data based on m3 (Prices case 3 simulation)
par(mfrow=c(1,2))
ROC_plot(dt$power, dt$typeI, title="case 1", xlim=c(0,1), ylim=c(0,1))
ROC_plot(dt$power, dt$typeI, title="case 1, typeI < .05", xlim=c(0,.05), ylim=c(0,0.5)) 


#######################
## Table & plots 
## case 2
#######################
## generate Prices et al(2006) table 1
alp<-c(seq(0, .05, by=.01))
ls=prices_table(alp, m2)
ls
## generate ROC curve
alp<-c(seq(0,.00001, by=.00000001),seq(0.00001,.001, by=.00001),seq(0.001, 1, by=.001))
dt=power_typeI(alp, m2) ## obtain power and typeI-error data based on m3 (Prices case 3 simulation)
par(mfrow=c(1,2))
ROC_plot(dt$power, dt$typeI, title="case2 B=200", xlim=c(0,1), ylim=c(0,1))
ROC_plot(dt$power, dt$typeI, title="case2 B=200, typeI < .05", xlim=c(0,.05), ylim=c(0,1)) 




## case 3 & 4 
##################################################################################
# admix population
####################################################################################
m3=test_case3(num_loop<-400, relative_risk<-1.5, alpha<-.05/150, nsnps<-50, nind<-1000, n1<-500, r<-2)
m4=test_case3(num_loop<-200, relative_risk<-1.5, alpha<-.05/150, nsnps<-50, nind<-1000, n1<-500, r<-3)
save(m3, file="case 3(400).RData")
save(m4, file="case 4(400).RData")


#######################
## Table & plots 
## case 3
#######################
## generate Prices et al(2006) table 1
alp<-c(seq(0, .05, by=.01))
ls=prices_table(alp, m3)
ls
## generate ROC curve
alp<-c(seq(0,.00001, by=.00000001),seq(0.00001,.001, by=.00001),seq(0.001, 1, by=.001))
dt=power_typeI(alp, m3) ## obtain power and typeI-error data based on m3 (Prices case 3 simulation)
par(mfrow=c(1,2))
ROC_plot(dt$power, dt$typeI, title="case3 B=400", xlim=c(0,1), ylim=c(0,1))
ROC_plot(dt$power, dt$typeI, title="case3 B=400, typeI < .05", xlim=c(0,.05), ylim=c(0,1)) 


#######################
## Table & plots 
## case 4
#######################
## generate Prices et al(2006) table 1
alp<-c(seq(0, .05, by=.01))
ls=prices_table(alp, m4)
ls
alp<-c(seq(0,.00001, by=.00000001),seq(0.00001,.001, by=.00001),seq(0.001, 1, by=.001))
dt2=power_typeI(alp, m4) 
par(mfrow=c(1,2))
ROC_plot(dt2$power, dt2$typeI, title="case4 r=3", xlim=c(0,1), ylim=c(0,1))
ROC_plot(dt2$power, dt2$typeI, title="case4 r=3, typeI < .05", xlim=c(0,.05), ylim=c(0,1))  

