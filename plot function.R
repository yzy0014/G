## prices_table: 
## generate table 1 from Prices et al. (2006) paper.
## Ex:
## alp<-c(seq(0, 1, by=.001))
## m3=test_case3(num_loop<-100, relative_risk<-1.5, alpha<-.05/150, nsnps<-50, nind<-1000, n1<-500, r<-2)
## ls=prices_table(alp, m3)
## ls

prices_table<-function(alp, model){
	alp=alp;
	proportion=matrix(NA, nrow=4, ncol=3)
for(i in 1: length(alp)){
	proportion = power_cal(model, alp[i])
   if(i==1) {
 			ls=list(proportion)
   	}else if(i>1){
   		   	ls=c(ls, list(proportion))
   	} 	 		
}
return(ls)
}

## power_typeI
## to generate power data and typeI error data based on p-value results
## alp: a vector of significant levels; 
## model: an object of test_case function.
## EX:
## alp<-c(seq(0, 1, by=.001))
## m3=test_case3(num_loop<-100, relative_risk<-1.5, alpha<-.05/150, nsnps<-50, nind<-1000, n1<-500, r<-2)
## dt=power_typeI(alp, m3)

power_typeI<-function(alp, model){
	alp=alp;
	model=model;
power=matrix(NA, nrow=4, ncol=length(alp))
typeI=matrix(NA, nrow=4, ncol=length(alp))
for(i in 1: length(alp)){
	proportion = power_cal(model, alp[i])
	power[,i] = proportion[,3]
 	typeI[,i] = rowSums(proportion[,1:2])/2
	}
	rownames(power)=c("null", "pca", "lpca", "poly")
	rownames(typeI)=c("null", "pca", "lpca", "poly")
	return(list(power=power, typeI=typeI));
	}


## ROC_plot
## to generate ROC curve based on: power data, typeI error data, title, and size of the graph.
## EX:
## alp<-c(seq(0, 1, by=.001))
## m3=test_case3(num_loop<-100, relative_risk<-1.5, alpha<-.05/150, nsnps<-50, nind<-1000, n1<-500, r<-2)
## dt=power_typeI(alp, m3)
## ROC_plot(dt, title="case3 B=100", xlim=c(0,1), ylim=c(0,1)) ## to generate full ROC curve
## ROC_plot(dt, title="case3 B=100, typeI < .05", xlim=c(0,.05), ylim=c(0,0.5)) ## to generate partial ROC curve

ROC_plot<-function(power, typeI,title, xlim, ylim){
	title=title;
	xlim=xlim;
	ylim=ylim;
	plot( typeI[1,], power[1,], type='b', pch=1, cex=.5, col=1,xlim=c(xlim), ylim=c(ylim),main=title, xlab="False positive-rate", ylab="True positive-rate")  ## black; null model
	lines(typeI[2,], power[2,],type='b', pch=2, cex=.5,col=2) ## red; pca
	lines(typeI[3,], power[3,],type='b', pch=3, cex=.5, col=3) ## green;logistic pca 
	lines(typeI[4,], power[4,],type='b', pch=4, cex=.5,col=4) ## blue; polychoric pca
	legend('bottomright', c("null", "pca", "lpca", "poly") , 
   lty=1, col=c(1:4), bty='n', cex=.75) 
}