######################################################################
library("coin")
library("mvtnorm")
library("logisticPCA")
library("psych")
library("mlogit")
## source("prices(2006) case 1.R")
######################################################################

test_case2 <- function(num_loop, relative_risk, alpha, nsnps, nind, n1, nd, a) {

	fst = 0.01
	maffilter = 0.9

	############################################
	## simulation parameters 
############################################

	B = num_loop
	R = relative_risk
	alp = alpha
	nsnps = nsnps
	nind = nind
	n1=n1
	nd=nd
	a=a

	chi_null = matrix(NA, nrow = B, ncol = 3 * nsnps)
	chi_poly = matrix(NA, nrow = B, ncol = 3 * nsnps)
	chi_lpca = matrix(NA, nrow = B, ncol = 3 * nsnps)
	chi_pca = matrix(NA, nrow = B, ncol = 3 * nsnps)

	for (i in 1:B) {
		print(i)

		###########################
		##
## requires "princes (2006) case_1.R" to run.
##
##
###########################
	
dataset <- case_2(nsnps, nind, n1, nd, R, a)
pop = dataset$dt
ls=apply(pop, 2, var, na.rm=TRUE) != 0
flag=length(which(ls==FALSE))
while(flag){
	dataset <- case_2(nsnps, nind, n1, nd, R, a)
	pop = dataset$dt
	ls=apply(pop, 2, var, na.rm=TRUE) != 0
	flag=length(which(ls==FALSE))
}		
		l = dataset$disease
		popID=dataset$popID

		#########################
		## logistic PCA
#########################    
lpca = logisticPCA(binary_pop(pop), k = 1, m = 2, main_effects = TRUE)
		x.lpca = lpca$PCs

		#########################
		# polychoric PCA
#########################
Poly_cor <- polychoric(pop)$rho
		out <- principal(Poly_cor, nfactors = 1, scores = TRUE)
		x.poly <- factor.scores(pop, out)$scores

		#########################
		# continous PCA
#########################

		cpca <- principal(cor(pop), nfactor = 1, scores = TRUE)
		x.pca <- factor.scores(pop, cpca)$scores



		###################################
		# initialize testing
###################################
		 l.pca = glm(l ~ x.pca, family = binomial("logit"))$residuals
		 l.poly = glm(l ~ x.poly, family = binomial("logit"))$residuals
		 l.lpca = glm(l ~ x.lpca, family = binomial("logit"))$residuals

#l.pca = lm(l ~ x.pca)$residuals
#l.poly = lm(l ~ x.poly)$residuals
#l.lpca = lm(l ~ x.pca)$residuals

		for (j in 1:ncol(pop)) {
			
###########################################
## Armitage trend test
##
###########################################
## print(j)

			SNP = pop[, j]
			data=data.frame(cbind(SNP, x.pca,x.poly, x.lpca))
			row.names(data)=c(1:nind)
			colnames(data)=c("SNP","x.pca","x.poly","x.lpca")
	
			SNP.pca = mlogit(SNP ~ 0 | x.pca, data=data, shape="wide", family="multinomial")$residuals
			SNP.lpca = mlogit(SNP ~ 0 | x.lpca, data=data, shape="wide", family="multinomial")$residuals
			SNP.poly = mlogit(SNP ~ 0 | x.poly, data=data, shape="wide", family="multinomial")$residuals
			
			chi_null[i, j] = pvalue(independence_test(l ~ SNP, teststat = "quad"))
			chi_pca[i, j] = pvalue(independence_test(l.pca ~ SNP.pca, teststat = "quad"))
			chi_lpca[i, j] = pvalue(independence_test(l.lpca ~ SNP.lpca, teststat = "quad"))
			chi_poly[i, j] = pvalue(independence_test(l.poly ~ SNP.poly, teststat = "quad"))

		}

	}

	return(list(pNA=chi_null, pPCA=chi_pca, pLPCA=chi_lpca, pPOLY=chi_poly));
}



