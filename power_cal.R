power_cal<-function(model, alp){
	
	###########################################

	alp = alp 
	
	chi_null=model$pNA
	chi_pca=model$pPCA
	chi_lpca=model$pLPCA
	chi_poly=model$pPOLY

	B=nrow(chi_null)
	
	#########################################################	

	sig_chi_null = matrix(0, nrow = B, ncol = 3 * nsnps)
	sig_chi_poly = matrix(0, nrow = B, ncol = 3 * nsnps)
	sig_chi_lpca = matrix(0, nrow = B, ncol = 3 * nsnps)
	sig_chi_pca = matrix(0, nrow = B, ncol = 3 * nsnps)

	sig_chi_null[which(chi_null < alp)] = 1
	sig_chi_poly[which(chi_poly < alp)] = 1
	sig_chi_lpca[which(chi_lpca < alp)] = 1
	sig_chi_pca[which(chi_pca < alp)] = 1

	cat_log = matrix(NA, 4, 3)
	cat_chi = matrix(NA, 4, 3)
	
	if (B==1){
		###########################################
		## first col: snps of category1; 
		## second col: snps of category2; 
		## thrid col: snps of category3
		###########################################
   cat_chi[1, ] = c(sum( (sig_chi_null[1:50]))/50, sum( (sig_chi_null[51:100]))/50, sum( (sig_chi_null[
		101:150]))/50) ## first row: null
	cat_chi[2, ] = c(sum( (sig_chi_pca[1:50]))/50, sum( (sig_chi_pca[51:100]))/50, sum( (sig_chi_pca[
		101:150]))/50) ## second row: pca
	cat_chi[3, ] = c(sum( (sig_chi_lpca[1:50]))/50, sum( (sig_chi_lpca[51:100]))/50, sum( (sig_chi_lpca[ 
		101:150]))/50) ## third row: lpca
	cat_chi[4, ] = c(sum( (sig_chi_poly[1:50]))/50, sum( (sig_chi_poly[51:100]))/50, sum( (sig_chi_poly[ 
		101:150]))/50) ## forth row: polychoric pca

	} else if(B>1){
		###########################################
		## first col: snps of category1; 
		## second col: snps of category2; 
		## thrid col: snps of category3
		###########################################
cat_chi[1, ] = c(sum(colMeans(sig_chi_null[, 1:50]))/50, sum(colMeans(sig_chi_null[, 51:100]))/50, sum(colMeans(sig_chi_null[, 
		101:150]))/50) ## first row: null
	cat_chi[2, ] = c(sum(colMeans(sig_chi_pca[, 1:50]))/50, sum(colMeans(sig_chi_pca[, 51:100]))/50, sum(colMeans(sig_chi_pca[, 
		101:150]))/50) ## second row: pca
	cat_chi[3, ] = c(sum(colMeans(sig_chi_lpca[, 1:50]))/50, sum(colMeans(sig_chi_lpca[, 51:100]))/50, sum(colMeans(sig_chi_lpca[, 
		101:150]))/50) ## third row: lpca
	cat_chi[4, ] = c(sum(colMeans(sig_chi_poly[, 1:50]))/50, sum(colMeans(sig_chi_poly[, 51:100]))/50, sum(colMeans(sig_chi_poly[, 
		101:150]))/50) ## forth row: polychoric pca

	}

	colnames(cat_chi) = c("random SNPs", "differentiated SNPs", "casual SNPs")
	rownames(cat_chi) = c("null", "pca", "lpca", "poly")

	return(cat_chi);
}




	#sig_log_null=matrix(0, nrow=B, ncol=3*nsnps)
#sig_log_poly=matrix(0, nrow=B, ncol=3*nsnps)
#sig_log_lpca=matrix(0, nrow=B, ncol=3*nsnps)
#sig_log_pca=matrix(0, nrow=B, ncol=3*nsnps)
	## with bonferroni adjust; .05/150
#sig_log_null[which(log_null<alp)]=1
#sig_log_poly[which(log_poly<alp)]=1
#sig_log_lpca[which(log_lpca<alp)]=1
#sig_log_pca[which(log_pca<alp)]=1
	#cat_log = matrix(NA, 4, 3)
	###########################################
	## first col: snps of category1; 
## second col: snps of category2; 
## thrid col: snps of category3
###########################################
#cat_log[1,]=c(sum(colMeans(sig_log_null[,1:50]))/50, sum(colMeans(sig_log_null[,51:100]))/50,sum(colMeans(sig_log_null[,101:150]))/50)  ## first row: null
#cat_log[2,]=c(sum(colMeans(sig_log_pca[,1:50]))/50, sum(colMeans(sig_log_pca[,51:100]))/50,sum(colMeans(sig_log_pca[,101:150]))/50)   ## second row: pca
#cat_log[3,]=c(sum(colMeans(sig_log_lpca[,1:50]))/50, sum(colMeans(sig_log_lpca[,51:100]))/50,sum(colMeans(sig_log_lpca[,101:150]))/50)  ## third row: lpca
#cat_log[4,]=c(sum(colMeans(sig_log_poly[,1:50]))/50, sum(colMeans(sig_log_poly[,51:100]))/50,sum(colMeans(sig_log_poly[,101:150]))/50)  ## forth row: polychoric pca

#colnames(cat_log) = c("random SNPs", "differentiated SNPs", "casual SNPs")
#rownames(cat_log) = c("null", "pca", "lpca", "poly")
#list(cat_log = cat_log, cat_chi = cat_chi)

