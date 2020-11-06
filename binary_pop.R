######################################################################
#  translating SNP dataset into Binary
#
######################################################################
binary_pop<-function(pop){
n=nrow(pop)
p=ncol(pop)
popb=matrix(0,n,2*p)
popb[,2*1:p-1]=pop
	for(i in 1:n)
	{
		rank=which(popb[i,]==2)
		if (length(rank)>0){
			for(j in 1:length(rank)){
				point=rank[j]
				popb[i, point:(point+1)]=c(1,1)
			}
		}
	}
	return(popb);
}
