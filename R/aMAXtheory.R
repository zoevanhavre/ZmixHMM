
 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

aMAXtheory<-function( K ,  Amin=0.0000001){
	p=1;d=1
	r1<-((K-1)*(K-2+d)+Amin*K)*((K-1)*(d+K-2)+(d/2))
    	r2<- r1/((d/2)-Amin*(K^2-3*K+4))
	return((r2-(K-p)*Amin)/p)
   	}

#  aMAXtheory<-function(amin, K, ref=c("Rousseau", "Gassiat")){
# 	if(ref=="Gassiat"){ simonSays<-K*((K-1)-1+1)-(K-1)*amin+0.01}
# 	if(ref=="Rousseau"){ simonSays<-(K-1)*(1+K-2+amin)*(1+(1/(.5-amin*(K-1)) ))-amin*(K-1)+0.01}

# 	return(simonSays)

# }
