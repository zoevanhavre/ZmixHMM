 #' gibbsHMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

trimit_HMM<-function(Out=Out, nEnd=EndSize){
				nmax<-length(Out$MAP)
				mu<-Out$Means[c(nmax-nEnd+1):nmax,]    
				sig<-Out$Trans[c(nmax-nEnd+1):nmax,]
				ps<-Out$q0[c(nmax-nEnd+1):nmax,]
				Loglike<-Out$MAP[c(nmax-nEnd+1):nmax]
				zs<-Out$States[c(nmax-nEnd+1):nmax,]
				K0<-Out$K0[c(nmax-nEnd+1):nmax ,]
	list("Mu"= mu,"Q"=sig, "q0" = ps, "MAP"=Loglike, "K0"=K0, Z=zs, YZ=Out$YZ)	
						}