#' Processing output of label switching function
#'
#' explain here
#' @param output of HMM label switching function, computes parameter estimates and some error values. 
#' @keywords postprocessing
#' @export
#' @examples
#' #nope

Zagg_HMM<-function(USout, .Y = Y) {
.par <- melt(USout$Pars, id.vars = c("Iteration", "k"))
theta <- aggregate(value ~  factor(k)+variable , mean, data = .par)

		theta<-aggregate( value~factor(k)+variable, mean ,data=.par)
		mu<-round(aggregate( value~factor(k)+variable, mean ,data=.par)[,3], 2)
		ci<-round(aggregate( value~factor(k)+variable, quantile,c(0.025, 0.975) ,data=.par)[,3],2)
	#	thetaCI<-cbind( theta[,c(1,2)] , "value"=paste( mu, "(", ci[,1] , "," ,ci[,2] ,")", sep=" " ))
		thetaCI<-cbind( theta[,c(1,2)] , "Mean"= mu, "CI.025"=ci[,1] , "CI.975"=ci[,2] )
names(thetaCI)[1:2]<-c("k","Parameter")
K <- max(.par$k)
Zhat<- factor( apply(USout$Z, 2,maxZ))#[-((length(.Y))+1)]
Zemu <- as.numeric(Zhat)
.Mus <- theta$value[theta$variable == "mu"]
for (i in 1:length(Zemu)) { Zemu[i] <- .Mus[as.numeric(Zhat[i])]  }
MSE <- sum((.Y - Zemu)^2)
MAE <- sum(abs(.Y - Zemu))
list(theta = theta, thetaCI=thetaCI, Zpred = Zhat, MSE = MSE, MAE = MAE)
				}