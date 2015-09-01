#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope

Zhmm_PP<-function( run , burn=1000, prep=1000, isSim=TRUE,  simlabel="sim", MakePlots=TRUE){
Y<-run$YZ$Obs
Grun<-TrimThin(run, burn, Thin=1)

targetK0<-Grun$K0
K<-dim(run$Means)[2]
K0<-as.numeric(names(table(targetK0 )))
	n<-length(Y)
	#K0<-run$K0[burn:length(run$MAP)]

	K0estimates<-vector("list", length(K0))
	Zestimates<-vector("list", length(K0))
	# distribution of $K_0$
	emptyK<-run$K0[burn:length(run$MAP)]

	# parameters by stationary dist
    both<-cbind(melt(run$Means[-c(1:burn),])[,3], melt(run$q0[-c(1:burn),])[,3])
	color <- as.factor(melt(run$Means[-c(1:burn),])[,2])
	raincol<-rainbow(length(levels(color))); 	levels(color)<-raincol
	trancol<-sapply(c(1:length(color)), function(x) adjustcolor(color[x], alpha.f=both[x,2]))
	#minmaxMEANS<-c(min(both[both[,2]>minq4PLOT ,1]), max(both[both[,2]>minq4PLOT ,1]))

# PLOT 1
slices <- prop.table(table((factor(emptyK, levels=c(1:K)))))
Lab.palette <- colorRampPalette(rainbow(K*3, alpha=.3), space = "Lab")
if(MakePlots==TRUE){
pdf(file=paste(simlabel, "_MCMC.pdf",sep="") ,width=4, height=6, pointsize=8)
par(mfrow=c(2,1))
#pdf(file=paste(simlabel, "_K0.pdf",sep='') ,width=4, height=3, pointsize=8)
		barplot(slices, ylim=c(0,1), main="Distribution of number of occupied states", xlab="Number of occupied states", ylab="Proportion of iterations")
	abline(h=seq(0, 1, .05), lwd=0.5, col='LightGrey')
	#plot(both, col=trancol, xlim=minmaxMEANS, xlab="Mean", ylab="Stationary Distribution", bg='grey', main="Posterior Samples")
	#if(is.na(trueValues)==FALSE){	points(trueValues, pch=7, cex=2)}
#dev.off()


#pdf(file=paste(simlabel, "_SURF.pdf",sep='') ,width=4, height=4, pointsize=8)
	#smoothScatter(both, colramp = Lab.palette, main="", nrpoints = 0, xlab=expression(lambda), ylab=expression(mu[Q]))
	smoothScatter(both, colramp = Lab.palette, main="Posterior density", nrpoints = 0, ylab="Stationary distribution", xlab="State means")
	#plot(both, col=trancol, xlim=minmaxMEANS, xlab="Mean", ylab="Stationary Distribution", bg='grey', main="Posterior Samples")
	#if(is.na(trueValues)==FALSE){	points(trueValues, pch=7, cex=2)}
dev.off()
}

# unswitch seperate groups and make indi plots
	p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$K0)/dim(Grun$q0)[1]), "RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
# OFF BY ONE
#SteadyScore$K0
#FinalStates
		grunK0<-Grun
		# split data by K0
for ( .K0 in 1:length(K0)){
		.iterK0<- c(na.omit(c(1:dim(Grun$q0) [1])[targetK0==K0[.K0]]))
	if(p_vals$PropIters[.K0]>0.05){
		grunK0$Mu<-	Grun$Means[.iterK0,]
		grunK0$Q<-	Grun$Trans[.iterK0,]
		grunK0$q0<-	Grun$q0[.iterK0,]
		grunK0$MAP<-Grun$MAP[.iterK0]
		grunK0$Z<-	Grun$States[.iterK0,]
		grunK0$K0<-	Grun$K0[.iterK0]

		## 2. unswitch
		grunK0us<-Zswitch_hmm(grunK0, 0.05 )
		Zetc<-Zagg_HMM(grunK0us, Y)
		K0estimates[[.K0]]<-cbind(Zetc$thetaCI, "K0"=K0[.K0])
		Zestimates[[.K0]]<-Zetc$Zpred

		q0.mean = aggregate(q0 ~ k, grunK0us$Pars, mean)
		mu.mean = aggregate(mu ~ k, grunK0us$Pars, mean)

		p1<-ggplot(data=grunK0us$Pars, aes(y=q0, x=factor(k))) + geom_boxplot(fill='lightblue',outlier.size=0.5)+ylab("Stationary dist.")+xlab("State")  +  theme(legend.position = "none")+theme_bw()+			ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(round(p_vals$PropIters[.K0],2)), sep="")), atop("Stationary distribution"))))+ geom_text(data =q0.mean, aes(label=signif(q0,4)),size=4,  col='red',vjust = 1)
		p2<-ggplot(data=grunK0us$Pars, aes(y=mu, x=factor(k))) + geom_boxplot( fill='lightblue',outlier.size=0.5)+ylab("Mean")+xlab("State") +  theme(legend.position = "none")+theme_bw()+ggtitle( bquote(atop(italic(paste("Results for K=", .(K0[.K0]), sep="")), atop("Means"))))+geom_text(data =mu.mean, aes(label=signif(mu,4)),size=4,  col='red',vjust = 1)
		#Allocations
		p3<-HmmAllocationPlot(outZ=grunK0us$Z, myY=Y)
		## 4. Predict replicates
		postPredTests<-PostPredFunk(grunK0us,Zetc, Y, prep, simlabel)
		p4<-postPredTests$ggp
		# clusters:
		plotyz<-data.frame("X"=1:n, "Y"=Y, "Post_Z"=Zetc$Zpred)
		p5<-ggplot(plotyz, aes(y=Y,x=X ))+geom_point(aes(colour=Post_Z),size=3)+ geom_line(alpha = 1/4)+theme_bw()+  theme(legend.position = "none")+ggtitle("Posterior Allocations")+xlab("Time")
			#RAND
		if (isSim==TRUE){ p_vals$RAND[.K0]<- sum(run$YZ$States==Zetc$Zpred)/length(Zetc$Zpred)}
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-1-postPredTests$Concordance


		# print plot
		#pdf( file= paste(simlabel, "K_ ",K0[.K0] , ".pdf",sep='') ,width=10, height=5)
	 if(p_vals$PropIters[.K0]== max(p_vals$PropIters)){
	# if(p_vals$PropIters[.K0]>0.05){
if(MakePlots==TRUE){

	 	pdf( file= paste(simlabel, "K_",K0[.K0] , ".pdf",sep="") ,width=8, height=6, pointsize=8)
	 		print( wq::layOut(	
	 				list(p1, 	1, 1:2),
		        	list(p2, 	1, 3:4),
		         	list(p3,	1,5:6),
		         	list(p5, 	2,1:3),
		          	list(p4, 	2,4:6)))
		dev.off()
}
	}
		} # Close loop of >pmin
		} # Close loop over each K0






return(list(p_vals, K0estimates, Zestimates ))

}

