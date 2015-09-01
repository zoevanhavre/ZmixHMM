 #' ReplicateSimer2
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

ReplicateSimer2<-function(  N, n, Kfit=10, SimID, ITERATIONS, BURN,  AMAX, AMIN, PRIOR_TYPE, PTchain=20){
		#  STORE SIMULATIONS in a list
		simFunctionMorpher<-function(SimNumber){
			if(	SimNumber==1){ 	return(FunkSim1)
			}else if (SimNumber==2){	return(FunkSim2)
			}else if (SimNumber==3){	return(FunkSim4)
			}	}
		MorphingSIMULATE<-simFunctionMorpher(SimID)
		SIMS<-lapply( rep(n,N),  MorphingSIMULATE )

		# Clean up Gibbs for lyra...

Result.store<-data.frame("Replicate"=c(1:N), "SimID"=SimID, "n"=n,"AlphaMax"=AMAX, "Prior"=PRIOR_TYPE, "ModeK0"=0, "MeanfDist"=0, "MeanfDistMERGED"=0, "WorstMixed"=0)

for (.rep in 1:N){
My.Result<-gibbsHMM_Main(YZ=SIMS[[.rep]],K=Kfit, ,  M=ITERATIONS,  alphaMAX=AMAX, type= PRIOR_TYPE, alphaMin=AMIN, J=PTchain, SuppressAll="TRUE")

Result.store$ModeK0[.rep]<-as.numeric(names(sort(table(factor(My.Result$Track$K0[-c(1:BURN)])),decreasing=TRUE)[1]))
Result.store$MeanfDist[.rep]<-mean(My.Result$Track$f2Dist[-c(1:BURN)])
Result.store$MeanfDistMERGED[.rep]<-mean(My.Result$Track$f2Dist_Merged[-c(1:BURN)])
Result.store$WorstMixed[.rep]<-min(My.Result$Track$WorstMixProp[-c(1:BURN)])

write.csv(Result.store[1:.rep,], file=paste( "RepResult_Sim" ,SimID,"n",n, "Prior", PRIOR_TYPE, "Alpha", AMAX,"Iters",ITERATIONS,".csv", sep=""))
save(Result.store, file=paste( "RepResult_Sim" ,SimID,"n",n, "Prior", PRIOR_TYPE,"Alpha",round( AMAX,3) ,"Iters",ITERATIONS, ".RDATA", sep="_"))

Sys.sleep(0.1)
print(Result.store[1:.rep,])
Sys.sleep(0.1)
}

pdf( file= paste( "Sim" ,SimID,"n",n, "Prior", PRIOR_TYPE, "MaxAlpha", AMAX,"Iters",ITERATIONS, ".pdf", sep="") ,width=5, height=3)
	 		print( wq::layOut(
	 		list(ggplot(data=Result.store, aes(x=ModeK0))+geom_histogram(binwidth=1)+theme_bw()+xlab("K_0")+
	 			ggtitle(paste( "Sim" ,SimID ,"(n=" ,n, ") Prior", PRIOR_TYPE,  "Alpha", round( AMAX,3) , ": K_0"))  ,	1 ,	 1:2),
		    list(ggplot(data=Result.store, aes(y=MeanfDist, x= factor(1)))+geom_boxplot()+theme_bw()+ylab("Distance") ,	1,	 3)  )    )

dev.off()

	return(Result.store)
		}



