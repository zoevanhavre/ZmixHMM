 #' gibbsHMM_PT
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


gibbsHMM_PT_wDist_LYRAfinally<-function(YZ, M=2000, K=10 ,alphaMin=0.00001, alphaMAX=1,  type= 1, J=1,  SuppressAll=FALSE, TV=NULL, gridN=100 ){
  #____SET UP_________________________________________
  ifelse(class(YZ)=='data.frame',    Y<-YZ$Observed, Y<-YZ)
  n=length(Y) # sample size
  varknown<-1 # known variace
  #   mu0=0
  # make mu0 mid range of Y like in mixtures
  mu0=mean(Y)
  var0=100
  # INITIALIZE
  states0<-lapply(rep(K, J), function(x) makeStartSimpler(Y, x))
  TrackParallelTemp<-data.frame("Chain"=c(1:J),"NumTries"=rep(0,J), "NumSuccess"=rep(0,J))

  # final
  FinalMu<- matrix(nrow=M, ncol=K)
  FinalQ<- matrix(nrow=M, ncol=K*K)
  Finalq0<- matrix(nrow=M, ncol=K)
  FinalStates<- matrix(nrow=M, ncol=n+1)
  SteadyScore<-data.frame("Iteration"=c(1:M), "K0"=K)  #FIX

  #temp storage
  MU<- matrix(nrow=J, ncol=K)
  Q <- matrix(nrow=J, ncol=K*K)
  Qold<- matrix(nrow=J, ncol=K*K)
  q0 <- matrix(nrow=J, ncol=K)
  Z  <- matrix(nrow=J, ncol=n+1) #include 0 for initial state to be estimated too?

  PTsuccess  <- data.frame("Chain"=c(1:J),"Tries"=0,"Success"=0, "Ratio"=NA) #include 0 for initial state to be estimated too?
  MAP<-rep(0,M)  # KEEP TARGET ONLY
  f2now<-vector(length=M)
  fDist<-vector(length=M)
  f2now_MERGED<-vector(length=M)
  fDistMERGED<-vector(length=M)

      # ALPHA
    # alphaMin<-max( 1/(100*n), 0.000001)
    # if( Theory==FALSE){alphaMAX=1
    # } else { alphaMAX<-AlphaMax_Theory(K, alphaMin) }

  Alpha_lows<-c(alphaMAX, exp(seq(log(alphaMAX), log(alphaMin), length=J))[-1])
  #Store alphs for PT
  AllAlphas<- lapply(c(1:J), function(x) matrix(Alpha_lows[x],ncol=K, nrow=K)  )

  AllAlphas<-lapply(AllAlphas, function(x) {
  if (type==1){
  x[,1]<-alphaMAX      # make said column Amax
  }else if (type==2){
  diag(x)<-alphaMAX }
  return(x)})
  STORE_Alphas<-AllAlphas

  if (SuppressAll==FALSE) pb <- txtProgressBar(min = 0, max = M, style = 3)

  ptmZZZ <- proc.time()# REMOVE ME

  for (m in 1:M){
  if (type==3){
  AllAlphas<- lapply(c(1:J), function(x) matrix(Alpha_lows[x],ncol=K, nrow=K)  )
  if(sample( c(1, 0), size=1, prob=c(0.5, 0.5))==1){   # Put on diagonal
  AllAlphas<- lapply(AllAlphas,  function(x) {diag(x)<-alphaMAX ; return(x)})
  } else {
  AllAlphas<-lapply( AllAlphas, function(x) {x[,1]<-alphaMAX;return(x)})
  }}


  if(m %% 10==0){  if (SuppressAll=="FALSE"){
  Sys.sleep(0.1)
  setTxtProgressBar(pb, m)}}

  # OVER EACH CHAINS
  # Transition
  if (m==1) {nt<-lapply(c(1:J),function(x) CountTrans(states0[[x]], K))
  } else { nt<-  lapply(c(1:J), function(j) CountTrans(Z[j,],K))
  Qold<-Q
  Q0old<-q0}
  qnew<-lapply(c(1:J), function(j) sapply(c(1:K), function(x) rdirichlet(par=nt[[j]][x,]+AllAlphas[[j]][x,]) ) )
  qnew<-lapply(qnew, t)
  q0new <-  lapply(qnew, getq0_trycatch)
  if(m==1) {Q<-  t(sapply(qnew, function(x) as.vector(t(x))  ) )# Save values at this iteration
  q0<- t(sapply(qnew, getq0_trycatch) )
  Qold<-Q
  Q0old<-q0}
  if (m>1){
  Qtest<- lapply(c(1:J),function(j) Funkme_MetHast_Q(.q0new=q0new[[j]],.q0old=Q0old[j,],.z1=Z[j,1]) ) # NOT DONE! # NEED TO APPLY OVER right values, particularly Z's (1'st)
  Q<- lapply(c(1:J), function(j) { if (Qtest[[j]]=="TRUE") {return( qnew[[j]])
  }else { return( matrix(Qold[j,], nrow=K, byrow="TRUE") )} })
  q0<-t(sapply(c(1:J), function(j) { if (Qtest[[j]]=="TRUE") {return( q0new[[j]])
  }else { return(Q0old[j,]  )} }))

  Q<-t(sapply(Q, function(x) as.vector(t(x)) ))
  }

  # Means
  if (m==1) { sumNcount<-lapply(c(1:J), function(j) formu2(states0[[j]][-(n+1)],Y,K))
  } else {    sumNcount<- lapply(c(1:J), function(j) formu2(Z[j,-(n+1)], Y,K))}   # MAKE LIST over Chains (RIGHT Z's)
  MU<-t(sapply(c(1:J), function(j) apply(sumNcount[[j]], 1,  function (x)     rnorm(1, mean= ((mu0/var0)+(x[1]/varknown)) / ((1/var0)+(x[2]/varknown)), sd= sqrt( 1/( (1/var0) + (x[2]/varknown)))  ))))

  # States
  newZ<- lapply(c(1:J), function(j) UpdateStates( Y, Q[j,], MU[j,], initS= q0[j,], m))
  Z<- t(sapply(c(1:J) , function(x) newZ[[x]]$Z))
  MAP[m]<-newZ[[J]]$MAP

  ## PT move
  if(m>20 & runif(1)>.1){
  if(runif(1)>.5){chainset<- c(1:(J-1))[c(1:(J-1))%%2==0]   #evens
  } else {chainset<- c(1:(J-1))[c(1:(J-1))%%2!=0] }   #odds

  chainset<-cbind(chainset, chainset+1)

  TrackParallelTemp[TrackParallelTemp$Chain%in%as.numeric(chainset),2]<-TrackParallelTemp[TrackParallelTemp$Chain%in%as.numeric(chainset),2]+1

  if(class(chainset)=='numeric')chainset<-t(data.frame(chainset))
  SwitchMe<-mapply( function(x,y) parallelAcceptHMM_QV2(Q[x,], Q[y,], AllAlphas[[x]] ,AllAlphas[[y]] , Z[x,1], Z[y,1]), chainset[,1], chainset[,2])
  chainset<-chainset[SwitchMe==1,]
  if(class(chainset)=='numeric')chainset<-t(data.frame(chainset))

  TrackParallelTemp[TrackParallelTemp$Chain%in%as.numeric(chainset),3]<-TrackParallelTemp[TrackParallelTemp$Chain%in%as.numeric(chainset),3]+1

  #q0
  .q01<-q0[chainset[,1],]
  .q02<-q0[chainset[,2],]
  q0[chainset[,1],]<-.q02
  q0[chainset[,2],]<-.q01
  #Q
  .Q1<-Q[chainset[,1],]
  .Q2<-Q[chainset[,2],]
  Q[chainset[,1],]<-.Q2
  Q[chainset[,2],]<-.Q1
  #Mu
  .MU1<-MU[chainset[,1],]
  .MU2<-MU[chainset[,2],]
  MU[chainset[,1],]<-.MU2
  MU[chainset[,2],]<-.MU1
  #Z
  .Z1<-Z[chainset[,1],]
  .Z2<-Z[chainset[,2],]
  Z[chainset[,1],]<-.Z2
  Z[chainset[,2],]<-.Z1
  #Qold
  .Qold1<-Qold[chainset[,1],]
  .Qold2<-Qold[chainset[,2],]
  Qold[chainset[,1],]<-.Qold2
  Qold[chainset[,2],]<-.Qold1
  }

  # NOW FINAL BITS
  ## Compute density f2 at this iteration (L1 norm)
if(is.null(TV)==FALSE){
fDist[m]<-L1Norm(Y, TV ,NOW=list("z"=Z[J,],   "q0"=q0[J,], "Q"=Q[J,],  "mu"= MU[J,]) , gridN)
} else {fDist<-"NA"}



  Finalq0[m,]<-q0[J,]
  FinalQ[m,]<-Q[J,]
  FinalMu[m,]<-MU[J,]
  FinalStates[m,]<-Z[J,]
  SteadyScore$K0[m]<-sum(table(Z[J,])>0)
  } # end of iteration loop

  if (SuppressAll=="FALSE") close(pb)

  AnyMix<-0
  if( sum(TrackParallelTemp$NumSuccess > 5) ==J) AnyMix<-1  # if all chains mix a bit

  #find worst mixed chain :
  ChainProportionSuccess<- TrackParallelTemp$NumSuccess/TrackParallelTemp$NumTries
  WorstMixProp<-min(ChainProportionSuccess)

  #print(proc.time() - ptmZZZ)# REMOVE ME


TrackResults<-data.frame("K0"=SteadyScore$K0,"f2Dist"=fDist, "WorstMixProp"=WorstMixProp )
allResults<-list("Means"=FinalMu, "Trans"=FinalQ, "States"=FinalStates, "q0"=Finalq0, "YZ"=YZ, "MAP"=MAP, "K0"=SteadyScore$K0, "f2dens"=f2now, "f2Dist"=fDist,"TrackParallelTemp" =TrackParallelTemp, "Track"=TrackResults)
  return(allResults)
  }

