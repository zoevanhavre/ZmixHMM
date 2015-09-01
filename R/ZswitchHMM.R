#' PostProcessing function univ
#'
#' This function draws samples from a Wishart distu
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope


Zswitch_hmm<-function(GrunK0, PropMin=0.01 ){
			K<-dim(GrunK0$q0)[2]
			
			 # Pick Reference 
			wml<-which.max(GrunK0$MAP)
			Zref<-factor(GrunK0$Z[wml,], levels=1:K,ordered=FALSE) 
			.zr<-table(Zref)
			.zr<-as.numeric(names(.zr[as.numeric(.zr)!=0]))
			FinalOrderChoice<-.zr[order(GrunK0$Mu[wml,.zr])]		
			refComp<-c(GrunK0$q0[wml,FinalOrderChoice], GrunK0$Mu[wml,FinalOrderChoice])  
			Zref<- factor(Zref)

			# RENAME TO neat numbers
			levels(Zref)<- c(1:K)[order(FinalOrderChoice[1:length(.zr)])]
			Zref<-factor(as.numeric(as.character(Zref)))
			# storage dataframes:
			numK0now<-sum(table(Zref )>0) 
			numPar<-length(unlist(lapply(1:dim(GrunK0$Z)[1], rep, numK0now)))
			AllPars<-data.frame(matrix(0, ncol=4+numK0now, nrow=numPar) )
			AllPars[,1]<-unlist(lapply(1:dim(GrunK0$Z)[1], rep, numK0now))
			Zfixed<-GrunK0$Z
			#for each iteration
			for(.iter in 1:dim(GrunK0$q0)[1]){
				#Store current states
				Znow<-factor(GrunK0$Z[.iter,])    
				
				#identify potential candidate switches:
				CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
				getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
				ListCandi<- apply(CandiCells, 1, getCandi)
				
				# R stuff to make sure it deals with inputs correctly
				if(class(ListCandi)=='matrix'){
					ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
					Candies<-expand.grid(ListCandi)  # each row is a labelling
					names(Candies)<-row.names(CandiCells)   # RAAAAH
					} else if (class(ListCandi)=='numeric'){
					Candies<-ListCandi
					} else {
					Candies<-expand.grid(ListCandi)  # each row is a labelling
					names(Candies)<-row.names(CandiCells)   # RAAAAH
					}

				namesCandies<-names(Candies)
		# Catch if no appropriate seperation available 
				done<-0 
				if(class(Candies)=='data.frame'){
				if  ( max(sapply(apply(Candies, 1, unique), length))<length(row.names(CandiCells))){
				Candies<- matrix( as.numeric( row.names(CandiCells))[permutations(numK0now )], ncol=numK0now)	
				colnames(Candies)<-as.numeric(names(as.data.frame(CandiCells)))
					
				MinusRefPars_catch<-function(x) 	{ flp<- Candies[x,]   
					if(length(unique(flp))<length(flp)) {	 Inf 
					} else {sum(abs( (refComp	-  c(GrunK0$q0[.iter,flp], GrunK0$Mu[.iter,flp]))/refComp))	}}
				BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars_catch))  # find the best perm out of options
				BestOne<-Candies[BestOne,]
				flp<-as.numeric(BestOne)[as.numeric(names(BestOne))]
					# Allocations
				Znew<-Znow
				levels(Znew)<-as.numeric(names(BestOne))
				Zfixed[.iter,]<-as.numeric(as.character(Znew))
					# Parameters
				swQ<-matrix(GrunK0$Q[.iter,], nrow=K, byrow=T)[flp, flp]					
				combinePars<- cbind(.iter, 1:numK0now,  GrunK0$q0[.iter,flp],GrunK0$Mu[.iter,flp],swQ)
					AllPars[AllPars[,1]==.iter,]<- combinePars ;	done<-1
				}}

		# Normally:					
			if (done==0){
					MinusRefPars<-function(x) 	{	
						#problem here (below)
						flp<- na.omit(as.numeric(  row.names(CandiCells)[unlist(Candies[x,])]))
						if(length(unique(flp))<length(flp)) { Inf
						} else {
						flp1<-as.numeric(names(Candies))	
						flp2<-flp1[order(unlist(Candies[x,]))]				
						sum(abs( (refComp	-  c(GrunK0$q0[.iter,flp2], GrunK0$Mu[.iter,flp2]))/refComp))	
						}}
											
				if( sum( apply(CandiCells, 1, sum)) >  dim(CandiCells)[1] ){
					BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars))  # find the best perm out of options
					BestOne<-Candies[BestOne,]
					} else {BestOne<- Candies }   # chose this one if no comparing needed
				if(is.null(names(BestOne))) {names(BestOne)<-namesCandies}
	
				# Allocations
				Znew<-Znow; levels(Znew)<-as.numeric(BestOne)
				Zfixed[.iter,]<-as.numeric(as.character(Znew))
				
			newOrder<-order(as.numeric(BestOne), decreasing=FALSE)
			combinePars<-cbind(.iter,as.numeric(BestOne),  GrunK0$q0[.iter,as.numeric(names(BestOne))],GrunK0$Mu[.iter,as.numeric(names(BestOne))] )[newOrder,]
			swQ<-matrix(GrunK0$Q[.iter,], nrow=K, byrow=T)[as.numeric(names(BestOne)),as.numeric(names(BestOne))]
			if(class(swQ)!="numeric") swQ<-swQ[newOrder,newOrder]			
			combinePars<-cbind(combinePars,swQ)
			AllPars[AllPars[1]==.iter,]<- combinePars
			}		
		}
		
		Zhat<- factor( apply(Zfixed, 2,maxZ))   # sumarise Z
		names(AllPars)[1:4]<-c("Iteration", "k", "q0", "mu")
		return(list('Pars'=AllPars, 'Z'=Zfixed))
		}