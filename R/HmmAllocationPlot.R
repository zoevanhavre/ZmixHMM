#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope

HmmAllocationPlot<-function( outZ, myY){
			grr<-outZ[,order(myY)]
			grrTable<-data.frame("myY"=NULL, "k"=NULL, "Prob"=NULL)
			#maxK<- length(levels(grr))
			maxK<-max(grr)
			for (i in 1:length(myY)){
			#rr<-factor(grr[i], levels=1:maxK)
			rr<-factor(grr[,i], levels=1:maxK)
			grrTable<-rbind(grrTable,cbind(i,c(1:maxK), matrix(table(rr)/ length(rr) )))
			    }
			names(grrTable)<-c("myY", "k", "Prob")
				grrTable$k<-as.factor(grrTable$k)

			gp<-ggplot(grrTable, aes(x=myY, y=k, fill=Prob)) + 
			geom_tile()+
			theme(legend.position='right', plot.title=element_text(size=1,vjust=2))+
			ggtitle("Allocation Probabilities")+ 
			xlab("index of ordered y")+
			scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8","#253494" ))+
			theme_bw()
			#ggsave( plot=gp,  filename=paste( "Allocations_", plotfilename ,"K_",maxK, ".pdf",sep="") )
			gp
			}

	