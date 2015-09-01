#' gibbsHMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

makeStartSimpler<-function(Y,k=K){
          n<-length(Y)
          states0<-rep(0,n)

            # split data by quantiles &  compute allocations
            states0<-as.numeric(cut2(Y, m=5, g=k))
           
        # randomize label names
        states0<-as.factor(states0)
        levels(states0)<-sample(levels(states0), k)
        states0<-as.numeric(     as.character(states0))
            # initstates0<-sample(c(1:k), size=1, prob=rep(1/k,k))
            # states0<-c(initstates0, states0)
           
            return(states0)
          }