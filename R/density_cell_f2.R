 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



density_cell_f2<-function( y1, y2, q01, Q12, mu1, mu2  ){
return( q01*Q12*dnorm(y1, mean=mu1 )*dnorm(y2, mean=mu2))
}