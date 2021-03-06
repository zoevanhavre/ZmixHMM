% Generated by roxygen2 (4.0.2): do not edit by hand
\name{gibbsHMM_Main}
\alias{gibbsHMM_Main}
\title{gibbsHMM_PT}
\usage{
gibbsHMM_Main(YZ, M = 2000, K = 10, alphaMin = 1e-05, alphaMAX = 1,
  type = 1, J = 1, SuppressAll = FALSE, TV = NULL, gridN = 100)
}
\arguments{
\item{x,}{alpha, log=False}
}
\description{
parallel tempering with a column prior - option to mix over column or stick to j=1
}
\examples{
dDirichlet(c(.1, .9), c(0.1,0.1))
}
\keyword{dirichlet}

