\name{lla}
\alias{lla}
\title{LLA approximation for CDF penalty
}
\description{
Linearly approximate a part of the objective function to greatly speed up computations.}
\usage{
lla(b.o,
    lmb.rho,
    bm_gm,
    nu,
    nstep.lla = 100L,
    eps.lla = 1E-6)
}
\arguments{
  \item{b.o}{Vector of sparse-solution.}
  \item{lmb.rho}{Lambda-rho ratio.}
  \item{bm_gm}{Vector of pseudo-solution}
  \item{nu}{Shape parameter of the penalty.}
  \item{nstep.lla}{Maximum number of iterations of the LLA-algorithm (if used).}
  \item{eps.lla}{Convergence threshhold of the LLA-algorithm (if used).}
}

\details{
The LLA approximation allows the computationally intensive part to be treated as a weighted LASSO (Tibshirani, 1996) problem. In this way the computational effort is significantly less while maintaining satisfactory accuracy of the results. See Zou and Li (2008). }
\value{

  \item{b}{Vector of the estimated sparse-solution.}
  \item{Conv}{Convergence check (0 if converged).}
  \item{nstep.lla}{Number of iterations done.}

}
\references{
Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B (Methodological), 58(1):267–288.

Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave penalized likelihood models. Annals of statistics, 36(4):1509}
