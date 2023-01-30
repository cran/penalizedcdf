\name{cdfPen.fit}
\alias{cdfPen.fit}
\title{Fitter function for CDF penalty
}
\description{
These are the fundamental computing algorithms that cdfPen invokes to estimate penalized linear models by varying lambda.}

\usage{
cdfPen.fit(b,
           b.tld,
           g,
           b.rho,
           H.rho,
           lmb.rho,
           nu,
           algorithm,
           nstep = 1E+5,
           eps = 1E-5,
           eps.lla = 1E-6,
           nstep.lla = 1E+5)
}

\arguments{

  \item{b}{Starting values of beta-vector.}
  \item{b.tld}{Starting values of sparse beta-vector.}
  \item{g}{Starting values of pseudo-variable.}
  \item{b.rho}{Ridge solution.}
  \item{H.rho}{Second part of ridge solution.}
  \item{lmb.rho}{Lambda-rho ratio.}
  \item{nu}{Shape parameter of the penalty. It affects the degree of the non-convexity of the penalty.}
  \item{algorithm}{Approximation to be used to obtain the sparse solution.}
  \item{nstep}{Maximum number of iterations of the global algorithm.}
  \item{eps}{Convergence threshold of the global algorithm.}
  \item{eps.lla}{Convergence threshold of the LLA-algorithm (if used).}
  \item{nstep.lla}{Maximum number of iterations of the LLA-algorithm (if used).}
}

\value{

  \item{b}{Estimated beta-vector.}
  \item{b.tld}{Estimated sparse beta-vector.}
  \item{g}{Final values of pseudo-variable.}
  \item{i}{Number of iterations.}
  \item{conv}{Convergence check status (0 if converged).}



}
\references{
Aggiungere Arxiv}
\author{
Daniele Cuntrera, Luigi Augugliaro, Vito Muggeo}

