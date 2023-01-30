\name{cdfPen}
\alias{cdfPen}
\title{Fit a Linear Model with with CDF regularization
}
\description{
Uses the CDF penalty to estimate a linear model with the maximum penalized likelihood. The path of coefficients is computed for a grid of values for the lambda regularization parameter.
}
\usage{
cdfPen(X,
       y,
       nu,
       lmb,
       nlmb = 100L,
       e = 1E-3,
       rho = 2,
       algorithm = c("lla", "opt"),
       nstep = 1E+5,
       eps = 1E-6,
       eps.lla = 1E-6,
       nstep.lla = 1E+5)
}

\arguments{
  \item{X}{Matrix of covariates, each row is a vector of observations. The matrix must not contain the intercept.}
  \item{y}{Vector of response variable.}
  \item{nu}{Shape parameter of the penalty. It affects the degree of the non-convexity of the penalty. If no value is specified, the smallest value               that ensures a single solution will be used.}
  \item{lmb}{A user-supplied tuning parameter sequence.}
  \item{nlmb}{number of lambda values; 100 is the default value.}
  \item{e}{The smallest lambda value, expressed as a percentage of maximum lambda. Default value is .001.}
  \item{rho}{Parameter of the optimization algorithm. Default is 2.}
  \item{algorithm}{Approximation to be used to obtain the sparse solution.}
  \item{nstep}{Maximum number of iterations of the global algorithm.}
  \item{eps}{Convergence threshold of the global algorithm.}
  \item{eps.lla}{Convergence threshold of the LLA-algorithm (if used).}
  \item{nstep.lla}{Maximum number of iterations of the LLA-algorithm (if used).}

}

\details{
We consider a local quadratic approximation of the likelihood to treat the problem as a weighted linear model.

The choice of value assigned to \eqn{\nu} is of fundamental importance: it affects both computational and estimation aspects. It affects the ''degree of non-convexity'' of the penalty and determines which of the good and bad properties of convex and non-convex penalties are obtained. Using a high value of \eqn{\nu} ensures the uniqueness of solution, but the estimates will be biased. Conversely, a small value of \eqn{\nu} guarantees negligible bias in the estimates. The parameter \eqn{\nu} has the role of determining the convergence rate of non-null estimates$: the lower the value, the higher the convergence rate. Using lower values of \eqn{\nu}, the objective function will have local minima.

}
\value{

  \item{coefficients}{The coefficients fit matrix. The number of columns is equal to nlmb, and the number of rows is equal to the number of                                 coefficients.}
  \item{lmb}{The vector of lambda used.}
  \item{e}{The smallest lambda value, expressed as a percentage of maximum lambda. Default value is .001.}
  \item{rho}{The parameter of the optimization algorithm used}
  \item{nu}{The shape parameters of the penalty used.}
  \item{X}{The design matrix.}
  \item{y}{The response.}
  \item{algorithm}{Approximation used}

}

\references{
Aggiungere Arxiv
}
\author{
Daniele Cuntrera, Luigi Augugliaro, Vito Muggeo}

\examples{

p <- 10
n <- 100
X <- cbind(1, matrix(rnorm(n * p), n , p))
b.s <- c(1, rep(0, p))
b.s[sample(2:p, 3)] <- 1
y <- drop(crossprod(t(X), b.s))
out <- cdfPen(X = X, y = y)}


