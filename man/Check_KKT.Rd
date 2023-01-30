\name{check_KKT}
\alias{check_KKT}
\title{Check on the condition of Karush-Kuhn-Tucker
}
\description{
Control over Karush-Kuhn-Tucker (Karush, 1939) conditions for the estimates obtained. }

\usage{
check_KKT(obj,
          intercept = TRUE)
}
\arguments{
  \item{obj}{Object to be checked.}
  \item{intercept}{Is the intercept used in the model?}

}


\value{
  \item{grd}{The value of gradient.}
  \item{hx}{The value of equality constraint.}
  \item{glob}{The global value of derivative (grd + hx).}
  \item{test}{Is the condition verified?}
  \item{lmb}{The values of lambda used in the model}


}
\references{
Karush, W. (1939). Minima of functions of several variables with inequalities as side constraints. M. Sc. Dissertation. Dept. of Mathematics, Univ. of Chicago.
}
\author{
Daniele Cuntrera, Luigi Augugliaro, Vito Muggeo
}
\examples{

p <- 10
n <- 100
X <- cbind(1, matrix(rnorm(n * p), n , p))
b.s <- c(1, rep(0, p))
b.s[sample(2:p, 3)] <- 1
y <- drop(crossprod(t(X), b.s))
out <- cdfPen(X = X, y = y)

KKT <- check_KKT(out)
plot(KKT$test)


}

