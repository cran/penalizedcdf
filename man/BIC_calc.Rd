\name{BIC_calc}
\alias{BIC_calc}

\title{
BIC calculator function}
\description{
Function that takes the resulting values of the estimated model as input, to compute BIC
}
\usage{
BIC_calc(X,
         b.tld,
         y,
         n)
}

\arguments{
  \item{X}{The covariates' matrix}
  \item{b.tld}{The estimated sparse-beta}
  \item{y}{The response variable}
  \item{n}{The number of observation}
}


\value{Returns the BIC value calculated for a single value of the tuning parameter.}

\examples{

p <- 10
n <- 100
X <- cbind(1, matrix(rnorm(n * p), n , p))
b.s <- c(1, rep(0, p))
b.s[sample(2:p, 3)] <- 1
y <- drop(crossprod(t(X), b.s))
out <- cdfPen(X = X, y = y)


(bic <- BIC_cdfpen(out))
plot(out$lmb, bic, "s")}
