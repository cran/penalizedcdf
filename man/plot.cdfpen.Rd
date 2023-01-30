\name{plot_cdfpen}
\alias{plot_cdfpen}

\title{Plot coefficients or BIC from a "cdfpen" object

}

\description{
Plot coefficient profile plot or BIC trend
}

\usage{
plot_cdfpen(object,
            ...)
}

\arguments{
  \item{object}{Object to be plotted.}
  \item{...}{Other graphical parameters to plot.}
}

\details{
A graph showing the BIC trend or profile of coefficients is displayed.}



\author{
Daniele Cuntrera, Luigi Augugliaro, Vito Muggeo}

\value{No return value}

\examples{

p <- 10
n <- 100
X <- cbind(1, matrix(rnorm(n * p), n , p))
b.s <- c(1, rep(0, p))
b.s[sample(2:p, 3)] <- 1
y <- drop(crossprod(t(X), b.s))
out <- cdfPen(X = X, y = y)

plot_cdfpen(out)         #Coefficients' path ~ lambda
plot_cdfpen(out, "l1")   #Coefficients' path ~ L1 norm
plot_cdfpen(out, "BIC")  #BIC ~ lambda }

