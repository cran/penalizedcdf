BIC_calc <- function(X, b.tld, y, n){

  fit <- X %*% b.tld
  edf <- sum(b.tld!=0)
  s2 <- sum((y-fit)^2)/(n-edf)

  ll <- sum(dnorm(y, fit, sqrt(s2*(n-edf)/n), TRUE))
  BIC <- -2*ll + log(n)*edf


}

BIC_cdfpen<-function(object){
  n<- length(object$y)
  nc<-ncol(object$coefficients)
  r<-NULL
  for(i in 1:nc) r[length(r)+1] <-BIC_calc(object$X, drop(object$coefficients[,i]), object$y, n)
  r
}
