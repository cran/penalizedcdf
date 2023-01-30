plot_path <- function(obj, lmb, coeff, type = c("path", "l1", "BIC"), ...){

  type <- match.arg(type)

  if(type == "path") matplot(lmb, t(coeff[-1,]), type = "l", ylab = expression(beta), xlab = expression(lambda), ...)
    else if(type == "l1") matplot(colSums(abs(coeff[-1,])), t(coeff[-1,]), type = "l", ylab = expression(beta), xlab = "L1 norm", ...)
    else if(type == "BIC") plot(lmb, BIC_cdfpen(obj), type="l", ylab = "BIC", xlab = expression(lambda))
}

plot_cdfpen<-function(object, ...){

  plot_path(object, object$lmb, object$coefficients, ...)

}
