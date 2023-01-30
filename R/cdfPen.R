cdfPen <- function(X, y, nu, lmb, nlmb = 100L, e = 1E-3, rho = 2, algorithm = c("lla", "opt"), nstep = 1E+5, eps = 1E-6,
                   eps.lla = 1E-6, nstep.lla = 1E+5) {

   algorithm <- match.arg(algorithm)
   n <- dim(X)[1L]
   p <- dim(X)[2L]

   if (missing(lmb)) {
      lmb.max <- max(abs(drop(crossprod(X, y - mean(y)))) / n)
      lmb <- seq(from = lmb.max * (1 + 1E-6), to = e * lmb.max, length = nlmb)
   } else {
      lmb <- sort(lmb, decreasing = TRUE)
      nlmb <- length(lmb)
   }

   if (missing(nu)) nu <- lmb * exp(-0.5) / rho
      else if(length(nu) == 1) nu <- rep(nu, length(lmb))

   B <- matrix(0, nrow = p, ncol = nlmb)

   # computing starting values
   H.rho <- crossprod(X) / n
   diag(H.rho) <- diag(H.rho) + rho
   H.rho <- solve(H.rho)
   b.rho <- drop(H.rho %*% crossprod(X, y) / n)

   b.tld <- c(mean(y), rep(0, p - 1))
   b <- b.tld
   g <- b.tld - rho^(-1) * drop(solve(H.rho) %*% (b.tld - b.rho))

   H.rho <- rho * H.rho
   for (h in seq_along(lmb)) {
      out <- cdfPen.fit(b = b, b.tld = b.tld, g = g, b.rho = b.rho,
                        H.rho = H.rho, lmb.rho = lmb[h] / rho,
                        nu = nu[h], algorithm = algorithm, nstep = nstep,
                        eps = eps, eps.lla = eps.lla, nstep.lla = nstep.lla)
      if(out$conv == 1) warning("cdfPen.fit does not converge")
      B[, h] <- out$b.tld
      b <- out$b
      b.tld <- out$b.tld
      g <- out$g
   }

   out <- list(coefficients = B, lmb = lmb, rho = rho, nu = nu, X = X, y = y, algorithm = algorithm)

   class(out)<-"cdfpen"
   out
}
