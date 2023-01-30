cdfPen.fit <- function(b, b.tld, g, b.rho, H.rho, lmb.rho, nu, algorithm, nstep = 1E+5, eps = 1E-5, eps.lla = 1E-6, nstep.lla = 1E+5) {
   # Description:
   # b, b.tld, g = Starting values
   # b.rho = soluzione ridge
   # H.rho = rho * {X^tX / n + diag(rho)}^(-1)
   # lmb.rho = lambda / rho
   # nu = weight
   # nstep = number of iterations
   # eps = threshold
   
   p <- length(b)
   conv <- 1
   for (i in seq_len(nstep)) {
      # step 1: update b
      b <- drop(b.rho + H.rho %*% (b.tld - g))
      # step 2: update b.tld
      b.tld[1L] <- b[1L] + g[1L] # calcolo semplificato dato che l'intercetta non è penalizzata
      if (algorithm == "opt") {
         for (m in 2:p) {
            bm_gm <- b[m] + g[m]
            if (abs(bm_gm) <= lmb.rho) {
               b.tld[m] <- 0
            } else {
               out.opt <- optimize(f = function(bm.tld) 0.5 * (bm_gm - bm.tld)^2 + lmb.rho * sqrt(2 * pi) * nu * pnorm(abs(bm.tld) / nu),
                                 interval = if(bm_gm > 0) c(0, 1E+6) else c(-1E+6, 0))
               b.tld[m] <- out.opt$minimum
            }
         }
      }
      else if (algorithm == "lla") {
         bm_gm <- b[2:p] + g[2:p]
         out <- lla(b.o = b.tld[2:p], lmb.rho = lmb.rho, bm_gm = bm_gm, nu = nu, nstep.lla = nstep.lla, eps.lla = eps.lla)
         if (out$conv == 1) stop("lla does not converge")
         b.tld[2:p] <- out$b
      }
      # step 3: update g
      g <- g + (b - b.tld)
      # check convergence
      sqrt.ave.r2 <- sqrt(mean((b - b.tld)^2))
      sqrt.ave.r2
      if (sqrt.ave.r2 <= eps) {
         conv <- 0
         break
      }
   }
   
   out <- list(b = b, b.tld = b.tld, g = g, i = i, conv = conv)
   out
}
