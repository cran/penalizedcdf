lla <- function(b.o, lmb.rho, bm_gm, nu, nstep.lla = 100L, eps.lla = 1E-6) {
  conv <- 1
  for (i in seq_len(nstep.lla)) {
    w <- exp( - b.o^2 / (2 * nu^2))
    b.n <- S(bm_gm, lmb.rho, w)
    check <- sqrt(sum((b.n - b.o)^2))
    if (check <= eps.lla) {
      conv <- 0
      break
    }
    b.o <- b.n
  }
  out <- list(b = b.n, conv = conv, nstep.lla = i)
  out
}

S <- function(bm_gm, db, w) {
   ifelse (abs(bm_gm) <= db * w, 0, sign(bm_gm) * (abs(bm_gm) - db * w))
}