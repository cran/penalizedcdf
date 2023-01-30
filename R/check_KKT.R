check_KKT <- function(obj, intercept  = TRUE) {

  X <- obj$X
  y <- obj$y
  b <- obj$coefficients
  lmb <- obj$lmb
  nu <- obj$nu ; if(length(nu) == 1) nu <- rep(nu, dim(b)[2L])
  n <- dim(X)[1L]

  grd  <- matrix(ncol = dim(b)[2L], nrow = dim(b)[1L])
  hx   <- matrix(ncol = dim(b)[2L], nrow = dim(b)[1L])
  glob <- matrix(ncol = dim(b)[2L], nrow = dim(b)[1L])
  test <- matrix(ncol = dim(b)[2L], nrow = dim(b)[1L])

  colnames(grd) <- colnames(hx) <- colnames(glob) <- colnames(test) <- round(lmb,3)


  for(i in seq_len(dim(b)[2L])){
    r <- y - X %*% b[,i]
    grd[,i] <- -drop(crossprod(X, r)) / n
    hx[,i] <- lmb[i]*exp(-(b[,i]^2)/(2*nu[i]^2))*sign(b[,i])
    glob[,i] <- grd[,i] + hx[,i]
    test[,i] <- ifelse(b[,i] == 0 , abs(grd[,i])<=lmb[i], abs(glob[,i])<=1E-3)
  }

  if(intercept) test[1,] <- abs(grd[1,]) <= 1E-3

  out <- list(grd = grd, hx = hx, glob = glob, test = test, lmb = lmb)
  out

}
