recons_fd = function(X){
  
  if (inherits(X,"scores_pred")) {
    vari=X$vari
    a=scores(X)
    mean_coef = X$SFD[[vari]]$fpca$meanfd$coefs
    nr=nrow(a)
    mean_coef = matrix(rep(mean_coef,nr),ncol = nr)
    coef = (X$SFD[[vari]]$fpca$harmonics$coef %*% t(a)) + mean_coef
    result = fd(coef, X$SFD[[vari]]$fpca$harmonics$basis)
  }else{
    stop("Wrong class scores_pred")
  }
  return(result)
}

