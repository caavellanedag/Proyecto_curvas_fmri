scores=function(X){

     if (inherits(X,"SpatFD") ) {

          puntaje=list()
          puntajes=list()
          for(j in 1:length(X)){


               puntaje[[j]]=X[[j]]$fpca$scores
               rownames(puntaje[[j]])=X[[j]]$coordsnames
               puntajes[[j]]=as.data.frame(puntaje[[j]])
          }
     } else if(inherits(X,"scores_pred")){
          puntajes = X$scores_pred
     } else {
          stop("Wrong class of X object. It must be of class SpatFD or scores_pred")
     }
     return(puntajes)
}
