KS_scores=function(SFD, newcoords,model,name=NULL,fill.all=NULL){

     # Validation --------------------------------------------------------------

     #all
     if(missing(SFD)){
          stop("Missing SFD")
     }
     if (missing(newcoords)){
          stop("Missing new coords")
     }
     if(missing(model)){
          stop("Missing model")
     }
     #SFD
     if(!inherits(SFD,"SpatFD")){
          stop("SFD must be an object SpatFD")
     }

     #newcoords
     if(!(is.matrix(newcoords) || is.data.frame(newcoords))){
          stop("Wrong class of newcoords object")
     }else if(!all(apply(newcoords, c(1,2), is.numeric))){
          stop("Newcoords must be numeric data")
     }else if(any(is.na(newcoords))){
          stop("There is some NA value in newcoords")
     }

     # messages default values
     if(missing(name)){
          message("Using first variable by default")
     }
     if(missing(fill.all)){
          message("Using fill.all = TRUE by default")
     }

     #vari
     if(is.null(name)){
          name=1
     } else if ((is.character(name)&& length(name)==1)){
          if (length(which(names(SFD)==name))==1){
               name=which(names(SFD)==name)
          }else if (length(which(names(SFD)==name))==0){
               stop(paste(name,"doesn't not exists. Change name for an existing nameable name."))
          }else if (length(which(names(SFD)==name))==0){
               stop("There are more than one nameable with the same name")
          }
     }
     if ((is.null(name)  || !(is.numeric(name)&& length(name)==1))){
          stop("Wrong class of name object")
     }

     #fill.all
     if(is.null(fill.all)){
          fill.all=TRUE
     }else if ( !( ( isTRUE(fill.all) || isFALSE(fill.all) ) && length(fill.all)==1 ) ){
          stop("Wrong class of fill.all object")
     }


     #model
     if(!(inherits(model,"variogramModel") || inherits(model,"list"))){
          stop("Wrong class of model, model should be of class variogramModel or a list of them (use vgm of gstat package) ")
     }else if(inherits(model,"list") && !all(lapply(model,inherits,"variogramModel"))){
          stop("Wrong class of model, each element of list should be of class variogramModel (use vgm of gstat package)")
     }else if(inherits(model,"list") && (length(model)!=ncol(as.data.frame(SFD[[name]]$fpca$scores)))){
          stop("length of list of models must be equal to number of harmonics of the choosen variable ")
     }else if(inherits(model,"variogramModel") && !(fill.all || (ncol(as.data.frame(SFD[[name]]$fpca$scores))==1))){
          stop("If model is not a list and there are more than one nharm of that variable, then fill.all must be TRUE or you can create a list of models with the same number of harmonics")
     }


     # Kriging -----------------------------------------------------------------

     #scores
     puntaje=SFD[[name]]$fpca$scores
     rownames(puntaje)=SFD[[name]]$coordsnames
     puntajes=as.data.frame(puntaje)
     coordinates(puntajes)=SFD[[name]]$coords

     #fitting variogram
     if(fill.all==T){
          v=list()
          fv=list()
          for (i in 1:ncol(puntajes)){

               v[[i]]=variogram(puntajes[[i]]~1,puntajes)
               fv[[i]]=fit.variogram(v[[i]],model)
          }
     }else{
          v=list()
          fv=list()
          for (i in 1:ncol(puntajes)){

               v[[i]]=variogram(puntajes[[i]]~1,puntajes)
               fv[[i]]=fit.variogram(v[[i]],model[[i]])
          }
     }

     #newcoords
     colnames(newcoords)=c('x','y')
     coordinates(newcoords)=~x+y

     #kriging
     K=list()
     for (i in 1:ncol(puntajes)){
          K[[i]] <- krige(puntajes[[i]]~1,puntajes,newcoords, model = fv[[i]])
     }

     #prediction
     pred=K[[1]]$var1.pred
     if(ncol(puntajes)>1){
          for (i in 2:ncol(puntajes)){
               pred=cbind(pred,K[[i]]$var1.pred)
          }
     }
     pred=as.data.frame(pred)
     colnames(pred)[1]="V1"
     rownames(pred)=rownames(newcoords)

     #Output
     out=list(SFD=SFD,scores_pred=pred,name=name)
     class(out)="scores_pred"
     return(out)
}
