#################################################################
### Cargue de librerías
#################################################################

library(pacman)
p_load(tidyverse, openxlsx, here, haven, data.table, fda, fda.usc)
p_load(maptools, sgeostat, scatterplot3d, car, fields, gstat, geoR)       

#################################################################
### Funciones construidas para el ejercicio
#################################################################

lectura_datos <- function(.x){
  datos_sujeto <- read.xlsx(paste0("input/",.x))
  label_subject <- str_extract(.x,"\\d+")
  
  datos_sujeto <- datos_sujeto %>% as.data.table()
  
  datos_sujeto[,c("MEDICION","REPETICION"):=list(rep(1:228,(nrow(datos_sujeto)/228)),
                                                 rep(1:(nrow(datos_sujeto)/228),each=228))]
  
  datos_sujeto_1 <- datos_sujeto %>%
    melt(id.vars = c("FRECUENCIA","CLASE","MEDICION","REPETICION"),variable.name="ELECTRODO",value.name = "VALOR")
  datos_sujeto_1[,c("SUJETO","KEY"):=list(paste0("S",label_subject),
                                          paste0("S",label_subject,"_",REPETICION,"_",ELECTRODO))]
  #class<-datos_sujeto_1[,c("REPETICION","CLASE")] %>% unique() 
  #datos_sujeto_1 %>% ggplot()+geom_line(aes(MEDICION,VALOR,group=REPETICION,color=ELECTRODO))+facet_wrap(~CLASE,scales="free")
  return(datos_sujeto_1)
}

SelectionNumBasis <- function(x, Y_mat, Num_basis){
  spline_basis <- create.bspline.basis(rangeval = c(min(x), max(x)), nbasis = Num_basis)
  object_fd <- smooth.basis(y = Y_mat, fdParobj = spline_basis)
  predict_fd <- eval.fd(x, fdobj = object_fd$fd)
  SSE <- sum((predict_fd-Y_mat)^2)
  CME <- SSE/(length(x)-Num_basis)
  X <- eval.basis(x, spline_basis)
  Var_estimator <- sum(diag(X%*% solve(t(X) %*% X) %*% t(X)))*CME
  return(data.frame(Num_basis = Num_basis,
                    Bias_squared = SSE,
                    Var_estimator = Var_estimator,
                    MSE = Var_estimator + SSE))
}


Select_lambda <- function(x, Y_mat, lambda){
  lambda_10 <- 10^lambda
  spline_basis <- create.bspline.basis(rangeval=c(min(x), max(x)), nbasis = 50)
  fdParobj <- fdPar(fdobj = spline_basis, Lfdobj = 2, lambda = lambda_10)
  object_fd <- smooth.basis(y = Y_mat, fdParobj = fdParobj)
  #mean((eval.fd(1:100, fdobj=absorb_fd$fd) - t(absorp)[,.y])^2)
  data.frame(lambda = lambda, gcv = sum(object_fd$gcv))
}


smooth_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  predict_fd <- eval.fd(x, fdobj=object_fd$fd)  
  return(predict_fd)
}

fitted_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=50)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  return(object_fd)
}



# df_num_basis <- map_dfr(4:50,~SelectionNumBasis(x = 1:100, Y_mat = t(absorp), Num_basis = .x))
# source("src/Funciones.R")

#################################################################
### Lectura de las coordenadas
#################################################################

coordenadas <- fread("input/coor1.txt")
coordenadas <- coordenadas %>%
  cbind(data.frame(ELECTRODO = paste0("E",1:21))) %>% 
  mutate(z=as.numeric(as.character(str_replace_all(z,",","\\."))))

#################################################################
### Lectura de la información por cada uno de los sujetos,
## mediante la función "lectura_datos"
#################################################################
files_sujetos <- list.files(path="input",pattern=".xlsx")
list <- map(files_sujetos, lectura_datos)

#################################################################################
##### Gráfico de resultados para el k-ésimo individuo
#################################################################################

k <- 1
list[[k]] %>%
  filter(REPETICION == 1) %>% 
  ggplot() +
  geom_line(aes(MEDICION, VALOR, group = KEY)) +
  facet_grid(ELECTRODO ~ CLASE) 

#################################################################################
##### Gráfico de resultados para el k-ésimo individuo
#################################################################################

lista_wider <- list[[k]] %>% 
  filter(REPETICION == 1) %>% 
  dplyr::select(-FRECUENCIA) %>% 
  pivot_wider(id_cols = c("CLASE", "REPETICION", "ELECTRODO", "SUJETO"),
              names_from = MEDICION, values_from = VALOR, 
              names_prefix = "MEDICION_")
lista_wider_2 <- lista_wider %>% 
  dplyr::select(starts_with("MEDICION")) %>% 
  as.matrix()

#################################################################################
##### Gráfico del número de bases
#################################################################################

x <- 1:228
df_num_basis <- map_dfr(5:20, ~SelectionNumBasis(x = 1:228, Y_mat = t(lista_wider_2), .x))
df_num_basis %>% 
  melt(id.vars = "Num_basis", variable.name = "Variable", value.name = "Value") %>% 
  ggplot() + 
  geom_line(aes(x = Num_basis, Value, linetype = Variable, color = Variable), size = 1.1)+
  scale_color_brewer(palette="Set1")+theme_bw()+ylab("Cuadrado medio del error total")+
  xlab("Número de funciones base")+ylim(0,5)

#################################################################################
##### Se encuentra la estimación del parámetro lambda mediante gcv
#################################################################################

df_gcv <- map_dfr(seq(-5,5,by=0.01),~Select_lambda(x= x , Y_mat = t(lista_wider_2), lambda=.x)) 
lambda_optimo <- 10^df_gcv[which.min(df_gcv$gcv),"lambda"]

#################################################################################
##### Generación de curvas estimadas
#################################################################################

#smooth_fda(x, Y_mat = t(lista_wider_2), lambda = lambda_optimo) %>% dim
fitted_fda <- fitted_fda(x, Y_mat = t(lista_wider_2), lambda = lambda_optimo)


total_var_prop <- 0
k <- 1
while(total_var_prop < 0.85){
  pca_overall <- pca.fd(fitted_fda$fd, nharm = k)
  total_var_prop <- sum(pca_overall$varprop)
  k <- k + 1 
}

num_components <- dim(pca_overall$scores)[2]
total_var_prop

pca_overall %>% names
pca_overall$harmonics


t(lista_wider_2) %>% dim
nrow(coordenadas)
SpatFD(data = fitted_fda,
       coords = coordenadas[, 1:3],
       basis="Bsplines",
       nbasis = 4,
       lambda=lambda_optimo,
       nharm = 5,
       vp=0.85)
