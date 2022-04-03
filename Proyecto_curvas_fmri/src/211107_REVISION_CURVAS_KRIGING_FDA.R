#################################################################
### Cargue de librerías
#################################################################
rm(list = ls())
library(pacman)
p_load(tidyverse, openxlsx, here, haven, data.table, fda, fda.usc)
p_load(maptools, sgeostat, scatterplot3d, car, fields, gstat, geoR) 
p_load(automap)
source("src/SpatFD-main/R/KS_scores.R")
source("src/SpatFD-main/R/recons_fd.R")
source("src/SpatFD-main/R/scores.R")
source("src/SpatFD-main/R/SpatFD.R")
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
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=15)
  fdParobj <- fdPar(fdobj=spline_basis, Lfdobj=2, lambda=10^lambda_optimo)
  object_fd <- smooth.basis(y=Y_mat, fdParobj=fdParobj)
  predict_fd <- eval.fd(x, fdobj=object_fd$fd)  
  return(predict_fd)
}

fitted_fda <- function(x, Y_mat, lambda){
  spline_basis <- create.bspline.basis(rangeval=c(min(x),max(x)),nbasis=15)
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

mins_coords <- apply(coordenadas[, 1:3], 2, min)-1
max_coords <- apply(coordenadas[, 1:3], 2, max)+1

grid <- expand.grid(seq(mins_coords[1], max_coords[1], length.out = 20),
                    seq(mins_coords[2], max_coords[2], length.out = 20)
                    #,seq(mins_coords[3], max_coords[3], length.out = 20)
                    )



#################################################################
### Lectura de la información por cada uno de los sujetos,
## mediante la función "lectura_datos"
#################################################################
files_sujetos <- list.files(path="input",pattern=".xlsx")
list <- map(files_sujetos, lectura_datos)

#################################################################################
##### Gráfico de resultados para el k-ésimo individuo
#################################################################################

Num_subject <- 1
lista_wider <- list[[Num_subject]] %>% 
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
  xlab("Número de funciones base")

#################################################################################
##### Se encuentra la estimación del parámetro lambda mediante gcv
#################################################################################

df_gcv <- map_dfr(seq(-5,5,by=0.01),~Select_lambda(x= x , Y_mat = t(lista_wider_2), lambda=.x)) 
lambda_optimo <- 10^df_gcv[which.min(df_gcv$gcv),"lambda"]


Predict_functions <- function(Num_subject, Replica, lambda_optimo){
    # list[[Num_subject]] %>%
  #   filter(REPETICION == 1) %>% 
  #   ggplot() +
  #   geom_line(aes(MEDICION, VALOR, group = KEY)) +
  #   facet_grid(ELECTRODO ~ CLASE) 
  
  #################################################################################
  ##### Gráfico de resultados para el k-ésimo individuo
  #################################################################################
  
  lista_wider <- list[[Num_subject]] %>% 
    filter(REPETICION == Replica) %>% 
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
  # df_num_basis <- map_dfr(5:20, ~SelectionNumBasis(x = 1:228, Y_mat = t(lista_wider_2), .x))
  # df_num_basis %>% 
  #   melt(id.vars = "Num_basis", variable.name = "Variable", value.name = "Value") %>% 
  #   ggplot() + 
  #   geom_line(aes(x = Num_basis, Value, linetype = Variable, color = Variable), size = 1.1)+
  #   scale_color_brewer(palette="Set1")+theme_bw()+ylab("Cuadrado medio del error total")+
  #   xlab("Número de funciones base")+ylim(0,5)
  
  #################################################################################
  ##### Se encuentra la estimación del parámetro lambda mediante gcv
  #################################################################################
  
  # df_gcv <- map_dfr(seq(-5,5,by=0.01),~Select_lambda(x= x , Y_mat = t(lista_wider_2), lambda=.x)) 
  # lambda_optimo <- 10^df_gcv[which.min(df_gcv$gcv),"lambda"]
  
  #################################################################################
  ##### Generación de curvas estimadas
  #################################################################################
  
  #smooth_fda(x, Y_mat = t(lista_wider_2), lambda = lambda_optimo) %>% dim
  fitted_fda_data <- fitted_fda(x, Y_mat = t(lista_wider_2), lambda = lambda_optimo)
  
  
  total_var_prop <- 0
  k <- 1
  while(total_var_prop < 0.85){
    pca_overall <- pca.fd(fitted_fda_data$fd, nharm = k)
    total_var_prop <- sum(pca_overall$varprop)
    k <- k + 1 
  }
  
  Spatfd <- SpatFD(data = fitted_fda_data,
         coords = coordenadas[, 1:3],
         basis="Bsplines",
         nbasis = 15,
         lambda=lambda_optimo,
         nharm = k)
  
  #scores_fda <- scores(Spatfd)
  
  #KS_scores(SFD, newcoords,model,name=NULL,fill.all=NULL)
  #variogram(scores_fda[[1]]~1, locations = coordenadas[, 1:3])
  
  Predict_scores <- KS_scores(SFD = Spatfd, 
            newcoords = grid,
            model = vgm("Exp"),
            name=NULL,fill.all=NULL)
  
  Predicted_spat_fda <- recons_fd(Predict_scores)
  
  return(list(Spatfd = Spatfd, 
       Predict_scores = Predict_scores,
       Predicted_spat_fda = Predicted_spat_fda))
}

Results_s_1 <- map(1:300, ~Predict_functions(1, .x, lambda_optimo = lambda_optimo))

Results_s_1[[1]]$Predict_scores %>% names

Results_s_1[[1]]$Predict_scores$empirical_variogram




Results_s_2 <- map(1:300, ~Predict_functions(2, .x, lambda_optimo = lambda_optimo))
Results_s_3 <- map(1:300, ~Predict_functions(1, .x, lambda_optimo = lambda_optimo))
Results_s_4 <- map(1:300, ~Predict_functions(1, .x, lambda_optimo = lambda_optimo))
Results_s_5 <- map(1:300, ~Predict_functions(1, .x, lambda_optimo = lambda_optimo))
Results_s_6 <- map(1:300, ~Predict_functions(1, .x, lambda_optimo = lambda_optimo))
Results_s_7 <- map(1:300, ~Predict_functions(1, .x, lambda_optimo = lambda_optimo))
Results_s_8 <- map(1:300, ~Predict_functions(8, .x, lambda_optimo = lambda_optimo))
Results_s_9 <- map(1:300, ~Predict_functions(9, .x, lambda_optimo = lambda_optimo))
Results_s_10 <- map(1:300, ~Predict_functions(10, .x, lambda_optimo = lambda_optimo))
Results_s_11 <- map(1:300, ~Predict_functions(11, .x, lambda_optimo = lambda_optimo))
Results_s_12 <- map(1:300, ~Predict_functions(12, .x, lambda_optimo = lambda_optimo))
Results_s_13 <- map(1:300, ~Predict_functions(13, .x, lambda_optimo = lambda_optimo))
Results_s_14 <- map(1:300, ~Predict_functions(14, .x, lambda_optimo = lambda_optimo))
Results_s_15 <- map(1:300, ~Predict_functions(15, .x, lambda_optimo = lambda_optimo))



map(Results_s_1, ~.x$Predict_scores$variogram_model) # Esta línea extrae los variogramas de los diferentes componentes principales
Results_s_1[[1]]$Predict_scores$variogram_model

Results_s_1[[1]] %>% names
Results_s_1[[1]]$Spatfd$fitted_fda_data
Results_s_1[[1]]$Predicted_spat_fda
df_to_plot <- eval.fd(evalarg = 1:228, Results_s_1[[1]]$Predicted_spat_fda)
df_to_plot_1 <- df_to_plot %>% as.data.frame() %>% mutate(x = 1:228) %>%
  pivot_longer(names_to = "Curva", values_to = "Valor", cols = -(x)) %>% 
  cbind(grid)


# df_to_plot_1$Var1 <- df_to_plot_1$Var1 + abs(min(df_to_plot_1$Var1))
# df_to_plot_1$Var2 <- df_to_plot_1$Var2 + abs(min(df_to_plot_1$Var2))
# df_to_plot_1$Var3 <- df_to_plot_1$Var3 + abs(min(df_to_plot_1$Var3))


df_to_plot_2 <- df_to_plot_1 %>% filter(Var3 == -8.48 & x == 228)

df_to_plot_2 %>% ggplot(aes(Var1, Var2)) + geom_raster(aes(fill = Valor)) +
  scale_fill_gradient(low = "yellow", high = "red")

  