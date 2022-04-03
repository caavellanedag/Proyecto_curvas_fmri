#
# Autor(es): Camilo Avellaneda
# Mantenimiento: Camilo Avellaneda
# Fecha creación: 22/04/2021
#==============================================

if(!("input" %in% list.files())){dir.create("input")}
if(!("output" %in% list.files())){dir.create("output")}
if(!("src" %in% list.files())){dir.create("src")}
if(!("docs" %in% list.files())){dir.create("docs")}

require(pacman)
p_load(tidyverse,openxlsx,here,haven,data.table,fda,fda.usc,roahd,fdaoutlier)
p_load(maptools,sgeostat,scatterplot3d,car,fields,gstat,geoR)
 # cat("labels=function(x)format(x, big.mark = \".\", scientific = FALSE) ",
 #   file = "src/Funciones.R")
source("src/Funciones.R")
dir <- here::here()
coordenadas <- fread("input/coor1.txt")

coordenadas <- coordenadas %>%
  cbind(data.frame(ELECTRODO = paste0("E",1:21))) %>% 
  mutate(z=as.numeric(as.character(str_replace_all(z,",","\\."))))

files_sujetos <- list.files(path="input",pattern=".xlsx")

list <- map(1:21,lectura_datos)

matrix_to_fda <- map(list, turn_it_to_matrix_to_fda)

spline_basis <- create.bspline.basis(rangeval = c(1, 228), nbasis = 30)
fitted_basis <- smooth.basis(y = as.matrix(matrix_to_fda[[1]]), fdParobj = spline_basis)
matrix_to_fda[[1]] %>% dim
W.pca <- pca.fd(fitted_basis$fd, nharm=2)
plot(W.pca$harmonics, lwd=3)
sum(W.pca$varprop)
names(W.pca)
dim(W.pca$scores)
df_scores <- data.frame(W.pca$scores)
list[[1]]$df$REPETICION %>% table()
table(list[[2]]$df$SUJETO)
map_dfr(list,~.x$df)

as.geodata(df_scores,coords.col=4:6,data.col = 2)
coordenadas
sum(df_scores$x)
sum(df_scores$z)
