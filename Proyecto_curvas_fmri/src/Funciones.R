labels <- function(x)format(x, big.mark = ".", scientific = FALSE) 

lectura_datos <- function(.x){
  datos_sujeto <- read.xlsx(paste0("input/",files_sujetos[.x]))
  label_subject <- str_extract(files_sujetos[.x],"\\d+")
  
  datos_sujeto <- datos_sujeto %>% as.data.table()
  
  datos_sujeto[,c("MEDICION","REPETICION"):=list(rep(1:228,(nrow(datos_sujeto)/228)),
                                                 rep(1:(nrow(datos_sujeto)/228),each=228))]
  
  datos_sujeto_1 <- datos_sujeto %>%
    melt(id.vars = c("FRECUENCIA","CLASE","MEDICION","REPETICION"),variable.name="ELECTRODO",value.name = "VALOR")
  datos_sujeto_1[,c("SUJETO","KEY"):=list(paste0("S",label_subject),
                                                 paste0("S",label_subject,"_",REPETICION,"_",ELECTRODO))]
  #class<-datos_sujeto_1[,c("REPETICION","CLASE")] %>% unique() 
  #datos_sujeto_1 %>% ggplot()+geom_line(aes(MEDICION,VALOR,group=REPETICION,color=ELECTRODO))+facet_wrap(~CLASE,scales="free")
  return(list(df=datos_sujeto_1))
}


turn_it_to_matrix_to_fda <- function(.x){
  .x$df[,c("MEDICION","REPETICION","VALOR")] %>% 
    dcast(MEDICION~REPETICION,value.var="VALOR") %>% dplyr::select(-"MEDICION")
}
