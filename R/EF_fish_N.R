library(raster)
library(rgdal)
library(dplyr)

setwd("C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/04data/fish_spicies/occ2range4fish-master/")
# 
# data path
dir_GNMdata = '../GNM/Nconc_valerio'
# dir_result = '../result/pairwise/'
dir_map <- 'C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/05figure/01map'

# read raster file
rd_conc_fun <- function(dir,year){
  str_year <- as.character(year)
  val_raster <- raster(N_data_path[str_year])
  crs(val_raster) <- "EPSG:102100" 
  val_raster
}
# zonal statistics
zonal_stats_fun <- function(val_raster, region_raster){
  val_raster <- setExtent(val_raster, region_raster, keepres=TRUE)
  zonal_stats <- zonal(val_raster, region_raster,fun=mean)
  zonal_stats
}
# effect factor functions- for rasters
# ef_average_fun <- function(conc_current, PDF_current,conc_0,PDF_0){
#   # average EF
#   diff_conc <- overlay(conc_current,
#                     conc_0,
#                     fun=function(r1, r2){return(r1-r2)})
#   diff_pdf <- overlay(PDF_current,
#                    PDF_0,
#                    fun=function(r1, r2){return(r1-r2)})
#   ef_average_raster <- overlay(diff_pdf,
#                             diff_conc,
#                             fun=function(r1, r2){return(r1/r2)})
#   ef_average_raster
# }


ef_fun <- function(df_current,df_target,df_criteria){
  # merge current and target state of Nconc, and the cretiria, a, b
  df_average <- merge(df_current, df_target, by = c("x","y","FEOW_ID"), all= TRUE) %>%
    merge(.,df_criteria, by = "FEOW_ID", all= TRUE) %>%
    rename(Nconc.current=Nconc.x, Nconc.target = Nconc.y) %>%
    na.omit()
  # predict PDF of current and target state
  ssd_current <- expression(1 / (1 + exp((a - log10(Nconc.current))/b)))
  ssd_target <- expression(1 / (1 + exp((a - log10(Nconc.target))/b)))
  df_average  <- within(df_average,PDF.current <-eval(ssd_current))
  df_average <- within(df_average,PDF.target <-eval(ssd_target))  
  df_average["PDF.current"][df_average["Nconc.current"]==0] <- NA
  df_average["PDF.target"][df_average["Nconc.target"]==0] <- NA
  # predict average EF
  df_average <- within(df_average, 
                       EF_average <- 1000*(PDF.current-PDF.target)/(Nconc.current-Nconc.target))

  # predict marginal EF 
  ssd_deriv <- D(expression(1 / (1 + exp((a - log10(Nconc.current))/b))), "Nconc.current")
  df_average <- within(df_average, EF_marginal <- 1000*eval(ssd_deriv)) %>%
    select(x,y,FEOW_ID,Nconc.current,Nconc.target,PDF.current,PDF.target,Pseudo_r2,a,b,EF_average,EF_marginal) %>%
    filter(Pseudo_r2>=0.5)
  # df_average["EF_average"][is.na(df_average["EF_average"])] <- 0.0
  df_average
}

# ef_marginal_fun <- function(df_current,df_criteria){
#   # marginal EF
#   df_current <- merge(df_current, df_criteria, by = "FEOW_ID", all= TRUE)
#   ssd_deriv <- D(expression(1 / (1 + exp((a - Nconc) * b))), "Nconc")
#   df_current <- within(df_current, EF_marginal <- eval(ssd_deriv)) %>%
#     select(x,y,FEOW_ID,Nconc,Pseudo_r2,a,b,EF_marginal) %>%
#     filter(Pseudo_r2>=0.5)
#   df_current}


ef_to_raster_fun <- function(df_current){
  # dataframe of EF TO RASTER
  EF_raster <- df_current %>%
    select(x,y,EF) %>%
    rasterFromXYZ(.)
  crs(EF_raster) <- "EPSG:102100"
  EF_raster
}


# Nconc_current_ecoregion <- rd_conc_fun(N_data_path,2010) %>%
#   zonal_stats_fun(.,rst_ecoregion)
# margin_effect_fun <- function(df_current,df_criteria_ecoregion,df_criteria_biome){
# 
#   ### EF marginal
#   ef_marginal_ecoregion <- ef_marginal_fun(df_current,df_criteria_ecoregion) %>%
#     select(x,y,FEOW_ID,EF_marginal) %>%
#     rename(EF_marginal = EF_marginal)%>%
#     mutate(source = 1) %>%
#     na.omit()
#   # ef_marginal <- ef_marginal_ecoregion
#   ef_marginal_biome <- ef_marginal_fun(df_current,df_criteria_biome) %>%
#     select(x,y,FEOW_ID,EF_marginal) %>%
#     mutate(source = 2) %>%
#     rename(EF_marginal = EF_marginal) 
#   
#   # ef_marginal <- merge(ef_marginal_ecoregion, ef_marginal_biome, by= 'loc',all=TRUE) %>%
#   #   select(x.x, y.x, FEOW_ID.x, FEOW_ID.y, EF_marginal_ecoregion, EF_marginal_biome) %>%
#   #   rename(x = x.x, y = y.x)
#   ef_marginal <- rbind(ef_marginal_ecoregion,ef_marginal_biome)%>%
#     distinct(x, y, .keep_all =TRUE) %>%
#     na.omit()
#     
# 
#   # ef_marginal <- ef_marginal %>%
#   #   mutate(EF_marginal=NaN) %>%
#   #   within(.,EF_marginal <- ifelse(is.na(EF_marginal_ecoregion), 
#   #                                    EF_marginal_biome, 
#   #                                    EF_marginal_ecoregion))
#   
#   ###
#   EF_marginal_raster <- ef_marginal %>%
#     select(x,y,EF_marginal) %>%
#     rename(EF=EF_marginal)  %>%
#     ef_to_raster_fun()
#   
#   
#   EF_marginal_raster
# }


main_effect_fun <- function(df_current,df_target,df_criteria_ecoregion,df_criteria_biome){
  ### EF 
  ef_ecoregion <- ef_fun(df_current,df_target, df_criteria_ecoregion) 
  ef_biome <- ef_fun(df_current,df_target,df_criteria_biome) 
  
  ### EF with ecoregion SSD
  ef_ecoregion <-  ef_ecoregion %>%
    select(x,y,FEOW_ID,EF_average,EF_marginal,PDF.current,PDF.target) %>%
    mutate(source = 1) %>%
    na.omit()

  # EF with biome SSD
  ef_biome <-  ef_biome %>%
    select(x,y,FEOW_ID,EF_average,EF_marginal,PDF.current,PDF.target) %>%
    mutate(source = 2) %>%
  na.omit()
  # combine EF, ecoregion the first 
  ef_combine <- rbind(ef_ecoregion,ef_biome)%>%
    distinct(x, y, .keep_all =TRUE) %>%
  na.omit()

  # Set EF average as no value when PDF current > PDF target
  ef_combine['EF_average'][ef_combine["PDF.current"] < ef_combine["PDF.target"]] <- NA
  
  # write.csv(ef_combine,file = paste(out_fig,'/efcombine_parameter.csv',sep=''),na="NA")

  ################### average EF #########
  ### dataframe to raster
  EF_average_raster <- ef_combine %>%
    select(x,y,EF_average) %>%
    rename(EF=EF_average)  %>%
    ef_to_raster_fun()
  

  writeRaster(EF_average_raster, paste(dir_map,'/02ecoregion/EF_marginal/EF_average_2010_no0kg.tif',sep=''),
              overwrite=TRUE)
  ################### marginal EF #########
  ### dataframe to raster
  EF_marginal_raster <- ef_combine %>%
    select(x,y,EF_marginal) %>%
    rename(EF=EF_marginal)  %>%
    ef_to_raster_fun()
  writeRaster(EF_marginal_raster, paste(dir_map,'/02ecoregion/EF_marginal/EF_marginal_2010_no0kg.tif',sep=''),
              overwrite=TRUE)
  ################### marginal EF #########
  ### dataframe to raster
  PDF_current_raster <- ef_combine %>%
    select(x,y,PDF.current) %>%
    rename(EF=PDF.current)  %>%
    ef_to_raster_fun()
  writeRaster(PDF_current_raster, paste(dir_map,'/02ecoregion/EF_marginal/PDF_current_2010_no0.tif',sep=''),
              overwrite=TRUE)

}

# read raster of ecoregion ID
path_ecoregion <- paste(dir_map,'/02ecoregion/rst_ecocopy.asc',sep='')
rst_ecoregion <- raster(path_ecoregion)
crs(rst_ecoregion) <- "EPSG:102100"

###read cretiria (pseudo r2) and a/b for fitting curve of ecoregion/biomes
table_name <- read.csv( paste(dir_map,'/02ecoregion/feow_list_join_biome.csv',sep=''))
df_criteria_ecoregion <- merge(criteria_ecoregion_N, table_name, by.x= 'ID', by.y= 'Ecoregion',all.x)
df_criteria_biome <- merge(criteria_biome_N, table_name, by.x= 'ID', by.y= 'Biome',all.x)

### extract N concentration
# directory of N concentration
N_data_path <- list.files(path=dir_GNMdata,pattern ="asc", full.names = TRUE,recursive =TRUE)
Year_Nconc <- c(dir(path=dir_GNMdata))
names(N_data_path) <- Year_Nconc
# unify the rasters of Nconc and ecoregion ID
Nconc_current_ecoregion <- rd_conc_fun(N_data_path,2010)
Nconc_target_ecoregion <- rd_conc_fun(N_data_path,1900)

### raster to dataframe 
# current state
df_current <-  stack(list(Nconc_current_ecoregion,rst_ecoregion)) %>%
  raster::as.data.frame(.,xy=TRUE) %>%
  na.omit()
names(df_current) <- c("x","y","Nconc", "FEOW_ID")
df_current["Nconc"][df_current["Nconc"]<0.0001] <- 0.0

# desied target--1900
df_target <-  stack(list(Nconc_target_ecoregion,rst_ecoregion)) %>%
  raster::as.data.frame(.,xy=TRUE) %>%
  na.omit()
names(df_target) <- c("x","y","Nconc", "FEOW_ID")
df_target["Nconc"][df_target["Nconc"]<0.0001] <- 0.0

# calculate  EF and PDF_current
main_effect_fun(df_current,df_target,df_criteria_ecoregion,df_criteria_biome)


# effect factors ====

