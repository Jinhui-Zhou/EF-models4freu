# Jinhui Zhou
# Script for assign fish species occurrence data to pairwise(N concentration) data
library(raster)
library(rasterize)
library(dplyr)
library(conflicted)
library(sf)
library(terra)
library(rgdal)
# library(sp)
conflict_prefer('filter','dplyr')
conflict_prefer('select','dplyr')
conflict_scout()

# data path
# setwd('C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/04data/fish_spicies/occ2range4fish-master/')
dir_fishdata = '../occ2range4fish-master/proc/'
dir_GNMdata = '../GNM/Nconc_valerio'
dir_pair = '../result/pairwise/'

# read GNM files
N_data_path <- list.files(path=dir_GNMdata,pattern ="asc", full.names = TRUE,recursive =TRUE)
Year_Nconc <- c(dir(path=dir_GNMdata))
stack_N_conc <- stack(N_data_path)
crs(stack_N_conc) <- "EPSG:102100"
names(stack_N_conc) <- Year_Nconc # name the stack with respective years

# counting species richness
# Calculate record of 1970-most recent 2010
sperich_total <- occ_total_fb  %>%
  select(name,lon,lat,year) %>%
  mutate(lon_rd = round(lon*2)/2+0.25, lat_rd = round(lat*2)/2+0.25) %>%
  filter(year>=1970 & year <=2010) %>%
  distinct() 

# compile pairwise data
pair_total <- c()
for (i in Year_Nconc){ 
  year_GNM <- as.numeric(i)
  sperich_year <- sperich_total %>%
    filter(year==year_GNM)
  # location of species occurrence: to search the N concentration in the same cell
  loc_occ <- cbind(sperich_year$lon_rd,sperich_year$lat_rd)
  # extract N concentration
  N_conc <- c()
  Layer_GNM <- paste('X',year_GNM,sep='')
  extract_reseult<- raster::extract(stack_N_conc[[Layer_GNM]], loc_occ)
  N_conc <- as.list(extract_reseult)
  # Add a column to save N concentration
  pair_year <- sperich_year %>%
    mutate(N_conc)
  pair_total <- rbind(pair_total,pair_year)
}

names(pair_total)[names(pair_total)=='N_conc']='Nconc'
pair_total_NoNA <- pair_total %>%
  filter(!is.na(pair_total$Nconc))
#names(pair_total_NoNA) <- c('Longitude','Latitude','Year','Species_richness','Nconc')

# pair_write <- data.frame(c(unlist(pair_total_NoNA$Longitude)),
#                          c(unlist(pair_total_NoNA$Latitude)),
#                          c(unlist(pair_total_NoNA$Year)),
#                          c(unlist(pair_total_NoNA$Species_richness)),
#                          c(unlist(pair_total_NoNA$Nconc)))
# names(pair_write) <- c('Longitude','Latitude','Year','Species_richness','Nconc')

# transform pairwise data to points
# the points are too many to process, so we split it into every 10000 items
# pair_point1 <- st_as_sf(x = pair_total_NoNA[1:10000,], 
#                         coords = c("Longitude", "Latitude"),
#                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# loc_pnt1 <- pair_total_NoNA[1:10000,] %>%
#   select(Longitude,Latitude)

#coordinates(pair_total_NoNA) <- c("Longitude", "Latitude")
#proj4string(pair_total_NoNA) <- proj4string(4326)
# assign pairwise data to shapefile-ecoregion
# path of ecoregion maps
dir_map <- 'C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/05figure/01map'
path_ecoregion <- paste(dir_map,'/02ecoregion/rst_ecocopy.asc',sep='')
# read ecoregion map
rst_ecoregion <- raster(path_ecoregion)
crs(rst_ecoregion) <- "EPSG:102100"


# path_ecoregion <- paste(dir_map,'/02ecoregion/feow_hydrosheds-original.shp',sep='')
# read ecoregion map
# shp_ecoregion <- read_sf(path_ecoregion) %>%
#   st_as_sf() %>%
#   st_set_crs(3857) 
# # transfer shapefile to raster
# r = raster(shp_ecoregion,res=0.5)
# rst_ecoregion = rasterize(shp_ecoregion,r,'FEOW_ID',fun=min)
# crs(rst_ecoregion) <- "EPSG:102100"
# rst_ecoregion<-resample(rst_ecoregion,stack_N_conc)
# rst_ecoregion <- round(rst_ecoregion)
### distinguish those small ecoregions with less than 10 cells
cellsize_ecoregion = as.data.frame(freq(rst_ecoregion))
lst_ecoregion <- cellsize_ecoregion %>%
  filter(.$count>=10 & !is.na(.$value)) %>%
  select(FEOW_ID=value)


# location of species occurrence: to search the N concentration in the same cell
loc_occ <- cbind(pair_total_NoNA$lon_rd,pair_total_NoNA$lat_rd)
# extract ecoregion ID
ECO_FEOW_ID <- raster::extract(rst_ecoregion, loc_occ) %>%
  as.list()

# Add a column to ecoregion ID
pair_total_NoNA <- pair_total_NoNA %>%
  mutate(ECO_FEOW_ID) %>%
  filter(ECO_FEOW_ID %in% lst_ecoregion$FEOW_ID)
names(pair_total_NoNA)[names(pair_total_NoNA)=='ECO_FEOW_ID']='FEOW_ID'

# save the result of species-Nconcentration-Ecoregion as csv
pair_write <- data.frame(c(unlist(pair_total_NoNA$name)),
                         c(unlist(pair_total_NoNA$lon)),
                         c(unlist(pair_total_NoNA$lat)),
                         c(unlist(pair_total_NoNA$year)),
                         c(unlist(pair_total_NoNA$lon_rd)),
                         c(unlist(pair_total_NoNA$lat_rd)),
                         c(unlist(pair_total_NoNA$Nconc)),
                         c(unlist(pair_total_NoNA$FEOW_ID)))
names(pair_write) <- c('Species_name','Longitude','Latitude','Year','Lon_raster','Lat_raster','Nconc','FEOW_ID')
write.csv(pair_write,file = 'C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/05figure/01map/pair_write20220824.csv',na="NA",row.names=TRUE)
#### maximum N conc for species-----ecoregion--USE PYTHON---
# merge realm/major habitat info in pair data
table_realm_hab <- read.csv( 'C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/05figure/01map/02ecoregion/feow_list_join.csv')
pair_biome <- merge(pair_write, table_realm_hab, by= 'FEOW_ID',all = FALSE)
write.csv(pair_biome,file = 'C:/Users/DELL/Desktop/zhouj8/01Research/04Global_CF/05figure/01map/pair_wise_biome20220824.csv',na="NA",row.names=TRUE)

##### maximum N conc for species-----realm-habitat-ecoregion
species_survival <- function(pair_data,region='FEOW_ID'){
  pair_maxconc <- aggregate(pair_biome$Nconc, by = list(pair_biome$Species_name,pair_biome[[region]]), max)
  names(pair_maxconc) <- c('Species_name','Nconc',region)
  pair_sperich_conc <- aggregate(pair_maxconc$Species_name, by = list(pair_maxconc$Nconc,pair_maxconc[[region]]), length) %>%
    distinct()
  names(pair_sperich_conc) <- c(region,'Nconc','Species_maxcon')  
  # numbers of survival species at N concentration level
  num_survival <- pair_sperich_conc %>%
    group_by(pair_sperich_conc[[region]])%>%
    arrange(desc(Nconc)) %>%
    summarise(new=cumsum(Species_maxcon),.groups='drop') %>%
    select(new)

  data_survival <- pair_sperich_conc %>%
    arrange(pair_sperich_conc[[region]], desc(Nconc)) %>%
    mutate(Species_richness = num_survival)
  
  # print the regionalized stats for data pair
  print(as.data.frame(table(data_survival[[region]])))
  # count how many species can survive by using the maximal N concentration
  pair_write2 <- data.frame(c(unlist(data_survival[[region]])),
                            c(unlist(data_survival$Nconc)),
                            c(unlist(data_survival$Species_maxcon)),
                            c(unlist(data_survival$Species_richness)))
  names(pair_write2) <- c('ID','Nconc','Species_maxcon','Species_richness')
  return(pair_write2)
}
# read species name against N conc with realm and habitat information
pair_biome <- read.csv(paste(dir_map,'/pair_wise_biome20220824.csv',sep=''))
pair_biome_survival <- species_survival(pair_biome,region='Biome')
write.csv(pair_biome_survival,file = paste(dir_map,'/pair_survival_biome20220824.csv',sep=''),na="NA")

pair_ecoregion_survival <- species_survival(pair_biome,region='Ecoregion')
write.csv(pair_ecoregion_survival,file = paste(dir_map,'/pair_survival_Ecoregion20220824.csv',sep=''),na="NA")

pair_realm_survival <- species_survival(pair_biome,region='Realm')
write.csv(pair_realm_survival,file = paste(dir_map,'/pair_survival_Realm20220824.csv',sep=''),na="NA")

pair_habitat_survival <- species_survival(pair_biome,region='Major.Habitat.Type')
write.csv(pair_habitat_survival,file = paste(dir_map,'/pair_survival_Habitat20220824.csv',sep=''),na="NA")


# clean data
rm(pair_point1,pair.point,pair_eco_data,pair_point,pair_total,
   sperich_year,sperich_total,pair_year,pair_write,r,
   occ_total_rst,ECO_FEOW_ID,loc_occ,loc_pnt,N_conc,pts,pair_survival)
gc()




