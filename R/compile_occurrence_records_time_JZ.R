#Valerio Barbarossa
# stript that compiles the occurrence records from the different datasets

library(valerioUtils)
library(lubridate)
libinv(c('dplyr','vroom'))

# diagnostics file
file_diag <- 'filtering_occurrence_datasets_diag.log'

diag <- function(df,name_df = '',file_out = '',append = T){
  
  if(file_out == ''){
    cat('Dataset: ',name_df,'\n',
        'No. records: ',prettyNum(nrow(df),big.mark = ','),'\n',
        'No. species: ',prettyNum(length(unique(df$name)),big.mark = ','),'\n\n\n')
  }else{
    cat('Dataset: ',name_df,'\n',
        'No. records: ',prettyNum(nrow(df),big.mark = ','),'\n',
        'No. species: ',prettyNum(length(unique(df$name)),big.mark = ','),'\n\n\n',
        file = file_out,append = append)
  }
}

# read reference names that should be used to extract species from the datasets----------------------------
# tab <- read.csv('proc/names_fishbase.csv',stringsAsFactors = F)
fishbase <- read.csv('proc/names_fishbase_and_tedesco_validated20220513.csv',stringsAsFactors = F)


# functions to clean nomenclature--------------------------------------------------------------------------

# function that removes any trailing word with capital (e.g., Linneaus)
rm_author <- function(x){
  s <- strsplit(x,' ')
  if(length(s[[1]]) > 2){
    s3 <- strsplit(s[[1]][3],'')
    if(s3[[1]][1] %in% LETTERS){
      return(
        paste(unlist(s)[1:2],collapse=' ')
      )
    }else{
      return(x)
    }
  }else{
    return(x)
  }
}

# function that removes anything after , or ( and trailing capitals
clean_binomial <- function(names){
  cl <- lapply(strsplit(names,','),function(x) x[[1]][1]) %>% 
    unlist() %>%
    strsplit(.,'\\(') %>%
    lapply(.,function(x) x[[1]][1]) %>%
    trimws() %>%
    lapply(.,rm_author) %>%
    unlist
  return(cl)
}

# clean and filter datasets--------------------------------------------------------------------------------
# setup location of datasets
dir_data <- '../data/fish_databases/'

# ala.org.au
ala <- vroom(paste0(dir_data,'ala.org.au/Fishes-brief.csv'),delim = ',') %>%
  select(name = scientificName,lon = decimalLongitude,lat = decimalLatitude, Date_occ = eventDate) %>% # JZ-add a column of time
  .[!is.na(.$lon) & !is.na(.$name) & !is.na(.$Date_occ),] %>%  # JZ-remove those without time information
  mutate(year = year(Date_occ))%>%  # JZ-extract year information
  mutate(name = clean_binomial(name)) %>%
  filter(name %in% filter_names) %>%
  filter(year>=1970 & year <=2010) %>% # JZ-filter 1990-2000
  #filter(year>=1900) %>% # JZ-filter 1990-most recent
  select(name, lon, lat, year) %>% # JZ-only keep spatial and temporal information
  distinct()
diag(ala,'ala.org.au')
diag(ala,'ala.org.au',file_diag)

#boldsystems #not implemented

#fishnet2
fishnet <- vroom(paste0(dir_data,'fishnet2/fishnet2.csv')) %>%
  select(name = ScientificName,lon = Longitude,lat = Latitude, year = YearCollected) %>%
  .[!is.na(.$lon) & !is.na(.$name) & !is.na(.$year),] %>%
  mutate(name = clean_binomial(name)) %>%
  filter(name %in% filter_names) %>%
  #filter(year>=1900 & year <=2000) %>% # JZ-filter 1990-2000
  filter(year>=1970 & year <=2010) %>% # JZ-filter 1990-most recent
  distinct()

diag(fishnet,'fishnet2')
diag(fishnet,'fishnet2',file_diag)


#gbif
gbif_o <- vroom(paste0(dir_data,'gbif/actinopterygii.csv'))

# based on species
gbifs <- gbif_o %>%
  select(name = species,lon = decimalLongitude,lat = decimalLatitude, year = year) %>% # JZ-add a column of time
  .[!is.na(.$lon) & !is.na(.$name) & !is.na(.$year),] %>% # JZ-remove those without time information
  mutate(name = clean_binomial(name)) %>%
  filter(name %in% filter_names) %>%
  #filter(year>=1900 & year <=2000)  # JZ-filter 1990-2000
  filter(year>=1970 & year <=2010)  # JZ-filter 1990-most recent

# based on scientific name
gbifn <- gbif_o %>%
  select(name = scientificName,lon = decimalLongitude,lat = decimalLatitude, year = year) %>% # JZ-add a column of time
  .[!is.na(.$lon) & !is.na(.$name) & !is.na(.$year),] %>% # JZ-remove those without time information
  mutate(name = clean_binomial(name)) %>%
  filter(name %in% filter_names) %>%
  #filter(year>=1900 & year <=2000) # JZ-filter 1990-2000
  filter(year>=1970 & year <=2010)  # JZ-filter 1990-most recent

gbif <- bind_rows(gbifs,gbifn) %>%
  distinct()

# checknames <- gbif_o %>%
#   mutate(species = clean_binomial(species),scientificName = clean_binomial(scientificName))
# ss <- sort(table(checknames$species),decreasing = T)
# sn <- sort(table(checknames$scientificName),decreasing = T)

diag(gbif,'gbif')
diag(gbif,'gbif',file_diag)
rm(gbif_o) # unload gbif_o from memory


#portalbiodiversidade.icmbio.gov.br
bra <- vroom(paste0(dir_data,'portalbiodiversidade.icmbio.gov.br/portalbio_export_17-09-2019-10-21-11.csv')) %>%
  select(name = Especie,lon = Longitude,lat = Latitude, Date_occ = "Data do evento") %>% # JZ-add a column of time
  .[.$lon != "Acesso Restrito" & .$name != "Sem Informações" & .$Date_occ != "Sem Informações" & !is.na(.$Date_occ),] %>%  # JZ-remove those without time information
  mutate(date_oc=as.Date(substr(Date_occ,1,10),tryFormats = c("%d/%m/%Y")))%>% # JZ-extract date information
  mutate(year = year(date_oc)) %>% # JZ-extract year information
  mutate(name = clean_binomial(name)) %>%
  filter(name %in% filter_names) %>%
  #filter(year>=1900 & year <=2000) %>% # JZ-filter 1990-2000
  filter(year>=1970 & year <=2010)%>%  # JZ-filter 1990-most recent
  select(name, lon, lat, year) %>% # JZ-only keep spatial and temporal information
  distinct()
diag(bra,'portalbiodiversidade.icmbio.gov.br')
diag(bra,'portalbiodiversidade.icmbio.gov.br',file_diag)


#splink.org
splink <- vroom(paste0(dir_data,'splink.org/speciesLink_all_112728_20190917101630.txt')) %>%
  select(name = scientificname,lon = longitude,lat = latitude, year = yearcollected) %>% # JZ-remove those without time information
  .[!is.na(.$lon) & !is.na(.$name) & !is.na(.$year),] %>% # JZ exclude null year
  mutate(name = clean_binomial(name)) %>%
  filter(name %in% filter_names) %>%
  #filter(year>=1900 & year <=2000) %>% # JZ-filter 1990-2000
  filter(year>=1970 & year <=2010)%>%  # JZ-filter 1990-most recent
  select(name, lon, lat, year) %>% # JZ-only keep spatial and temporal information
  distinct()
diag(splink,'splink.org')
diag(splink,'splink.org',file_diag)

# combine records from different datasets------------------------------------------------------------------

# bind and filter occurrence records
occ <- rbind(ala,fishnet,gbif,bra,splink) %>%
  arrange(name) %>%
  distinct()

# there are still some lat and lon coordinates with letters in there e.g., 1044730E
# difficult to convert without knowing the CRS, filter out!
# create vector to report rows with letters in latitude or longitude

# understand unique symbols to take outsbatch 
unique_symbols <- strsplit(paste0(occ$lon,occ$lat),'') %>%
  do.call('c',.) %>%
  unique(.)

# exclude everything that is not numeric (numbers,.,-)
# checkout biogeo package also for more through corrections
to_exclude <- unique_symbols[!unique_symbols %in% c(as.character(0:9),'.','-')] 
lonlat_to_exclude <- lapply(strsplit(paste0(occ$lon,occ$lat),''),
                            function(x) sum(x %in% to_exclude)> 0 ) %>%
  do.call('c',.)

# and filter occ
cat('Removing ',prettyNum(sum(lonlat_to_exclude),big.mark = ','),' records with letters in coordinates\n',
    file = file_diag,append = T)

occ <- occ %>%
  filter(!lonlat_to_exclude) %>%
  mutate(name,lon = as.numeric(lon),lat = as.numeric(lat)) %>% # transform lat lon to numeric
  filter(!is.na(lon) & !is.na(lat)) %>% #make sure there are no NAs in the coords
  filter(lon >= -180 & lon <= 180) %>% # and no coords outside allowed boundaries
  filter(lat >= -90 & lat <= 90)


# merge with fishbase synonyms----------------------------------------------------------------------------

# read fishbase + tedesco synonym tab
syn_fb <- fishbase %>%
  as_tibble() %>%
  select(name_ref = name,name = name_synonym)

occ_fb_nosyn <- occ %>%
  filter(name %in% unique(syn_fb$name_ref))
occ_fb_syn <- inner_join(occ,syn_fb,by='name') %>% #select only records that are synonyms
  select(name = name_ref,lon,lat,year) %>% # and assign the fishbase name to the name
  distinct()

# ### count species that have synonyms only
# syn_fb2 <- syn_fb %>%
#   filter(name != name_ref) 
# occ_fb_syn <- occ %>%
#   filter(name %in% unique(syn_fb2$name)) %>% # and assign the fishbase name to the name
#   distinct()
# length(unique(occ_total_fb$name))

occ_total_fb <- bind_rows(occ_fb_nosyn,occ_fb_syn) %>%
  arrange(name) %>%
  distinct()

diag(occ_fb_nosyn,'merged occurrence records (no synonyms)',file_diag)
diag(occ_fb_syn,'additional merged occurrence records from fishbase synonyms',file_diag)
diag(occ_total_fb,'final occurrence records cleaned and checked for synonyms',file_diag)

vroom_write(occ_total_fb,'proc/compiled_occurrence_records_JZ.tsv.gz')
saveRDS(occ_total_fb,'proc/compiled_occurrence_records_JZ.rds')

# Jz-species counting
no_species <- occ_total_fb %>%
  select(name=name) %>%
  distinct()




# JZ-records that do not consider the year
no_loc_occ <- occ %>% select(name, lon, lat) %>% distinct()
# # Jz-species/occirrences counting by year
# as.data.frame(table(occ_total$year)) #occurence by year
# species_oc <- occ_total_fb %>%
#   select(name=name,year) %>%
#   distinct()
# as.data.frame(table(species_oc$year)) #species by year
#JZ-round the logitude/latitude to decrease the data item for map showing
# library(raster)
# # Calculate record of 1900-most recent 2017
occ_total_rd_19702010 <- occ_total_fb  %>%
  select(name,lon,lat,year) %>%
  mutate(lon_rd = round(2*lon)/2, lat_rd = round(2*lat)/2) %>%
  count(.$lon_rd,.$lat_rd)
head_rst <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, res=0.5, crs="+proj=longlat +datum=WGS84")

occ_total_rst_19702010 <- rasterize(occ_total_rd_19702010[, c('.$lon_rd', '.$lat_rd')], head_rst, occ_total_rd_19702010[, 'n'])
plot(occ_total_rst_19702010)

writeRaster(occ_total_rst_19702010, 'map/occurence_species1970-2010.asc', format = 'ascii',overwrite=TRUE) # write raster
# JZ- remove unecessory variables
rm(ala,bra,fishnet,
   gbif,gbifn,gbifs)
rm(occ_total, occ_syn,occ_total_rd_19002000,
   occ_total_rd_19002017,occ_total_rd_19702000,occ_total_rd_19702010,
   occ_total_rd_19702017,occ_total_rst_19002000,occ_total_rst_19002017,
   occ_total_rst_19702000,occ_total_rst_19702017)