#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Weyerhaeuser Watershed Delineation
#Date: 12/21/2020
#Coder: Nate Jones (cnjones7@ua.edu)
#Purpose: Initial watershed delineation for #AIMS project
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Required WBT Package -- See download directions here: https://github.com/giswqs/whiteboxR

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Setup workspace -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Clear Memory
remove(list = ls())

#Load relevant packages
library(tidyverse) #join the cult
library(readxl)
library(whitebox)
library(sf)
library(raster)
library(stars)
library(mapview)
library(htmlwidgets)

#Define data directories
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_weyerhaeuser//"
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#Define master CRS
p<-"+proj=utm +zone=16 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#Load relevant data
outlets<-st_read(paste0(data_dir,"flumes.shp"))
dem<-raster(paste0(data_dir,"dem.tif"))

#Plot in mapview for funzies
mapview(outlets, map.types=c("OpenTopoMap", "Esri.WorldImagery")) + mapview(dem)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Delineate Watersheds --------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Prep DEM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Export DEM to scratch workspace
writeRaster(dem, paste0(scratch_dir,"dem.tif"), overwrite=T)

#Smooth DEM
wbt_gaussian_filter(
  input = "dem.tif", 
  output = "dem_smoothed.tif",
  wd = scratch_dir)

#breach depressions
wbt_breach_depressions(
  dem =    "dem_smoothed.tif",
  output = "dem_breached.tif",
  fill_pits = F,
  wd = scratch_dir)

#2.2 Flow accumulation and direction rasters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Flow direction raster
wbt_d8_pointer(
  dem= "dem_breached.tif",
  output ="fdr.tif",
  wd = scratch_dir
)

#Flow accumulation raster
wbt_d8_flow_accumulation(
  input = "dem_breached.tif",
  output = "fac.tif",
  wd = scratch_dir
)

#2.3 Snap pour point ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Create Stream Layer
stream<-raster(paste0(scratch_dir,"fac.tif"))
stream[stream<10000]<-NA
writeRaster(stream, paste0(scratch_dir,"stream.tif"), overwrite=T)

#Write outlet points to the scratch dir
st_write(outlets, paste0(scratch_dir,"outlet.shp"), delete_dsn = T)

#Snap pour point
wbt_jenson_snap_pour_points(
  pour_pts = "outlet.shp", 
  streams = "stream.tif",
  output = "snapped.shp",
  snap_dist = 100,
  wd = scratch_dir)

#2.4 Delineat watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wbt_unnest_basins(
  d8_pntr = "fdr.tif",
  pour_pts = "snapped.shp", 
  output = "sheds.tif" ,
  wd=scratch_dir)

#load watershed and convert to polygon
sheds<-list.files(scratch_dir) %>% 
  as_tibble() %>% 
  dplyr::filter(stringr::str_detect(value,"tif")) %>% 
  dplyr::filter(stringr::str_detect(value, 'sheds'))
sheds<-lapply(
  FUN = function(x){
    r<-raster(paste0(scratch_dir, sheds[x,]))
    p<-r %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
    p}, 
  X = seq(1,nrow(sheds))) %>% bind_rows()

#2.5 Create stream vector ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stream<-raster(paste0(scratch_dir,"fac.tif"))
stream[stream<10000]<-NA
stream_10m<-raster::aggregate(stream, fact =10, fun=max, na.rm=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Plot ------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clean up shape 
sheds<-sheds %>% 
  #remove unneeded cols
  dplyr::select(-sheds_1.tif, -sheds_2.tif, -sheds_3.tif) %>% 
  #Calc area
  mutate(area_ha = st_area(.)) %>% 
  mutate(area_ha = as.numeric(area_ha)/10000) %>% 
  mutate(area_ha = round(area_ha)) %>% 
  filter(area_ha>0.1) %>% 
  arrange(-area_ha)

#plot
m<-mapview(
  sheds,
  #zcol="name",  
  alpha.regions=0.3,
  map.types=c("OpenTopoMap", "Esri.WorldImagery")) + 
mapview(dem) +
mapview(outlets) +
mapview(stream_10m)

#SAve map file
setwd("docs/")
mapshot(m, "cp_sheds.html")
