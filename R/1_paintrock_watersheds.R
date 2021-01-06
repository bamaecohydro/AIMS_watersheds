#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Paint Rock Watershed Delineation
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
library(whitebox)
library(sf)
library(raster)
library(stars)
library(mapview)
library(htmlwidgets)

#Define data directories
data_dir<-"C://WorkspaceR//PaintRock//spatial_data//I_data//"
scratch_dir<-"C://WorkspaceR//PaintRock//spatial_data//II_scratch//"
output<-"C://WorkspaceR//PaintRock//spatial_data//III_output//"

#Load DEM and pour points
dem<-raster(paste0(data_dir,"dem_10m.tif"))
pp<-st_read(paste0(data_dir, "watershed_outlets.shp"))
gps<-st_read(paste0(data_dir, "202012_GPS_Tracks.shp"))


#Plot in mapview for funzies
mapview(dem) + mapview(pp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Delineate Watersheds --------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Create function to create watershed shape~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fun<-function(n){
  
  #isolate pp
  pp<-pp[n,]

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

  #Create Stream Layer
  stream<-raster(paste0(scratch_dir,"fac.tif"))
  stream[stream<10000]<-NA
  writeRaster(stream, paste0(scratch_dir,"stream.tif"), overwrite=T)

  #Paste point points in scratch dir
  st_write(pp, paste0(scratch_dir,"pp.shp"), delete_dsn = T)

  #Snap pour point
  wbt_jenson_snap_pour_points(
    pour_pts = "pp.shp", 
    streams = "stream.tif",
    snap_dist = 100,
    output =  "snap.shp",
    wd= scratch_dir)

  #2.4 Delineat watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  wbt_watershed(
    d8_pntr = "fdr.tif",
    pour_pts = "snap.shp", 
    output = "sheds.tif" ,
    wd=scratch_dir)

  #load watershed raster into R env
  sheds<-raster(paste0(scratch_dir,"sheds.tif"))
  
  #Convert raster to vector
  sheds<- sheds %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #Add pp ID
  sheds$name <-pp$Name
  
  #export shape
  sheds
}

#2.2 run function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheds<-lapply(seq(1,nrow(pp)), fun) %>% bind_rows()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Plot ------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#unload package
detach("package:whitebox")
detach("package:stars")

#Clean up shape 
sheds<-sheds %>% 
  #remove unneeded cols
  dplyr::select(-sheds.tif) %>% 
  #Calc area
  mutate(area_ha = st_area(.)) %>% 
  mutate(area_ha = as.numeric(area_ha)/10000) %>% 
  filter(area_ha>0.1) %>% 
  arrange(-area_ha)

#rname rows
rownames(sheds)<-sheds$name

#plot
m<-mapview(
  sheds,
  zcol="name",  
  alpha.regions=0.3,
  map.types=c("OpenTopoMap")) +
mapview(gps, color="red")

#SAve map file
setwd("docs/")
mapshot(m, "sheds.html")
