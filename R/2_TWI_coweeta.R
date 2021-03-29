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
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_coweeta//"
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#Load DEM and pour points
dem<-raster(paste0(data_dir,"dem_10m"))
pp<-st_read(paste0(data_dir,"pp.shp"))

#Plot in mapview for funzies
mapview(dem) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Whole Basin Analysis----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Create stream -------------------------------------------------------------
#Export DEM to scratch workspace
writeRaster(dem, paste0(scratch_dir,"dem.tif"), overwrite=T)

#Define thresshold area to define stream
threshold=4046

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
  out_type= "cells",
  output = "fac.tif",
  wd = scratch_dir
)

#Create Stream Layer
wbt_extract_streams(
  flow_accum = "fac.tif",
  output = "stream.tif",
  threshold = threshold,
  wd = scratch_dir
)

#2.2 Estimate TWI --------------------------------------------------------------
#Create bianary stream layer
wbt_multiply(
  input1 = "stream.tif",
  input2 = "0",
  output = "stream_0.tif",
  wd= scratch_dir
)
wbt_add(
  input1 = 'stream_0.tif',
  input2 = "1",
  output = 'stream_1.tif',
  wd=scratch_dir
)
  
#Estimate slope
wbt_slope(
  dem ="dem_breached.tif",
  output = "slope.tif", 
  wd=scratch_dir
)

#estimate SCA
wbt_d8_flow_accumulation(
  input = "dem_breached.tif",
  output = "sca.tif",
  out_type="specific contributing area", 
  wd = scratch_dir
)

#Run TWI Function
wbt_wetness_index(
  sca    = "sca.tif",
  slope  = "slope.tif",
  output = "twi.tif",
  wd     = scratch_dir
)

#Estimate TWI for stream
wbt_multiply(
  input1 = "twi.tif",
  input2 = "stream_1.tif",
  output = "twi_stream.tif",
  wd     = scratch_dir
)

# 2.3 Estimate contributing area -----------------------------------------------
#Define area of grid cells
grid_area<-res(dem)[1]*res(dem)[2]

#Create catchment area grid
wbt_multiply(
  input1 = "fac.tif",
  input2 = grid_area, 
  output = "ca.tif",
  wd     = scratch_dir
)

#Create catchment area grid for stream
wbt_multiply(
  input1 = "ca.tif",
  input2 = "stream_1.tif",
  output = "ca_stream.tif",
  wd     = scratch_dir
)
  
#2.3 Convert stream grids to points----------------------------------------------
#Convert to polygon
wbt_raster_streams_to_vector(
  streams = "stream.tif",
  d8_pntr = "fdr.tif",
  output = "streams.shp",
  wd = scratch_dir
)

#Bring sream polygon into r environment
streams<-st_read(paste0(scratch_dir, "streams.shp"))
st_crs(streams)<-st_crs(dem@crs)

#Bring twi into R env
twi<-raster(paste0(scratch_dir,"twi.tif"))

#Convert to point
pnts<-
  st_line_sample(
      x=streams, 
      density = 0.01) %>% 
  st_coordinates(stream_pnt) %>% 
  as_tibble() %>% 
  st_as_sf(
    coords  = c('X','Y'), 
    crs = st_crs(dem@crs)) %>% 
  mutate(twi = extract(twi, .))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Watershed Delineation ------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Create function to create watershed shape~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fun<-function(n){
  
  #isolate pp
  pp<-pp[n,]
  
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
  stream[stream<2500]<-NA
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

#3.2 run function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheds<-lapply(seq(1,nrow(pp)), fun) %>% bind_rows()

#3.3 Add area data in ha ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sheds$ara_ha<-st_area(sheds, byid=T)
sheds$ara_ha<-as.numeric(paste0(sheds$ara_ha))/10000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 4: Export ----------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Save map file
setwd("docs/")
m<-mapview(sheds, map.types = c("Esri.WorldShadedRelief", "OpenStreetMap.DE")) + 
  mapview(streams) + 
  mapview(pnts, zcol = 'twi') 
mapshot(m, "twi_coweeta.html")

