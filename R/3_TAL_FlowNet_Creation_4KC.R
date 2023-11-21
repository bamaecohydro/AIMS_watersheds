#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: TAL Stream Network Creation
#Date: 11/218/2023
#Coder: Nate Jones (cnjones7@ua.edu)
#Purpose: Create streamnetwork for KC Zarek
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_Talladega//"

#Define data inputs
dem<-raster(paste0(data_dir,"elevation//lidar01m33085g5.tif"))
crop<-st_read(paste0(data_dir, "crop.shp"))
dem<-crop(dem,crop)
pp<-tibble(
    lat = c(33.76218), 
    lon = c(-85.5955)) %>% 
  st_as_sf(., coords = c("lon", "lat"), crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(., crs = st_crs(dem))

#Define variables
threshold <- 25000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Create flow net -------- ----------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Delineate stream network --------------------------------------------------
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
  wd = scratch_dir)

#Flow accumulation raster
wbt_d8_flow_accumulation(
  input = "dem_breached.tif",
  out_type= "cells",
  output = "fac.tif",
  wd = scratch_dir)

#Create Stream Layer
wbt_extract_streams(
  flow_accum = "fac.tif",
  output = "stream.tif",
  threshold = threshold,
  wd = scratch_dir)

#Create flow net
wbt_raster_streams_to_vector(
  streams = "stream.tif",
  d8_pntr = "fdr.tif", 
  output = "flownet.shp",
  wd = scratch_dir
)

#Bring flownet into R environment
flownet <- st_read(paste0(scratch_dir, "flownet.shp"), crs = st_crs(dem))

#2.2 Watershed Delineation -----------------------------------------------------
#Paste point points in scratch dir
st_write(pp, paste0(scratch_dir,"pp.shp"), delete_dsn = T)

#Snap pour point
wbt_jenson_snap_pour_points(
  pour_pts = "pp.shp", 
  streams = "stream.tif",
  snap_dist = 100,
  output =  "snap.shp",
  wd= scratch_dir)

#Delineat watersheds
wbt_watershed(
  d8_pntr = "fdr.tif",
  pour_pts = "snap.shp", 
  output = "sheds.tif" ,
  wd=scratch_dir)

#load watershed raster into R env
shed<-raster(paste0(scratch_dir,"sheds.tif"))

#Convert raster to vector
shed<- shed %>% st_as_stars() %>% st_as_sf(., merge = TRUE)

#2.3 Clip flow net to watershed ------------------------------------------------
flownet <- flownet[shed,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Export to workspace ---------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plot for funzies
mapview(shed) + mapview(flownet)

#Export flownet
st_write(shed   , paste0(output_dir, "TAL_shed_KC.shp"))
st_write(flownet, paste0(output_dir, "TAL_flownet_KC.shp"))
