#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: STIC Placement
#Date: 4/25/2021
#Coder: Nate Jones (cnjones7@ua.edu)
#Purpose: Create function to estimate initiaal distribution of STICs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_weyerhaeuser//"
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#Load DEM and pour points
dem<-raster(paste0(data_dir,"dem.tif"))
pp<-st_read(paste0(data_dir,"flumes.shp")) %>% 
  filter(Flumes=="Alt 3") %>% slice(n=1)

#Plot in mapview for funzies
mapview(dem) + mapview(pp)

#Add stock var
name<- 'Weyerhauser'
threshold<-10000
n_stics<-25

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Setup workspace ------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Semi code: 

#Delineate watershed
#Create stream netework
#create distribution of TWI and Aws
#Export point, stream network, and watershed
#Export as google earth points


# -------------------------------------------------------------------------------
#fun<-function(pp, dem, threshold=22500, name, n_stics = 20){
  
  #A. Delineate watershed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  
  #Estimate slope
  wbt_slope(
    dem ="dem_breached.tif",
    output = "slope.tif", 
    wd=scratch_dir
  )
  
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
  
  #Estimate distance to outlet
  wbt_distance_to_outlet(
    d8_pntr = 'fdr.tif',
    streams = "stream.tif",
    output =  'dist.tif', 
    wd =       scratch_dir
  )
  
  #Paste point points in scratch dir
  st_write(pp, paste0(scratch_dir,"pp.shp"), delete_dsn = T)
  
  #Snap pour point
  wbt_jenson_snap_pour_points(
    pour_pts = "pp.shp", 
    streams = "stream.tif",
    snap_dist = 100,
    output =  "snap.shp",
    wd= scratch_dir)
  
  #delineate watershed 
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
  sheds$name <- "name"
  
  #B. Stream Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Convert to polygon
  wbt_raster_streams_to_vector(
    streams = "stream.tif",
    d8_pntr = "fdr.tif",
    output = "streams.shp",
    wd = scratch_dir
  )
  
  #Bring stream polygon into r environment
  streams<-st_read(paste0(scratch_dir, "streams.shp"))
  st_crs(streams)<-st_crs(dem@crs)
  
  #Crop sheds to basin
  streams<-streams[sheds,]
  
  #Bring twi, fac, and slope into R env
  twi<-raster(paste0(scratch_dir,"twi.tif"))
  twi<-crop(twi, sheds)
  fac<-raster(paste0(scratch_dir,"fac.tif"))
  fac<-crop(fac, sheds)
  slope<-raster(paste0(scratch_dir,"slope.tif"))
  slope<-crop(slope, sheds)
  dem<-raster(paste0(scratch_dir,"dem_smoothed.tif"))
  dem<-crop(dem, sheds)
  dist<-raster(paste0(scratch_dir,"dist.tif"))
  dist<-crop(dist, sheds)
  
  #EStimate location of STIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Convert stream network to point w/ TWI and Aws
  pnts<-
    st_line_sample(
      x=streams, 
      density = 1) %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    st_as_sf(
      coords  = c('X','Y'), 
      crs = st_crs(dem@crs)) %>% 
    mutate(
      twi = extract(twi, .), 
      con_area_ha = extract(fac, .)/10000,
      elevation_m = extract(dem, .),
      slope = extract(slope, .), 
      dist = extract(dist, .))
  
  #Define 'scaled rank' based on twi and Aws
  pnts<-pnts %>% 
    mutate(scale_twi = scale(twi), 
           scale_dist = scale(1/dist), 
           scale_sum = scale_twi + scale_dist) %>% 
    arrange(scale_sum) %>% 
    mutate(scale_rank = seq(1, nrow(.)))
  
  #Select n points from accross distribution
  rank<-seq(1,nrow(pnts), length.out = n_stics + 1)
  rank<-rank[2:(n_stics + 1)] %>% round( 0)
  pnts<-pnts %>% filter(scale_rank %in% rank)
    
  #Plot for testing 
  streams %>% st_geometry() %>% plot(col="#232D4B")
  pnts %>% st_geometry() %>% plot(pch=21, col="grey30", bg="#E57200", add=T, cex = 2)
  
  #Temporary output
  st_write(streams, "C:\\WorkspaceR\\AIMS_watersheds\\data\\I_data_weyerhaeuser\\map\\streams.shp", 
           delete_layer=T)
  
  st_write(pnts, "C:\\WorkspaceR\\AIMS_watersheds\\data\\I_data_weyerhaeuser\\map\\stics.shp", 
           delete_layer=T)
  
  
  
  
  # --------------------------------------------------------------------------------
  # Next Steps 
  #    Remove points within 200 m of each other. One approach -- seelct one point
  #    from each cluster of colocated points, and delete the others in the cluster. 
  #    Then, using all of the sampling points, create 200 m buffer around each point, clip 
  #    from the master list. Then, select the next closes point thats remaining for points lost.
  #   
  #    Finsih function
  #    Export as kml/kmz
  # --------------------------------------------------------------------------------
  
  
  # ---------------------
  
  #2.6 export points
  #return(list(sheds, pnts, streams))
#}

#2.2 run function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# output<-fun(pp, dem, threshold=20000, name="Kanza")
# shed<-output[[1]]
# pnts<-output[[2]]
# streams<-output[[3]]
# export<-pnts %>% st_drop_geometry()

