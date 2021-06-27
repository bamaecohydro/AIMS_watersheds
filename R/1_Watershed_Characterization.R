#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Watershed Delineation
#Date: 6/25/2021
#Coder: Nate Jones (cnjones7@ua.edu)
#Purpose: Characterize Watershed Attributes
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
library(mapview); mapviewOptions(fgb = FALSE)

#Define data directories
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Create analysis function ----------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Function
fun<-function(
  scratch_dir, #Temp folder to store intermediates
  dem,         #Input digital elevation model
  pp,          #Watershed pourpoint
  threshold    #Threshold for stream length
  ){

  #2.1 Initial Raster Processing -------------------------------------------------
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
  
  #2.2 Watershed Delineation -----------------------------------------------------
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
  sheds<-raster(paste0(scratch_dir,"sheds.tif"))
  
  #Convert raster to vector
  sheds<- sheds %>% st_as_stars() %>% st_as_sf(., merge = TRUE)
  
  #Define area [ha]
  sheds$area_ha<-st_area(sheds, byid=T)
  sheds$area_ha<-as.numeric(paste0(sheds$area_ha))/10000
  sheds$area_ha<-round(sheds$area_ha, 0)
  
  #2.3 Estimate characteristics along catchment-----------------------------------
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
    
  #Convert  streams rater to polygon
  wbt_raster_streams_to_vector(
    streams = "stream.tif",
    d8_pntr = "fdr.tif",
    output = "streams.shp",
    wd = scratch_dir
  )
  
  #Bring stream polygon into r environment
  streams<-st_read(paste0(scratch_dir, "streams.shp"))
  st_crs(streams)<-st_crs(dem@crs)
  
  #Bring twi, fac, and slope into R env
  twi<-raster(paste0(scratch_dir,"twi.tif"))
  twi<-crop(twi, sheds)
  fac<-raster(paste0(scratch_dir,"fac.tif"))
  fac<-crop(fac, sheds)
  slope<-raster(paste0(scratch_dir,"slope.tif"))
  slope<-crop(slope, sheds)
  dem<-raster(paste0(scratch_dir,"dem_smoothed.tif"))
  dem<-crop(dem, sheds)
  
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
    mutate(
      twi = round(raster::extract(twi, .), 0), 
      #con_area_ha = round(raster::extract(fac, .)/10000, 0),
      elevation_m = round(raster::extract(dem, .), 0),
      slope = round(raster::extract(slope, .),0))

  #Clip to watershed
  pnts<-pnts[sheds,] 

  #2.4 Create and export output ------------------------------------------------
  #Create list
  output<-list(sheds, pnts, streams)
  
  #Export list
  output
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Southeast Region Watersheds--------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Coweeta -------------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_coweeta//"

#Define data inputs
dem<-raster(paste0(data_dir,"dem_10m"))
pp<-st_read(paste0(data_dir,"pp.shp"))[3,]
threshold=4046

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

# Export Map
mapshot(m, "docs//Coweeta.html")
st_write(shed, paste0(output_dir, "coweeta_shed.shp"))

#3.2 Paint Rock (Fanning Hollow) -----------------------------------------------
#Define Data Directory
data_dir<-data_dir<-"C://WorkspaceR//PaintRock//spatial_data//I_data//"

#Define data inputs
dem<-raster(paste0(data_dir,"dem_10m.tif"))
pp<-st_read(paste0(data_dir, "watershed_outlets.shp")) %>% 
  filter(Name=="Fanning Hollow")
threshold=1000

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//paintrock.html")
st_write(shed, paste0(output_dir, "paintrock_shed.shp"))

#3.3 Weyerhaeuser (Shambley Creek) -----------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_weyerhaeuser//"

#Define data inputs
dem<-raster(paste0(data_dir,"dem.tif"))
pp<-tibble(
  y=c(32.98465), 
  x=c(-88.01226)) %>% 
  st_as_sf(., 
           coords = c("x", "y"), 
           crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(., crs = st_crs(dem@crs))
threshold<-20000

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//weyerhaeuser.html")
st_write(shed, paste0(output_dir, "weyerhaeuser.shp"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 4: Great Plains Region Watersheds--------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.1 Blue Creek (Oklahoma) -----------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_blue//"

#Define data inputs
p<-"+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
dem<-raster(paste0(data_dir,"dem.tif"))
dem<-projectRaster(dem, crs=p)
pp<-tibble(
  x = 34.455846,
  y = -96.664498) %>% 
  st_as_sf(., coords = c("y","x"), crs=4326) %>% 
  st_transform(., crs = st_crs(p))
threshold<-4604

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//blue.html", selfcontained=T)
st_write(shed, paste0(output_dir, "blue.shp"))

#4.2 King Creek (Konza) -----------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_kanza//"

#Define data inputs
p<-"+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
dem<-raster(paste0(data_dir,"GIS201.tif"))
dem<-projectRaster(dem, crs=p)
pp<-tibble(
  x = 39.093064,
  y = -96.588494) %>% 
  st_as_sf(., coords = c("y","x"), crs=4326) %>% 
  st_transform(., crs = st_crs(p))
threshold<-4604

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//konza.html", selfcontained=T)
st_write(shed, paste0(output_dir, "konza.shp"), delete_layer = T)

#4.3 Youngmeyer Ranch ----------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_youngmeyer//"

#Define data inputs
p<-"+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
dem<-raster(paste0(data_dir,"dem.tif"))
dem<-projectRaster(dem, crs=p)
pp<-tibble(
  x = 37.532459,
  y = -96.471819) %>% 
  st_as_sf(., coords = c("y","x"), crs=4326) %>% 
  st_transform(., crs = st_crs(p))
threshold<-4604

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//young.html", selfcontained=T)
st_write(shed, paste0(output_dir, "young.shp"), delete_layer = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 5: Mountain West Region Watersheds--------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#5.1 Dry Creek -----------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_dry_creek//"

#Define data inputs
#Load DEM and pour points
dem<-raster(paste0(data_dir,"DCEW-DEMclip.tif"))
pp<-st_read(paste0(data_dir, "pp.shp"))
threshold<-1000

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//drycreek.html", selfcontained=T)
st_write(shed, paste0(output_dir, "drycreek.shp"), delete_layer = T)

#5.2 Gibson Jack ---------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_gibson//"

#Define data inputs
#Load DEM and pour points
dem<-raster(paste0(data_dir,"gibson_lidar_1m_UTM_dem.tif"))
pp<-st_read(paste0(data_dir, "gage.shp")) %>% 
  st_transform(., st_crs(dem@crs))
threshold<-10000

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//gibson.html", selfcontained=T)
st_write(shed, paste0(output_dir, "gibson.shp"), delete_layer = T)

#5.3 Johnston Draw -------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_johnston_draw//"

#Define data inputs
#Load DEM and pour points
dem<-raster(paste0(data_dir,"rcew_DEM_3m_filled.tif"))
pp<-tibble(
  x = 518285, 
  y = 4774255 ) %>% 
  st_as_sf(., coords = c("x","y"), crs = st_crs(dem@crs))
threshold<-1000

#Run model
output<-fun(scratch_dir, dem, pp, threshold)
shed<-output[[1]]
pnts<-output[[2]]
streams<-output[[3]]

#Create interactive map
m<-mapview(
  shed,
  alpha.regions=0.3) +
  mapview(streams) +
  mapview(pnts, zcol='twi') 
m

#export
mapshot(m, "docs//johnston.html", selfcontained=T)
st_write(shed, paste0(output_dir, "johnston.shp"), delete_layer = T)

