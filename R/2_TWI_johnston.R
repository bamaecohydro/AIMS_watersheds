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
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_johnston_draw//"
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#Load DEM and pour points
dem<-raster(paste0(data_dir,"rcew_DEM_3m_filled.tif"))
pp<-tibble(
    x = 518285, 
    y = 4774255 ) %>% 
  st_as_sf(., coords = c("x","y"), crs = st_crs(dem@crs))
  
#Plot in mapview for funzies
mapview(dem) + mapview(pp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Spatial Analysis --------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.1 Create function to create watershed shape~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fun<-function(pp, dem, threshold=1000, name){
  
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

  #Paste point points in scratch dir
  st_write(pp, paste0(scratch_dir,"pp.shp"), delete_dsn = T)

  #Snap pour point
  wbt_jenson_snap_pour_points(
    pour_pts = "pp.shp", 
    streams = "stream.tif",
    snap_dist = 100,
    output =  "snap.shp",
    wd= scratch_dir)

  #2.4 Delineate watersheds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
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
  sheds$name <-name
  
  #2.5 Estimate TWI and Contributing Area ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  wbt_multiply(
    input1 = "twi_stream.tif",
    input2 = "sheds.tif",
    output = "twi_stream.tif",
    wd     = scratch_dir
  )
  
  #Estiamte SCA for stream
  wbt_multiply(
    input1 = "sca.tif",
    input2 = "stream_1.tif",
    output = "sca_stream.tif",
    wd     = scratch_dir
  )
  wbt_multiply(
    input1 = "sca_stream.tif",
    input2 = "sheds.tif",
    output = "sca_stream.tif",
    wd     = scratch_dir
  )
  
  # #Export to gloval env
  # twi_stream<-raster(paste0(scratch_dir,"twi_stream.tif"))
  # return(twi_stream)
  
  #Convtert twi_stream to points
  wbt_raster_to_vector_points(
    input  = "twi_stream.tif",
    output = "twi_pnts.shp",
    wd     = scratch_dir
  )
  
  #Add fac value
  grid_area<-res(dem)[1]*res(dem)[2]
  wbt_multiply(
    input1 = "fac.tif",
    input2 = grid_area, 
    output = "ca.tif",
    wd     = scratch_dir
  )
  wbt_multiply(
    input1 = "ca.tif",
    input2 = "stream_1.tif",
    output = "ca_stream.tif",
    wd     = scratch_dir
  )
  wbt_multiply(
    input1 = "ca_stream.tif",
    input2 = "sheds.tif",
    output = "ca_stream.tif",
    wd     = scratch_dir
  )
  wbt_raster_to_vector_points(
    input  = "ca_stream.tif",
    output = "ca_pnts.shp",
    wd     = scratch_dir
  )
  
  #Join spatially
  twi<-st_read(paste0(scratch_dir,"twi_pnts.shp"))
  sca<-st_read(paste0(scratch_dir,"ca_pnts.shp"))
  pnts<-st_join(twi, sca) %>% 
    rename(twi = VALUE.x,
           con_area_m2 = VALUE.y) %>% 
    mutate(con_area_ha = con_area_m2/10000) %>% 
    dplyr::select(twi, twi, con_area_ha)
  export<- pnts %>% 
    st_drop_geometry() 
  
  #Reduce DEM to shed
  wbt_multiply(
    input1 = "dem_breached.tif",
    input2 = "sheds.tif", 
    output = "dem_clip.tif",
    wd     = scratch_dir
  )
  dem_clip<-raster(paste0(scratch_dir,"dem_clip.tif"))
  
  #2.6 export points
  return(list(sheds, pnts, export, dem_clip))
}

#2.2 run function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output<-fun(pp, dem, name="dry_creek")
shed<-output[[1]]
pnts<-output[[2]]
export<-output[[3]]
dem_clip<-output[[4]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Plot ------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Create map ----------------------------------------------------------------
area_max<-max(pnts$con_area_ha, na.rm=T)
pnts<-pnts %>% mutate(area_prop =con_area_ha/area_max)
st_crs(pnts)<-st_crs(shed)

#plot
m<-mapview(
    shed,
    alpha.regions=0.3,
    map.types=c("OpenTopoMap")) +
  #mapview(dem) +
  mapview(pnts, zcol='twi') +
  mapview(pnts, zcol='area_prop')

m

#3.2 Create Plot ---------------------------------------------------------------
p<-export %>% 
  as_tibble() %>% 
  ggplot(aes(y=twi, x = con_area_ha)) +
    geom_point(pch=19, col="grey30", alpha=70) + 
    theme_bw() + 
      xlab("Contributing Area [ha]") +
      ylab("TWI") + 
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) 
p

#3.3 Export --------------------------------------------------------------------
#Save map file
setwd("docs/")
mapshot(m, "twi_johnston.html")

png("twi_jonston.png", height = 3, width = 3.5, units = "in", res=100)
p
dev.off()
