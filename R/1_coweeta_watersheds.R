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
output<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#Load DEM and pour points
dem<-raster(paste0(data_dir,"dem.tif"))
sheds <- st_read(paste0(data_dir,"coweeta_subwatersheds.shp"))

#Plot in mapview for funzies
mapview(dem) + mapview(sheds)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Plot ------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#unload package
detach("package:whitebox")
detach("package:stars")

#Clean up shape 
sheds<-sheds %>% 
  #remove unneeded cols
  dplyr::select(WATERSHED) %>% 
  #Calc area
  mutate(area_ha = st_area(.)) %>% 
  mutate(area_ha = as.numeric(area_ha)/10000) %>% 
  arrange(-area_ha)

#plot
m<-mapview(
  sheds,
  zcol="WATERSHED",  
  alpha.regions=0.3,
  map.types=c("OpenTopoMap")) 

#SAve map file
setwd("docs/")
mapshot(m, "coweeta_sheds.html")
