#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: STIC Location
#Date: 4/18/2022
#Coder: Nate Jones (cnjones7@ua.edu)
#Purpose: STIC Map for super synoptic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: Setup workspace -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear Memory
remove(list = ls())

#Load relevant packages
library(tidyverse) #join the cult
library(ggsn) # for scalebar
library(patchwork)
library(vegan)
library(whitebox)
library(sf)
library(raster)
library(stars)
library(mapview); mapviewOptions(fgb = FALSE)
library(htmlwidgets)
library(BAMMtools)

#Define data directories
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_tal_updated//"
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#Read in points
pnts<-st_read(paste0(data_dir, "allsites.shp"))
flow_net<-st_read(paste0(data_dir, "TAL_stream_network.shp"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: Create Map ------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m<-mapview(flow_net) + mapview(pnts)
mapshot(m, "docs/all_points.html")
