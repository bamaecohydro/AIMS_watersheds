#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: STIC Placement
#Date: 7/5/2021
#Coder: Nate Jones (cnjones7@ua.edu)
#Purpose: Define initial stic placement
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
scratch_dir<-"C://WorkspaceR//AIMS_watersheds//data//II_scratch//"
output_dir<-"C://WorkspaceR//AIMS_watersheds//data//III_output//"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 2: STIC Function ---------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Initiate function ---------------------------------------------------------
stic_fun<-function(
  dem, #Digital elevation model for site
  pp,  #watershed outlet location
  threshold, #Rough threshold to initiate stream channel
  n_stics = 20, #number of stics
  buffer_dist = 25, #dist between stics
  hydro_pnts, #approximate location of hydrologic infrastructure
  stic_clip, #areas where stics cannot be locatated
  scratch_dir #Scratch directory
  ){

#2.1 Delineate Watershed -------------------------------------------------------
#Clear scratch dir
sapply(list.files(scratch_dir, full.names = T, recursive = T), unlink)

#Set random seed
set.seed(10)
  
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

#2.2 Stream network delineation ------------------------------------------------
#Create Stream Layer
wbt_extract_streams(
  flow_accum = "fac.tif",
  output = "stream.tif",
  threshold = threshold,
  wd = scratch_dir
)

#Convert  streams raster to polygon
wbt_raster_streams_to_vector(
  streams = "stream.tif",
  d8_pntr = "fdr.tif",
  output = "streams.shp",
  wd = scratch_dir
)

#Bring stream polygon into r environment
streams<-st_read(paste0(scratch_dir, "streams.shp"))
st_crs(streams)<-st_crs(dem@crs)

#Crop streams to sheds
streams<-streams[sheds,]

#2.3 Characterize points along stream ------------------------------------------
#distance to outlet
wbt_distance_to_outlet(
  d8_pntr = 'fdr.tif',
  streams = "stream.tif",
  output = 'dst2out.tif',
  wd = scratch_dir
)

#estimate SCA
wbt_d8_flow_accumulation(
  input = "dem_breached.tif",
  output = "sca.tif",
  out_type="specific contributing area", 
  wd = scratch_dir
)

#Run slope function
wbt_slope(
  dem ="dem_breached.tif",
  output = "slope.tif", 
  wd=scratch_dir
)

#Run TWI Function
wbt_wetness_index(
  sca    = "sca.tif",
  slope  = "slope.tif",
  output = "twi.tif",
  wd     = scratch_dir
)

#Bring rasters into R env
twi<-raster(paste0(scratch_dir,"twi.tif"))
fac<-raster(paste0(scratch_dir,"fac.tif"))
dist2out<-raster(paste0(scratch_dir, "dst2out.tif"))

#Convert stream network to point w/ TWI and Aws information
pnts<-
  st_line_sample(
    x=streams, 
    density = res(dem)[1]) %>%  
  st_coordinates() %>% 
  as_tibble() %>% 
  sf::st_as_sf(
    coords  = c('X','Y'), 
    crs = st_crs(dem@crs)) %>% 
  dplyr::rename(StreamReach = L1) %>% 
  mutate(
    dist2out = raster::extract(dist2out, .),
    twi = raster::extract(twi, .), 
    con_area_ha = raster::extract(fac, .)/10000,
    pid = seq(1, nrow(.)))  # a unique identifier for each point

#Crop to stream
pnts<-pnts[sheds,]

#2.4 Subset points into bins ---------------------------------------------------
#Create 10 'bins' based drainage area
n_groups_area <- 10
bins_area <- getJenksBreaks(pnts$con_area_ha, n_groups_area+1)
pnts$bin_1<-cut(pnts$con_area_ha, breaks=bins_area) %>% as.numeric(.)

#create function to subset bins based on TWI
sub_bin_fun<-function(n_bin){
  
  # calculate number of sites to sample
  sites_per_bin <- round(n_stics/n_groups_area,0)
 
  #Define sub bins
  temp<-pnts %>% 
    st_drop_geometry() %>% 
    filter(bin_1 == n_bin) %>% 
    mutate(
      bin_2 = cut(con_area_ha, breaks=sites_per_bin), 
      bin_2 = as.numeric(bin_2)) %>% 
    dplyr::select(pid, bin_2)
  
  #export temp tibble
  temp
}

#execute function
sub_bin<-lapply(seq(1, n_groups_area), sub_bin_fun) %>% bind_rows()

#Join to pnts
pnts<-left_join(pnts, sub_bin) 

#Identify unique groups
group<-pnts %>% 
  st_drop_geometry() %>% 
  dplyr::select(bin_1, bin_2) %>% 
  unique() %>% drop_na() %>% 
  mutate(group = seq(1, nrow(.)))
pnts<-left_join(pnts, group); remove(group)

#2.5 Create function to remove points along stream -----------------------------
trim_fun<-function(pids){
  
  #create inner function to identify points to crop based on river distance
  inner_fun<-function(n){
    #Isolate points based on pids
    temp<-pnts %>% filter(pid==pids[n]) %>% 
      mutate(dist2out = if_else(is.na(dist2out), 0, dist2out))
    
    #Identify points within buffer dist along stream
    output<-pnts %>% 
      st_drop_geometry() %>% 
      filter(StreamReach == temp$StreamReach) %>% 
      filter(dist2out <= temp$dist2out+buffer_dist,
             dist2out >= temp$dist2out - buffer_dist) %>% 
      select(pid) %>% 
      pull()
    
    #Export
    output
  }
  
  #Apply inner funciton and export!
  lapply(seq(1, length(pids)), inner_fun) %>% unlist
}

# 2.6 Remove hydro points ------------------------------------------------------
#create function to identify closest point
snap_fun<-function(n){

  #Identify feature
  output<-hydro_pnts %>% 
    slice(n) %>% 
    st_buffer(., dist=100) %>% 
    st_intersection(pnts, .) %>% 
    mutate(dist = st_distance(., hydro_pnts %>% slice(n))) %>% 
    slice_min(dist) %>% 
    st_drop_geometry() %>% 
    select(pid) %>% 
    pull()
  
  #Export output
  output
}

#Identify points to remove
hydro_pnt_pids<-lapply(seq(1, nrow(hydro_pnts)), snap_fun) %>% unlist()

#remove from points
crop<-trim_fun(hydro_pnt_pids)
pnts_cropped<-pnts %>% 
  dplyr::filter(!(pid %in% crop)) 

#2.7 Remove areas where we can't place stics -----------------------------------
stic_clip<-pnts_cropped %>% 
  st_intersection(., stic_clip) %>% 
  st_drop_geometry() %>% 
  select(pid) %>% pull()

pnts_cropped<-pnts_cropped %>% 
  filter(!(pid %in% stic_clip))

#2.8 STIC Placement ------------------------------------------------------------
#Randomly select point from first group
stics<-pnts_cropped %>% 
  st_drop_geometry() %>% 
  arrange(dist2out) %>% 
  filter(group == 1) %>% 
  sample_n(1) 

#Remove points within buffer length along stream reach
crop<-trim_fun(stics$pid)
pnts_cropped<-pnts_cropped %>% dplyr::filter(!(pid %in% crop))

#identify remaining groups
grp<-pnts_cropped %>% 
  st_drop_geometry() %>% 
  select(group) %>% 
  distinct() %>% 
  drop_na() %>% 
  pull() 
  
#Iterate through all remaining groups
i<-2
while(i < length(grp)){
  temp<-pnts_cropped %>% 
    st_drop_geometry() %>% 
    arrange(dist2out) %>% 
    filter(group == grp[i]) %>% 
    sample_n(1)
  crop<-trim_fun(temp$pid)
  pnts_cropped<-pnts_cropped %>% dplyr::filter(!(pid %in% crop))
  stics<-bind_rows(stics, temp)
  pnts_cropped %>% 
    st_drop_geometry() %>% 
    select(group) %>% 
    distinct() %>% 
    drop_na() %>% 
    pull() 
  i<-i+1
  print(i)
}

#filter points to stics
stics<-stics %>% select(pid) %>% pull()
stics<-pnts %>% filter(pid %in% stics)

#2.9 Export STICS --------------------------------------------------------------
list(stics, streams, sheds)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Watershed Analysis ----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Shambly --------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_weyerhaeuser//"

#Load DEM
dem<-raster(paste0(data_dir,"dem_crop.tif"))

#Load watershed outlet location
pp<-tibble(
    y=c(32.98465), 
    x=c(-88.01226)) %>% 
  st_as_sf(., 
    coords = c("x", "y"), 
    crs = '+proj=longlat +datum=WGS84 +no_defs') %>% 
  st_transform(., crs = st_crs(dem@crs)) 

#Load existing hydro infrastructure
flumes<-st_read(
    paste0(data_dir,"flumes.shp"), 
    crs =st_crs(" +proj=utm +zone=16 +datum=WGS84 +units=m +no_defs")) %>% 
  rename(site = Flumes) %>% 
  filter(site == "GR1" | site == "GR2") 
ltm<-st_read(paste0(output_dir, "ltm_shambly.shp")) %>% 
  mutate(site = paste0('shambly_',seq(1,nrow(.))))
pp <- pp %>% mutate(site='super_sensor')
hydro_pnts<-bind_rows(flumes, ltm, pp)
  
#Load no-go zone shapefiles for stics
stic_clip<-st_read(paste0(data_dir,"stic_clip.shp"))

#run function
output<-stic_fun(dem, pp, threshold = 50000, 
  n_stics = 20, buffer_dist = 25,
  hydro_pnts, stic_clip, scratch_dir 
)

#Create interactive plots
m<-mapview(output[[1]]) + 
  mapview(output[[2]]) + 
  mapview(hydro_pnts,  col.regions = 'orange')
m
mapshot(m, "docs/stics_shambly.html")

#Export Points
st_write(output[[1]], paste0(output_dir, "stics_shambly.shp"), append=FALSE)
st_write(output[[2]], paste0(output_dir, "streams_shambly.shp"), append = F)
st_write(output[[3]], paste0(output_dir, "shed_shambly.shp"), append = F)

#3.2 Paint Rock --------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_paintrock//"

#Load DEM
dem<-raster(paste0(data_dir,"dem_mask_10m.tif"))

#Load watershed outlet location
pp<-st_read(paste0(data_dir, "watershed_outlets.shp")) %>% 
  filter(Name == "Fanning Hollow") %>% 
  st_transform(., crs = st_crs(dem@crs))

#Load existing hydro infrastructure
ltm<-st_read(paste0(output_dir, "ltm_paintrock.shp")) %>% 
  mutate(site = paste0('paintrock_',seq(1,nrow(.))))
pp <- pp %>% mutate(site='super_sensor')
hydro_pnts<-bind_rows(ltm, pp)

#Load no-go zone shapefiles for stics
stic_clip<-st_read(paste0(data_dir,"stic_clip.shp"))

#run function
output<-stic_fun(dem, pp, threshold = 1000, 
                 n_stics = 20, buffer_dist = 25,
                 hydro_pnts, stic_clip, scratch_dir 
)

#Create interactive plots
m<-mapview(output[[1]]) + 
  mapview(output[[2]]) + 
  mapview(hydro_pnts,  col.regions = 'orange')
m
mapshot(m, "docs/stics_paintrock.html")

#Export Points
st_write(output[[1]], paste0(output_dir, "stics_paintrock.shp"), append = F)
st_write(output[[2]], paste0(output_dir, "streams_paintrock.shp"), append = F)
st_write(output[[3]], paste0(output_dir, "shed_paintrock.shp"), append = F)

#3.3 Talladega --------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_Talladega//"

#Load DEM
dem<-raster(paste0(data_dir,"elevation//lidar01m33085g5.tif"))
  crop<-st_read(paste0(data_dir, "crop2.shp"))
  dem<-crop(dem,crop)

#Load watershed outlet location
pp<-st_read(paste0(data_dir,"pp.shp"))
  
#Load existing hydro infrastructure
ltm<-st_read(paste0(output_dir, "ltm_talladega.shp")) %>% 
  mutate(site = paste0('talladega_',seq(1,nrow(.))))
pp <- pp %>% mutate(site='super_sensor')
hydro_pnts<-bind_rows(ltm, pp)

#run function
output<-stic_fun(dem, pp, threshold = 35000, 
                 n_stics = 20, buffer_dist = 25,
                 hydro_pnts, scratch_dir 
)

#Create interactive plots
m<-mapview(output[[1]]) + 
  mapview(output[[2]]) + 
  mapview(hydro_pnts,  col.regions = 'orange')
m
mapshot(m, "docs/stics_talladega.html")

#Export Points
st_write(output[[1]], paste0(output_dir, "stics_talladega.shp"), append = F)
st_write(output[[2]], paste0(output_dir, "streams_talladega.shp"), append = F)
st_write(output[[3]], paste0(output_dir, "shed_talladega.shp"), append = F)

