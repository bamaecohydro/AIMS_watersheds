#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Long term monitoring station placement
#Date: 7/10/2021
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
#Step 2: Sampling Location Function Function -----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampling_station_fun<-function(
  dem, #Digital elevation model for site
  pp,  #watershed outlet location
  threshold, #Rough threshold to initiate stream channel
  scratch_dir, #Scratch directory
  main_stem_stn = c(0.1,0.25,0.5, 0.75)
){

#2.1 Delineate Watershed -------------------------------------------------------
#Clear scratch dir
sapply(list.files(scratch_dir, full.names = T, recursive = T), unlink)
  
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
  output = "streams.tif",
  threshold = threshold,
  wd = scratch_dir
)

#Paste point points in scratch dir
st_write(pp, paste0(scratch_dir,"pp.shp"), delete_dsn = T)

#Snap pour point
wbt_jenson_snap_pour_points(
  pour_pts = "pp.shp", 
  streams = "streams.tif",
  snap_dist = 100,
  output =  "snap.shp",
  wd= scratch_dir)

#Delineat watersheds
wbt_watershed(
  d8_pntr = "fdr.tif",
  pour_pts = "snap.shp", 
  output = "shed.tif" ,
  wd=scratch_dir)

#load watershed raster into R env
shed<-raster(paste0(scratch_dir,"shed.tif"))

#Convert raster to vector
shed<- shed %>% st_as_stars() %>% st_as_sf(., merge = TRUE)

#Define area [ha]
shed$area_ha<-st_area(shed, byid=T)
shed$area_ha<-as.numeric(paste0(shed$area_ha))/10000
shed$area_ha<-round(shed$area_ha, 0)

#Write shed to workspace
st_write(shed, paste0(scratch_dir, "shed.shp"), delete_dsn=T)

#2.2 Stream network characterization -------------------------------------------
#2.2.1 Create Streams Vector file ----------------------------------------------
#Convert  streams raster to polygon
wbt_raster_streams_to_vector(
  streams = "streams.tif",
  d8_pntr = "fdr.tif",
  output = "streams.shp",
  wd = scratch_dir
)

#Bring stream polygon into r environment
streams<-st_read(paste0(scratch_dir, "streams.shp"))
st_crs(streams)<-st_crs(dem@crs)

#Crop to watershed
streams<-streams[shed,]

#2.2.2 Estimate stream characteristics -----------------------------------------
#distance to outlet
wbt_distance_to_outlet(
  d8_pntr = 'fdr.tif',
  streams = "streams.tif",
  output = 'dist2out.tif',
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

#2.2.3 Convert stream layer to points-------------------------------------------
#Bring rasters into R env
twi<-raster(paste0(scratch_dir,"twi.tif"))
  twi<-crop(twi, shed)
  twi<-mask(twi, shed)
fac<-raster(paste0(scratch_dir,"fac.tif"))
  fac<-crop(fac, shed)
  fac<-mask(fac, shed)
  fac<-fac*res(fac)[1]*res(fac)[2]
dist2out<-raster(paste0(scratch_dir, "dist2out.tif"))
  dist2out<-crop(dist2out, shed)
  dist2out<-mask(dist2out, shed)
  dist2out<-dist2out-raster::cellStats(dist2out, min)

#Define reach
streams<-streams %>% 
  arrange(FID) %>% 
  mutate(StreamReach = seq(1, nrow(.))) 
  
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
    drain_area_ha = raster::extract(fac, .)/10000,
    pid = seq(1, nrow(.)))    # a unique identifier for each point

#Crop pnts to shed
pnts<-pnts[shed,]

#Remove points within 25 m of juncture 
pnts_junctures<-pnts %>% 
  group_by(StreamReach) %>% 
  slice_max(drain_area_ha, n=1, with_ties = F) %>%
  st_buffer(., dist=25) %>% 
  st_intersection(pnts,.) %>% 
  st_drop_geometry %>% 
  select(pid) %>% 
  pull()
pnts<-pnts %>% filter(!(pid %in% pnts_junctures))

#2.2.4 Define network connectivity----------------------------------------------
streams$drainage_area<-raster::extract(fac,streams, fun=max, na.rm=T)

#Create function to define conn by reach
con_fun<-function(n){
  
  #Identify reach of interest
  reach<-streams[n,]
  
  #Define reach id
  StreamReach<-reach$StreamReach
  
  #Identify intersecting reaches
  streams_temp<-streams[reach,]
  
  #Identify downstream reach
  Flows_To<-streams_temp %>% 
    #Kill geo
    st_drop_geometry() %>% 
    #Remove study reach
    filter(drainage_area>reach$drainage_area) %>% 
    #Identify reach w/ largest ws area
    slice_max(drainage_area, n=1, with_ties = F) %>% 
    #select uid 
    select(StreamReach) %>% pull()
  
  #export point
  tibble(StreamReach, Flows_To)
}

#Apply function and join to streams fun
con<-lapply(seq(1,nrow(streams)), con_fun) %>% bind_rows() %>% distinct()
streams <- streams %>% left_join(., con)

#Remove drainage area estimate based (sometimes this is innacurate due to raster resolution)
streams<-streams %>% select(-drainage_area)

# #2.2.5 Estimate drainage area of each reach ------------------------------------
#Identify drainage area based on points
drainage_area<-pnts %>%
  #Drop geomoetry info
  st_drop_geometry() %>%
  #group by stream reach
  group_by(StreamReach) %>%
  #Select point with min downstream dist
  slice_max(drain_area_ha, with_ties = F) %>%
  #Select
  select(StreamReach, drain_area_ha)

#Join to streams layer
streams<-left_join(streams, drainage_area)
  
#2.3 Identify sampling points---------------------------------------------------
#2.3.1 Define main stem (longest connected stream length) ----------------------
#Identify headwater reach of main stem
main_stem<-pnts %>% 
  st_drop_geometry() %>% 
  slice_max(dist2out, n=1, with_ties = F) %>% 
  select(StreamReach) %>% 
  slice(n=1) %>% 
  pull()

#Use while loop to identify downstream reaches
n<-1
while(length(main_stem)==n){
  temp<-streams %>% 
    st_drop_geometry() %>% 
    filter(StreamReach %in% main_stem) %>% 
    select(Flows_To) %>% pull()
  main_stem<-c(main_stem, temp) %>% unique()
  n<-n+1
}
main_stem<-streams %>% filter(StreamReach %in% main_stem)

#2.3.2 Identify tributary stations ---------------------------------------------
#What reaches intersect the main stem?
tribs<-streams[main_stem,]

#Crop out main stem
tribs<-tribs %>% 
  filter(!(StreamReach %in% main_stem$StreamReach))

#Select top two tribs
tribs<-tribs %>% 
  slice_max(drain_area_ha, n=2, with_ties = F)

#Define sampling point pid
trib_pnts<-pnts %>% 
  st_drop_geometry() %>% 
  filter(StreamReach %in% tribs$StreamReach) %>% 
  group_by(StreamReach) %>% 
  slice_max(drain_area_ha, n=1, with_ties = F) %>% 
  select(pid) %>% pull()

#2.3.3 Define stations upstream of tribs ---------------------------------------
# #define upstream tribs
# upstream_trib<-main_stem[tribs,] 
# upstream_trib<-upstream_trib %>% 
#   filter(!(StreamReach %in% tribs$StreamReach)) %>% 
#   slice_min(drain_area_ha, n=1, with_ties = F)
# #Identify points
# upstream_trib_pnts<-pnts %>% 
#   st_drop_geometry() %>% 
#   filter(StreamReach %in% upstream_trib$StreamReach) %>% 
#   slice_max(drain_area_ha, n=1, with_ties = F) %>% 
#   select(pid) %>% 
#   pull()

#2.3.4 Define stations along main stem -----------------------------------------
#Define main stem points
main_stem_pnts <- pnts[shed,] %>% 
  st_drop_geometry() %>% 
  filter(StreamReach %in% main_stem$StreamReach) %>% 
  select(pid, drain_area_ha) 

#Subset based on % of watershed area
ws_area<-(st_area(shed)/10000) %>% paste %>% as.numeric()
ws_area<-c(ws_area*main_stem_stn[1], 
           ws_area*main_stem_stn[2], 
           ws_area*main_stem_stn[3], 
           ws_area*main_stem_stn[4])
main_stem_pnts<-main_stem_pnts %>% 
  mutate(
    q25 = abs(drain_area_ha-ws_area[1]),
    q50 = abs(drain_area_ha-ws_area[2]),
    q75 = abs(drain_area_ha-ws_area[3]),
    q90 = abs(drain_area_ha-ws_area[4])) %>% 
  select(-drain_area_ha) %>% 
  pivot_longer(-pid, names_to = 'quant', values_to = 'diff') %>% 
  group_by(quant) %>% 
  slice_min(diff,n=1,with_ties=F) %>% 
  ungroup() %>% 
  select(pid) %>% 
  pull()

#2.3.5 Export points -----------------------------------------------------------
#Combine points
output<-pnts %>% filter(pid %in% c(trib_pnts, main_stem_pnts))
output<-st_as_sf(output)

#Export
list(shed, streams, output)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 3: Watershed Analysis ----------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Shambly --------------------------------------------------------------
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

#Run function
output<-sampling_station_fun(dem, pp, threshold = 50000, scratch_dir)

#Map
shed<-output[[1]]
streams<-output[[2]]
stations<-output[[3]]
m<-mapview(shed) + mapview(streams) + mapview(stations) + mapview(pp)
m
mapshot(m, "docs/lts_shambly.html")

#Export points 
st_write(stations, paste0(output_dir, "lts_shambly.shp"))

#3.2 Paint Rock ----------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_paintrock//"

#Define data inputs
dem<-raster(paste0(data_dir,"dem_mask_10m.tif"))
pp<-st_read(paste0(data_dir, "watershed_outlets.shp")) %>% 
  filter(Name == "Fanning Hollow") %>% 
  st_transform(., crs = st_crs(dem@crs))

#Run function
output<-sampling_station_fun(dem, pp, threshold = 1000, scratch_dir)

#Map
shed<-output[[1]]
streams<-output[[2]]
stations<-output[[3]]
m<-mapview(shed) + mapview(streams) + mapview(stations) + mapview(pp)
m
mapshot(m, "docs/lts_paintrock.html")

#Export points 
st_write(stations, paste0(output_dir, "lts_paintrock.shp"))

#3.3 Talladega------------------------------------------------------------------
#Define Data Directory
data_dir<-"C://WorkspaceR//AIMS_watersheds//data//I_data_Talladega//"

#Define data inputs
dem<-raster(paste0(data_dir,"elevation//lidar01m33085g5.tif"))
  crop<-st_read(paste0(data_dir, "crop2.shp"))
  dem<-crop(dem,crop)
pp<-st_read(paste0(data_dir,"pp.shp"))

#Run function
output<-sampling_station_fun(dem, pp, threshold = 35000, scratch_dir)

#Map
shed<-output[[1]]
streams<-output[[2]]
stations<-output[[3]]
m<-mapview(shed) + mapview(streams) + mapview(stations) + mapview(pp)
m
mapshot(m, "docs/lts_talladega.html")

#Export points 
st_write(stations, paste0(output_dir, "lts_talladega.shp"))

