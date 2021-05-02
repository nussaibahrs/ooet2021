# This script divides the world in different ocean basins to allow for 
# per-basis analyses.

# Load libraries ----------------------------------------------------------

library(raster)
library(chronosphere)

# Create and save ocean basins file ------------------------------------------------

world <- chronosphere::fetch("NaturalEarth")

atlantic1 <- extent(c(                
  xmin = -92,
  xmax = 10.64,
  ymin = 11.78,
  ymax = 79.38 
)) 

atlantic2<- extent(c(                
  xmin = -65.52,
  xmax = 26.35,
  ymin = -90,
  ymax = 15.41 
)) 

a1 <- as(atlantic1, 'SpatialPolygons')  
a2 <- as(atlantic2, "SpatialPolygons")

atlantic <- raster::union(a1,a2)

pacific1 <- extent(c(                
  xmin = -155,
  xmax = -65.52,
  ymin = -84.11,
  ymax = 10.89 
)) 

pacific2 <- extent(c(                
  xmin = -155,
  xmax = -96.88,
  ymin = 10,
  ymax = 71.14 
)) 

p1 <- as(pacific1, 'SpatialPolygons')  
p2 <- as(pacific2, "SpatialPolygons")
westpacific <- raster::union(p1,p2)

pacific3<- extent(c(                
  xmin = 107.32,
  xmax = 178.71,
  ymin =0,
  ymax = 66.5 
)) 
pacific4<- extent(c(                
  xmin = 142.63,
  xmax = 178.71,
  ymin = -86.43,
  ymax = 10 
)) 

pacific5 <- extent(c(                
  xmin = -180.55,
  xmax = -155,
  ymin = -86.43,
  ymax = 73.46 
))

p3 <- as(pacific3, 'SpatialPolygons')  
p4 <- as(pacific4, "SpatialPolygons")
p5 <- as(pacific5, "SpatialPolygons")

eastpacific <- raster::union(p5, raster::union(p3,p4))


indian1 <- extent(c(                
  xmin = 24.41,
  xmax = 141.86,
  ymin = -67.89,
  ymax = -30.82 
)) 

indian2<- extent(c(                
  xmin = 25.18,
  xmax = 107.32,
  ymin = -35.45,
  ymax = 24.79 
)) 

indian3 <- extent(c(                
  xmin = 102.71,
  xmax = 133.42,
  ymin = -33.14,
  ymax = 0.46 
)) 


i1 <- as(indian1, 'SpatialPolygons')  
i2 <- as(indian2, "SpatialPolygons")
i3 <- as(indian3, "SpatialPolygons")

indian <- raster::union(raster::union(i1,i2), i3)

basins <- list(atlantic=atlantic, 
               westpacific=westpacific, 
               eastpacific=eastpacific, 
               indian=indian)
save(basins, file=file.path("output", "ocean_basins.RData"))
    