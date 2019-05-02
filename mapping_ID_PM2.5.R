#########################################################################
# Taking cleaned, averaged ODIN data and generate maps and animations.
# Developed for the project "Air Quality Measurements in Gisborne - 2018"
# Authors:
#         Gustavo Olivares
#########################################################################
#########################################################################
##### Load relevant packages #####
library(librarian) # To more flexibly manage packages
shelf(readr,
      reshape2,
      automap,
      raster,
      gstat,
      sp,
      rgdal,
      ggmap,
      ggplot2,
      scales,
      gstat,
      RNetCDF,
      RJSONIO,
      curl,
      base64enc,
      zoo,
      openair,
      stringi,
      viridis,
      dplyr,
      RColorBrewer,
      purrr,
      magick)
# Set the data folder ####
data_path <- path.expand("~/data/ODIN_IDAHO/2019/WORKING/")
plot_path <- path.expand("~/data/ODIN_IDAHO/2019/WORKING/PLOTS/")

# Load data
load(paste0(data_path,'alldata.RData'))
load(paste0(data_path,'alldataTAVG.RData'))
# Focus only on Kellogg
all.data <- subset(all.data,ODINsn != "SD0006")
all.data <- subset(all.data,ODINsn != "SD0067")
all.data <- subset(all.data,ODINsn != "SD0060")
all.data.tavg <- subset(all.data.tavg,ODINsn != "SD0006")
all.data.tavg <- subset(all.data.tavg,ODINsn != "SD0067")
all.data.tavg <- subset(all.data.tavg,ODINsn != "SD0060")
# Some useful constants
proj4string_latlon <- CRS('+init=epsg:4326')
# Idaho Transverse Mercator - NAD83
proj4string_IDTM <- CRS("+proj=tmerc +lat_0=42 +lon_0=-114 +k=0.9996 +x_0=2500000 +y_0=1200000 +ellps=GRS80 +units=m +no_defs")
# Assign coordinates to the dataframes
coordinates(all.data) <- ~ lon + lat
proj4string(all.data) <- proj4string_latlon
coordinates(all.data.tavg) <- ~ lon + lat
proj4string(all.data.tavg) <- proj4string_latlon

centre_lat <- mean(all.data.tavg$lat)
centre_lon <- mean(all.data.tavg$lon)

## Prepare the map to plot animations #####

# Get the basemap
ca <- get_googlemap(
  c(lon=centre_lon,lat=centre_lat),
  zoom=15,
  scale=2,
  color="bw",
  key = "AIzaSyACi3pNvPQTxZWx5u0nTtke598dPqdgySg",
  maptype = "terrain")

# Re-project to IDTM #####
all.data.tavg <- spTransform(all.data.tavg,proj4string_IDTM)
print("Starting the interpolation")

#Setting the  prediction grid properties #####
cellsize <- 100 #pixel size in projection units (IDTM, i.e. metres)
min_x <- all.data.tavg@bbox[1,1] - cellsize#minimun x coordinate
min_y <- all.data.tavg@bbox[2,1] - cellsize #minimun y coordinate
max_x <- all.data.tavg@bbox[1,2] + cellsize #mximum x coordinate
max_y <- all.data.tavg@bbox[2,2] + cellsize #maximum y coordinate

x_length <- max_x - min_x #easting amplitude
y_length <- max_y - min_y #northing amplitude

ncol <- round(x_length/cellsize,0) #number of columns in grid
nrow <- round(y_length/cellsize,0) #number of rows in grid

grid <- GridTopology(cellcentre.offset=c(min_x,min_y),cellsize=c(cellsize,cellsize),cells.dim=c(ncol,nrow))

#Convert GridTopolgy object to SpatialPixelsDataFrame object. #####
grid <- SpatialPixelsDataFrame(grid,
                               data=data.frame(id=1:prod(ncol,nrow)),
                               proj4string=proj4string_IDTM)


# Get rid of NA containing rows
all.data.tavg <- subset(all.data.tavg,!is.na(PM2.5))
all_dates <- sort(unique(all.data.tavg$date))
valid_dates <- FALSE * (1:length(all_dates))
# limits for colorscales #####
cmin <- min(all.data.tavg$PM2.5)
cmax <- max(all.data.tavg$PM2.5) * 0.5
## Interpolate and plot #####
ndates <- length(all_dates)
breaks <- as.numeric(quantile((1:ndates),c(0,0.5,1), type = 1))
nbreaks <- length(breaks)
i <- 0
for (d_slice in (1:ndates)){
  c_data <- subset(all.data.tavg,subset = (date==all_dates[d_slice]))
  
  if (length(unique(c_data$ODINsn))<2){
    next
  }
  valid_dates[d_slice] <- TRUE
  
  surf.idw <- idw(PM2.5 ~ 1,newdata = grid, locations = c_data, idp = 1,na.action = na.omit)
  surf.idw$timestamp <-d_slice
  proj4string(surf.idw) <- proj4string_IDTM
  
  
  if (i==0){
    
    to_rast.idw <- surf.idw
    r0.idw <- rasterFromXYZ(cbind(surf.idw@coords,surf.idw$var1.pred))
    crs(r0.idw) <- proj4string_IDTM
    raster_cat.idw<- r0.idw
    
    i <- 1
  }
  else {
    
    to_rast.idw <- surf.idw
    r0.idw <- rasterFromXYZ(cbind(surf.idw@coords,surf.idw$var1.pred))
    names(r0.idw) <- as.character(all_dates[d_slice])
    crs(r0.idw) <- proj4string_IDTM
    raster_cat.idw<- addLayer(raster_cat.idw,r0.idw)
    
  }
  rtp <- rasterToPolygons(projectRaster(r0.idw,crs = proj4string_latlon))
  
  points <- data.frame(spTransform(c_data,proj4string_latlon))
  
  # Build the animation
  map_out <- ggmap(ca) + geom_polygon(data = rtp,aes(x = long, y = lat, group = group,
                                                     fill = rep(rtp[[1]], each = 5)),
                                      size = 0,
                                      alpha = 0.85) +
    scale_fill_gradient(low="white",
                        high="red",
                        limits=c(0, cmax),
                        name = "PM2.5",
                        oob=squish) +
    geom_point(data=points,
               aes(x=lon,y=lat),
               colour = "black") +
    ggtitle(paste(as.character(all_dates[d_slice]),
                  "UTC"))
  
  ggsave(filename=paste0(plot_path,
                         'idw/',
                         format(all_dates[d_slice],
                                format = "%Y-%m-%d %H:%M"),
                         '.png'),
         plot=map_out,
         width=6,
         height=6,
         units = "in")
}
save('raster_cat.idw',file = paste0(data_path,'raster_cat.idw.RData'))

print("Done with interpolating ...")

raster_cat_idw_LL <- projectRaster(raster_cat.idw,crs = proj4string_latlon)

save(list = c('raster_cat_idw_LL'),file = paste0("raster_odin_LL_IDW.RData"))

# Write NetCDF files ####
# IDW
lat_dim <- unique(coordinates(raster_cat_idw_LL)[,2])
lon_dim <- unique(coordinates(raster_cat_idw_LL)[,1])
tim_dim <- all_dates[valid_dates==1]
nc.idw <- create.nc(paste0(path.expand(plot_path),"odin_idw.nc"))
# Dimensions specifications
dim.def.nc(nc.idw, "time", unlim=TRUE)
dim.def.nc(nc.idw, "latitude",length(lat_dim))
dim.def.nc(nc.idw, "longitude",length(lon_dim))
# Variable specifications
var.def.nc(nc.idw,"time","NC_INT","time")
att.put.nc(nc.idw,"time","units","NC_CHAR","seconds since 1970-01-01 00:00:0.0")
att.put.nc(nc.idw,"time","long_name","NC_CHAR","time")

var.def.nc(nc.idw,"latitude","NC_FLOAT","latitude")
att.put.nc(nc.idw,"latitude","units","NC_CHAR","degrees_north")
att.put.nc(nc.idw,"latitude","long_name","NC_CHAR","latitude")
att.put.nc(nc.idw,"latitude","standard_name","NC_CHAR","latitude")

var.def.nc(nc.idw,"longitude","NC_FLOAT","longitude")
att.put.nc(nc.idw,"longitude","units","NC_CHAR","degrees_east")
att.put.nc(nc.idw,"longitude","long_name","NC_CHAR","longitude")
att.put.nc(nc.idw,"longitude","standard_name","NC_CHAR","longitude")

var.def.nc(nc.idw,"pm2p5","NC_FLOAT",c("longitude","latitude","time"))
att.put.nc(nc.idw,"pm2p5","units","NC_CHAR","ug m**-3")
att.put.nc(nc.idw,"pm2p5","long_name","NC_CHAR","Mass concentration of PM2.5 ambient aerosol particles in air")
att.put.nc(nc.idw,"pm2p5","standard_name","NC_CHAR","mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air")
att.put.nc(nc.idw,"pm2p5","cell_methods","NC_CHAR","time: mean (interval: 15 minutes)")
att.put.nc(nc.idw,"pm2p5","missing_value","NC_FLOAT",-999.9)

# Global attributes
att.put.nc(nc.idw,"NC_GLOBAL","title","NC_CHAR","PM2.5 interpolated surface (Inverse Square Distance)")
att.put.nc(nc.idw,"NC_GLOBAL","Conventions","NC_CHAR","CF-1.7")
att.put.nc(nc.idw,"NC_GLOBAL","Institution","NC_CHAR","NIWA (National Institute of Water and Atmospheric Research, Auckland, New Zealand)")
att.put.nc(nc.idw,"NC_GLOBAL","project_id","NC_CHAR","CONA - 2018")
att.put.nc(nc.idw,"NC_GLOBAL","history","NC_CHAR",paste0(format(max(all.data.tavg$date),format = "%Y%m%d"),
                                                         " Data generated and formatted"))
att.put.nc(nc.idw,"NC_GLOBAL","comment","NC_CHAR","Data for visualisation only")

# Load data
var.put.nc(nc.idw,"latitude",lat_dim)
var.put.nc(nc.idw,"longitude",lon_dim)
var.put.nc(nc.idw,"time",as.numeric(tim_dim))
rast_data <- getValues(raster_cat_idw_LL)[,(1:length(tim_dim))]
dim(rast_data) <- c(length(lon_dim),
                    length(lat_dim),
                    length(tim_dim))
var.put.nc(nc.idw,"pm2p5",rast_data)

# Close the file and save
close.nc(nc.idw)

## Create MP4 video ####
system(paste0("ffmpeg -f image2 -r 6 -pattern_type glob -i '",
              path.expand(plot_path),
              "idw/*.png' ",
              path.expand(plot_path),
              "idw/",
              format(min(all.data.tavg$date),
                     format = "%Y%m%d"),"_",
              format(max(all.data.tavg$date),
                     format = "%Y%m%d"),
              ".mp4"))

## Upload data ####

RCurl::ftpUpload(paste0(path.expand(plot_path),"odin_idw.nc"),
                 "ftp://ftp.niwa.co.nz/incoming/GustavoOlivares/Idaho/odin_idw_idaho_PM2.5.nc")

RCurl::ftpUpload(paste0(path.expand(plot_path),
                        "idw/",
                        format(min(all.data.tavg$date),
                               format = "%Y%m%d"),"_",
                        format(max(all.data.tavg$date),
                               format = "%Y%m%d"),
                        ".mp4"),
                 "ftp://ftp.niwa.co.nz/incoming/GustavoOlivares/Idaho/odin_idw_idaho_PM2.5.mp4")

