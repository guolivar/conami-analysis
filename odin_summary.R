#'---
#'title: "ODIN-SD Idaho 2019 Summary"
#'author: Gustavo Olivares
#'---
#'This script parses and cleans the ODIN-SD data from their memory cards to then calculate summary statistics 
#'and plot draft maps of PM~10~ and PM~2.5~
#'

#'
#' ## Prepare libraries
#'  
library(librarian) # To more flexibly manage packages
shelf(readr,
      openair,
      automap,
      raster,
      gstat,
      sp,
      rgdal,
      ggmap,
      ggplot2,
      scales)
#' ## Set constants
data_path <- path.expand("~/data/ODIN_IDAHO/2019/RAW/FromPortal/")
data_path.out <- path.expand("~/data/ODIN_IDAHO/2019/WORKING/")
files_list <- dir(data_path,pattern = 'SD')
# Define time average for output
tavg <- '1 hour'
#' ## Load data
#'
#'  Get devices locations
#'  
odin_locations <- readr::read_delim(paste0(data_path,"odin_locations.txt"), 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
#'
#'  Cycle through the folders to work with all the DATA.TXT files
#'  

for (i in (1:length(files_list))){
  file <- files_list[i]
  print(file)
  odin.data <- readr::read_delim(paste0(data_path,file),
                                 delim = ',',
                                 skip = 1,
                                 col_names = c('date',
                                               'PM1',
                                               'PM2.5',
                                               'PM10',
                                               'Temperature',
                                               'RH'))
  # Add serial number to the dataframe
  odin.data$ODINsn <- substr(file,1,6)
  # Add location to the dataframe
  # Find the location of the current ODIN
  odin_id <- which(odin_locations$ODIN == substr(file,1,6))
  odin.data$lat <- odin_locations$lat[odin_id]
  odin.data$lon <- -odin_locations$lon[odin_id]

  # Construct the ALLDATA frame
  
  if (i == 1){
    # This is the first iteration so we just copy the "odin.data" dataframe
    all.data <- odin.data
    all.data.tavg <- timeAverage(odin.data,avg.time = tavg)
    all.data.tavg$ODINsn <- odin.data$ODINsn[1]
  } else {
    # We already have "all.data" so we need to append the current "odin.data"
    all.data <- rbind(all.data,odin.data)
    tmp1 <- timeAverage(odin.data,avg.time = tavg)
    tmp1$ODINsn <- odin.data$ODINsn[1]
    all.data.tavg <- rbind(all.data.tavg,tmp1)
    # Remove all.data to clean for next iteration
    rm(odin.data)
  }
}

#'
#' ## Summary statistics
#'
#'Note that the statistics are calculated on 10 minutes averages and that the campaign
#'went from January 29th 16:00 UTC until March 1st 15:50 UTC
#'
#'The units for the parameters are:
#' * PM~2.5~ [$mu$g/m^3^]
#' * PM~10~ [$mu$g/m^3^]
#' * Temperature [celsius]
#' * RH [%]
#' 

# Calculate the summary table for each unit
summary_mean <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data, FUN = mean)
summary_max <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data, FUN = max)
summary_min <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data, FUN = min)
summary_sd <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data, FUN = sd)
summary_N <- aggregate(cbind(Temperature, RH,PM2.5, PM10) ~ODINsn, all.data, FUN = length)

#' 
#' ## Average concentrations
#' 
print(format(summary_mean,digits = 1))
#' 
#' ## Maximum concentrations
#' 
print(format(summary_max,digits = 1))
#' 
#' ## Minimum concentrations
#' 
print(format(summary_min,digits = 1))
#' 
#' ## Standard deviation
#' 
print(format(summary_sd,digits = 1))

#'
#' ## Time series
#' 

#'
#' ### PM~2.5~
#' 
ggplot(data = all.data, aes(x=date)) +
  geom_line(aes(y=PM2.5,colour = ODINsn))
#'
#' ### PM~10~
#' 
ggplot(data = all.data, aes(x=date)) +
  geom_line(aes(y=PM10,colour = ODINsn))
#'
#' ### Temperature
#' 
ggplot(data = all.data, aes(x=date)) +
  geom_line(aes(y=Temperature,colour = ODINsn))
#'
#' ### Relative Humidity
#' 
ggplot(data = all.data, aes(x=date)) +
  geom_line(aes(y=RH,colour = ODINsn))

ggplot(data = all.data, aes(x=RH)) +
  geom_point(aes(y=PM2.5,colour = ODINsn))

ggplot(data = all.data, aes(x=Temperature)) +
  geom_point(aes(y=PM2.5,colour = ODINsn))


#'
#' ## Average Maps
#' 

# Some useful constants
proj4string_NZTM <- CRS('+init=epsg:2193')
proj4string_latlon <- CRS('+init=epsg:4326')
# Calculate the summaries for cold and warm periods
temperature_threshold <- 5
summary_mean_map_cold <- aggregate(cbind(Temperature, RH, PM2.5, PM10, lon, lat) ~ODINsn,
                                   subset(all.data,Temperature < temperature_threshold),
                                   FUN = mean)
summary_mean_map_warm <- aggregate(cbind(Temperature, RH, PM2.5, PM10, lon, lat) ~ODINsn,
                                   subset(all.data,Temperature >= temperature_threshold),
                                   FUN = mean)

#' COLD ( < 5 C)
# We need to remove the ODINs not in Kellogg
summary_mean_map_cold <- subset(summary_mean_map_cold,ODINsn != "SD0006")
summary_mean_map_cold <- subset(summary_mean_map_cold,ODINsn != "SD0067")
summary_mean_map_cold <- subset(summary_mean_map_cold,ODINsn != "SD0060")
# Assign coordinates to the dataframe
coordinates(summary_mean_map_cold) <- ~ lon + lat
proj4string(summary_mean_map_cold) <- proj4string_latlon

# Get the basemap
centre_lat <- mean(summary_mean_map_cold$lat)
centre_lon <- mean(summary_mean_map_cold$lon)
ca <- get_googlemap(
  c(lon=centre_lon,lat=centre_lat),
  zoom=15,
  scale=2,
  color="bw",
  key = "AIzaSyACi3pNvPQTxZWx5u0nTtke598dPqdgySg")

#' ### PM~2.5~
ggmap(ca) + 
  geom_point(data=as.data.frame(summary_mean_map_cold),aes(x=lon,y=lat,colour = PM2.5),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map_cold$PM2.5)),
                          name = "PM2.5", oob=squish)

#' ### PM~10~
ggmap(ca) + 
  geom_point(data=as.data.frame(summary_mean_map_cold),aes(x=lon,y=lat,colour = PM10),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map_cold$PM10)),
                          name = "PM10", oob=squish)

#' ### PM~coarse~
ggmap(ca) + 
  geom_point(data=as.data.frame(summary_mean_map_cold),aes(x=lon,y=lat,colour = PM10 - PM2.5),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map_cold$PM10 - summary_mean_map_cold$PM2.5)),
                          name = "PMcoarse", oob=squish)

#' WARM ( <= 5 C)
# We need to remove the ODINs not in Kellogg
summary_mean_map_warm <- subset(summary_mean_map_warm,ODINsn != "SD0006")
summary_mean_map_warm <- subset(summary_mean_map_warm,ODINsn != "SD0067")
summary_mean_map_warm <- subset(summary_mean_map_warm,ODINsn != "SD0060")
# Assign coordinates to the dataframe
coordinates(summary_mean_map_warm) <- ~ lon + lat
proj4string(summary_mean_map_warm) <- proj4string_latlon

# Get the basemap
centre_lat <- mean(summary_mean_map_warm$lat)
centre_lon <- mean(summary_mean_map_warm$lon)
ca <- get_googlemap(
  c(lon=centre_lon,lat=centre_lat),
  zoom=15,
  scale=2,
  color="bw",
  key = "AIzaSyACi3pNvPQTxZWx5u0nTtke598dPqdgySg")

#' ### PM~2.5~
ggmap(ca) + 
  geom_point(data=as.data.frame(summary_mean_map_warm),aes(x=lon,y=lat,colour = PM2.5),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map_warm$PM2.5)),
                          name = "PM2.5", oob=squish)

#' ### PM~10~
ggmap(ca) + 
  geom_point(data=as.data.frame(summary_mean_map_warm),aes(x=lon,y=lat,colour = PM10),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map_warm$PM10)),
                          name = "PM10", oob=squish)

#' ### PM~coarse~
ggmap(ca) + 
  geom_point(data=as.data.frame(summary_mean_map_warm),aes(x=lon,y=lat,colour = PM10 - PM2.5),size = 5) +
  scale_colour_continuous(low="white", high="red",limits=c(0, max(summary_mean_map_warm$PM10 - summary_mean_map_warm$PM2.5)),
                          name = "PMcoarse", oob=squish)



# Save the "all.data" dataframe
save(all.data,file = paste0(data_path.out,'alldata.RData'))
save(all.data.tavg,file = paste0(data_path.out,'alldataTAVG.RData'))
data.output <- all.data.tavg[,c('date','PM2.5','PM10','Temperature','RH','ODINsn')]

