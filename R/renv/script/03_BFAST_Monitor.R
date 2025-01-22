library(terra)
library(bfast)
library(zoo)
library(lubridate)
library(dplyr)

getwd()

modis_dir <- "/Users/osako/github/aeo-project/python/data/modis_images"

load_and_sort_tifs <- function(directory) {
  tif_files <- list.files(path = directory, pattern = "*.tif$", full.names = TRUE)
  file_dates <- ymd(basename(tif_files))
  sorted_indices <- order(file_dates)   
  list(files = tif_files[sorted_indices], dates = file_dates[sorted_indices])
}

modis_data <- load_and_sort_tifs(modis_dir)
modis_rasters <- rast(modis_data$files)

sorted_indices <- order(modis_data$dates)
file_dates_sorted <- modis_data$dates[sorted_indices]
rasters_sorted <- modis_rasters[[sorted_indices]]
dates <- file_dates_sorted

timeser <- function(val_array, time_array) {
  z <- zoo(val_array, time_array)
  yr <- as.numeric(format(time(z), "%Y"))
  jul <- as.numeric(format(time(z), "%j"))
  delta <- min(unlist(tapply(jul, yr, diff)))
  zz <- aggregate(z, yr + (jul - 1) / delta / 23)
  tso <- as.ts(zz)
  return(tso)
}

bfmRaster = function(pixels) {
  tspx <- timeser(pixels, dates)
  bfm <- bfastmonitor(tspx, 
                      response ~ trend + harmon, 
                      order = 1, 
                      start = c(2020,1))
  return(c(bfm$breakpoint, bfm$magnitude))
}


bfm_results <- app(rasters_sorted, bfmRaster)
names(bfm_results) <- c('time of break', 'magnitude of change')

plot(bfm_results)

writeRaster(bfm_results, "/Users/osako/github/aeo-project/R/renv/data/bfm_first.tif", overwrite=TRUE)
