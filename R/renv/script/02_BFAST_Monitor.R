library(terra)
library(bfast)
library(zoo)
library(lubridate)
library(dplyr)

getwd()

modis_dir <- "/Users/osako/github/aeo-project/python/data/modis_images"
landsat_dir <- "/Users/osako/github/aeo-project/python/data/landsat_images"

load_and_sort_tifs <- function(directory) {
  tif_files <- list.files(path = directory, pattern = "*.tif$", full.names = TRUE)
  file_dates <- ymd(basename(tif_files))
  sorted_indices <- order(file_dates)   
  list(files = tif_files[sorted_indices], dates = file_dates[sorted_indices])
}

modis_data <- load_and_sort_tifs(modis_dir)
landsat_data <- load_and_sort_tifs(landsat_dir)

modis_rasters <- rast(modis_data$files)
landsat_rasters <- rast(landsat_data$files)

combined_rasters <- c(landsat_rasters, modis_rasters)
combined_dates <- c(landsat_data$dates, modis_data$dates)
sorted_indices <- order(combined_dates)
file_dates_sorted <- combined_dates[sorted_indices]
combined_rasters_sorted <- combined_rasters[[sorted_indices]]
dates <- file_dates_sorted

timeser <- function(val_array, time_array) {
  z <- zoo(val_array, time_array)
  yr <- as.numeric(format(time(z), "%Y"))
  jul <- as.numeric(format(time(z), "%j"))
  delta <- min(unlist(tapply(jul, yr, diff)), na.rm = TRUE)
  zz <- aggregate(z, yr + (jul - 1) / delta / 16)
  tso <- as.ts(zz)
  return(tso)
}

bfmRaster <- function(pixels) {
  df <- data.frame(date = dates, EVI = pixels)
  df <- df %>% filter(EVI >= -10 & EVI < 10)
  
  
  df_monthly <- df %>%
    mutate(year_month = floor_date(date, "month")) %>%
    group_by(year_month) %>%
    summarize(EVI_mean = mean(EVI, na.rm = TRUE)) %>%
    ungroup()
  
  ts_monthly <- try(ts(df_monthly$EVI_mean, 
                       start = c(year(min(df_monthly$year_month)), 
                                 month(min(df_monthly$year_month))), 
                       frequency = 12))
  
  if (inherits(ts_monthly, "try-error")) {
    return(c(NA, NA))
  }
  
  bfm <- try(bfastmonitor(ts_monthly,
                          response ~ trend + harmon,
                          order = 2,
                          start = c(2001, 1)))
  
  if (inherits(bfm, "try-error")) {
    return(c(NA, NA))
  }
  
  return(c(bfm$breakpoint, bfm$magnitude))
}


e <- ext(combined_rasters_sorted)
subset_raster <- crop(combined_rasters_sorted,  ext(102.0,102.08,19.99,20.03))
bfmR_subset <- app(subset_raster, bfmRaster)
plot(bfmR_subset)

bfmR <- app(combined_rasters_sorted, bfmRaster)
names(bfmR) <- c('time of break', 'magnitude of change')
plot(bfmR)

writeRaster(bfmR, "/Users/osako/Downloads/bfm_result.tif", overwrite=TRUE)
