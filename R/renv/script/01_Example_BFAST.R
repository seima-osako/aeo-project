# pkgTest is a helper function to load packages and install packages only when they are not installed yet.
pkgTest <- function(x)
{
  if (x %in% rownames(installed.packages()) == FALSE) {
    install.packages(x, dependencies= TRUE)
  }
  library(x, character.only = TRUE)
}

neededPackages <- c("zoo", "bfast", "terra", "raster", "leaflet", "MODISTools")
for (package in neededPackages){pkgTest(package)}
# Utility function to create time series object from a numeric vector
# val_array: data array for one single pixel (length is number of time steps)
# time_array: array with Dates at which raster data is recorded (same length as val_array)
timeser <- function(val_array, time_array) {
  z <- zoo(val_array, time_array) # create zoo object
  yr <- as.numeric(format(time(z), "%Y")) # extract the year numbers
  jul <- as.numeric(format(time(z), "%j")) # extract the day numbers (1-365)
  delta <- min(unlist(tapply(jul, yr, diff))) # calculate minimum time difference (days) between observations
  zz <- aggregate(z, yr + (jul - 1) / delta / 23) # aggregate into decimal year timestamps
  (tso <- as.ts(zz)) # convert into timeseries object
  return(tso)
}

# Downloading the NDVI data, starting from 2000-01-01
VI <- mt_subset(product = "MOD13Q1",
                site_id = "nl_gelderland_loobos",
                band = "250m_16_days_NDVI",
                start = "2000-01-01",
                end = "2022-03-22",
                km_lr = 2,
                km_ab = 2,
                site_name = "testsite",
                internal = TRUE,
                progress = TRUE)


# Downloading the pixel reliability data, starting from 2000-01-01
QA <- mt_subset(product = "MOD13Q1",
                site_id = "nl_gelderland_loobos",
                band = "250m_16_days_pixel_reliability",
                start = "2000-01-01",
                end = "2022-03-22",
                km_lr = 2,
                km_ab = 2,
                site_name = "testsite",
                internal = TRUE,
                progress = TRUE)

VI_r <- mt_to_terra(df = VI)
QA_r <- mt_to_terra(df = QA)

## clean the data
# create mask on pixel reliability flag set all values <0 or >1 NA
m <- QA_r
m[(QA_r < 0 | QA_r > 1)] <- NA # continue working with QA 0 (good data), and 1 (marginal data)

# apply the mask to the NDVI raster
VI_m <- mask(VI_r, m, maskvalue=NA, updatevalue=NA)

# plot the first image
plot(m,1) # plot mask

plot(VI_m,1) # plot cleaned NDVI raster

px <- 78 # pixel number; adjust this number to select the center pixel
tspx <- timeser(unlist(VI_m[px]),as.Date(names(VI_m), "%Y-%m-%d")) # convert pixel "px" to a time series
plot(tspx, main = 'NDVI') # NDVI time series cleaned using the "reliability information"

bfm1 <- bfastmonitor(tspx, response ~ trend + harmon, order = 3, start = c(2018,3)) # Note: the third observation in 2018 marks the transition from 'history' to 'monitoring'
plot(bfm1)

dates <- as.Date(names(VI_m), "%Y-%m-%d")

# here we define the function that we will apply across the brick using the `app` function:
bfmRaster = function(pixels)
{
  tspx <- timeser(pixels, dates) # create a timeseries of all pixels
  bfm <- bfastmonitor(tspx, response ~ trend + harmon, order = 3, start = c(2019,1)) # run bfast on all pixels
  return(c(bfm$breakpoint, bfm$magnitude)) 
}

# apply the function to each raster cell (time series)
# Optionally you can supply an argument cores=n (where n is the number of cores on your computer) for a potential speed boost
bfmR <- app(VI_m, bfmRaster)
names(bfmR) <- c('time of break', 'magnitude of change')
plot(bfmR) # resulting time and magnitude of change