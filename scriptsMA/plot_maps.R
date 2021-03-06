library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

nc <- nc_open('ERA5_July2003_2mtemp.nc')

# load package
library(sp)
library(ncdf4)

# extract variable name, size and dimension
v <- nc$var[[1]]
size <- v$varsize
dims <- v$ndims
nt <- size[dims]              # length of time dimension
lat <- nc$dim$latitude$vals   # latitude position
lon <- nc$dim$longitude$vals  # longitude position

# read sst variable
r<-list()
for (i in 1:nt) {
  start <- rep(1,dims)     # begin with start=(1,1,...,1)
  start[dims] <- i             # change to start=(1,1,...,i) to read    timestep i
  count <- size                # begin with count=(nx,ny,...,nt), reads entire var
  count[dims] <- 1             # change to count=(nx,ny,...,1) to read 1 tstep
  
  dt<-ncvar_get(nc, varid = 't2m', start = start, count = count)
  
  # convert to raster
  r[i]<-raster(dt)
}

# create layer stack with time dimension
r<-stack(r)

# transpose the raster to have correct orientation
rt<-t(r)
extent(rt)<-extent(c(range(lon), range(lat)))


madridRome <- data.frame(row.names=1:2)
madridRome[1:2, 'x'] <- c((360-3.6919444), 12.4964)
madridRome[1:2, 'y'] <- c(40.41889, 41.9028)
coordinates(madridRome) <- ~ x + y
 
# plot the result
#spplot(rt)
windows(12,8)
spplot(rt, sp.layout = list("sp.points", modeldataxy, pch = 16, cex = 2, col = "black"))

