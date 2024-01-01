library(ncdf4)

nc <- nc_open("/work/bb0839/b380126/echam-6.3.01/ape-boundaries-sc/L192/T63_sst_aqua_Qobs.nc",
              write = TRUE)
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
sst <- ncvar_get(nc, "sst")

sst.pert <- array(rep(outer(lon, lat, function(lon, lat) -5 * exp(-0.5 * ((lon - 275)^2 + (lat + 20)^2) / 100)),
                      12),
                  dim(sst))
sst <- sst + sst.pert
ncvar_put(nc, "sst", sst)
nc_close(nc)
