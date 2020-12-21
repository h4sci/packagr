
rm(list = ls())

# PACKAGES ----------------------------------------------------------------

library(timeseriesdb)
library(kofdata)
library(DatastreamDSWS2R)
library(zoo)
library(tsbox)
library(data.table)
library(openxlsx)
library(forecast)
library(tempdisagg)
library(tibbletime)
library(dplyr)
library(xts)

# FUNCTIONS ---------------------------------------------------------------

# datastream request function
ds2ts <- function(dskeys, wrn = TRUE, pwfile = "ds.txt"){
  
  library(DatastreamDSWS2R)
  creds <- unlist(strsplit(scan(pwfile, quiet = T, what = "")[3], ":", fixed = T))
  mydsws <- dsws$new(username = creds[4], password = creds[5])
  
  out <- lapply(dskeys, function(x) mydsws$timeSeriesRequest(instrument = x, 
                                                             startDate = "1960-01-01",
                                                             endDate = "3Y"))
  names(out) <- dskeys
  
  return(out)
}

daily2weekly <- function(x){
  
  idx <- plyr::round_any(x = as.numeric(format(time(x), "%Y")) + 
                           (as.numeric(format(time(x), "%m"))-1)/12 + 
                           as.numeric(format(time(x), "%d"))/365,
                         accuracy = 1/48,
                         f = floor)
  
  ts_weekly <- as.ts(aggregate(x = x,
                               by = idx,
                               FUN = mean,
                               na.rm=T))
  ts_weekly[is.nan(ts_weekly)] <- NA
  ts_weekly
  
}


# Import Metadata ---------------------------------------------------------

# get metadata & keys
metadata <- openxlsx::read.xlsx("code/in/Indicators_UK.xlsx", sheet = 6) # sheet 6: test dataset, sheet = 5 total dataset
metadata[,6] <- convertToDate(metadata[,6])
save(metadata, dat_final, file = "code/out/metadata_UK.Rda")



# IMPORT DATA: DATASTREAM -------------------------------------------------

inds <- seq(as.Date("1960-01-01"), as.Date(Sys.Date()), by = "day")

# (YOU NEED THE SAME ds.txt FILE AS IN RSTUDIO SERVER IN THE PROJECT DIRECTORY)
keys_ds <- metadata$keys
ts_ds <- ds2ts(keys_ds)

out_ts <- list()
ts_ds_adj <- lapply(names(ts_ds), function(ix){
  
  out_raw <- ts_ds[[ix]]
  freq <- metadata[which(metadata$keys == ix),"frequency"]
  
  if(freq == 4){
    
    out_ts <- ts(out_raw,
                 start = as.numeric(as.yearqtr(time(out_raw))[1]),
                 frequency = freq)
    
  } else if(freq == 12){

    out_ts <- ts(out_raw,
                 start = as.numeric(as.yearmon(time(out_raw))[1]),
                 frequency = freq)
    
  } else if(freq == 365){
    # Aggregate Daily to Monthly Variables
    ts <- xts(out_raw, as.Date(time(out_raw)), "%Y-%m-%d")
    ts_m = apply.monthly(ts, last)
    out_ts <- ts(ts_m,
                 start = as.numeric(as.yearmon(time(ts_m))[1]),
                 frequency = freq)
    
  # } else if(freq == 48){
  # 
  #   out_ts <- zoo(x = as.vector(out_raw), order.by = time(out_raw))
  #   idx <- plyr::round_any(x = as.numeric(format(time(out_raw), "%Y")) +
  #                            (as.numeric(format(time(out_raw), "%m"))-1)/12 +
  #                            as.numeric(format(time(out_raw), "%d"))/365,
  #                          accuracy = 1/48,
  #                          f = floor)
  #   out_ts <- as.ts(aggregate(x = out_ts,
  #                             by = idx,
  #                             FUN = mean,
  #                             na.rm=T))

  }
  
  return(out_ts)
  
}); names(ts_ds_adj) <- names(ts_ds)


# TRANSFORMATION ----------------------------------------------------------

# run appropriate transformations
dat <- c(ts_ds_adj)

# doesn't work with daily variables that are aggregated to monthly
dat_adj = lapply(metadata$keys, function(ix){
  # Transform according to order of integration
  if(metadata[which(metadata$keys == ix),"transform"] == "none"){

    out = dat[[ix]]

  } else if(metadata[which(metadata$keys == ix),"transform"] == "logdiff"){

    out = diff(log(dat[[ix]]))

  } else if(metadata[which(metadata$keys == ix),"transform"] == "yoydt"){

    x_adj <- diff(log(dat[[ix]]), lag = frequency(dat[[ix]]))
    out = x_adj - zoo::rollmeanr(x_adj, k = 3*frequency(dat[[ix]]), na.rm=T)

  } else if(metadata[which(metadata$keys == ix),"transform"] == "dt"){

    out = dat[[ix]] - zoo::rollmeanr(dat[[ix]], k = 3*frequency(dat[[ix]]), na.rm=T)

  }

  return(out)

}); names(dat_adj) <- metadata$keys

# start sample in 1990
dat_final <- lapply(dat_adj, function(x) na.trim(window(x, start = 1990)))


# group into stocks and flows
types <- sapply(names(dat_final), function(ix) metadata[which(metadata$keys == ix),"flow"])
data = list("flows" = dat_final[which(types == 1)],
           "stocks" = dat_final[which(types == 0)])

# save to file
save(data, dat_final, file = "code/out/data_UK.Rda")


