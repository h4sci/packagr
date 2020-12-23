#' Prepare time series data for usage
#'
#' \code{prepare_data} transforms the time series into a usable format for time
#' series analysis.
#'
#' The function conducts a normalization by demeaning and dividing by the
#' standard deviation of the individual time series.
#'
#' @param inventory Inventory of the corresponding data.
#' @param target Defines the target variable.
#' @inheritParams create_inventory
#' @return A dataframe
#' @import stats
#' @importFrom dplyr left_join
#' @importFrom zoo na.trim
#' @examples
#' Xmat <- prepare_data(flows = data$flows, stocks = data$stocks, inventory, target = target)
#' @export
#'
prepare_data <- function(flows, stocks, inventory, target){

  data <- c(flows, stocks)

  # standardize data
  data_std <- lapply(inventory$key, function(ix){

    (data[[ix]] - inventory[which(inventory$key == ix),"mean"])/
      inventory[which(inventory$key == ix),"sd"]

  }); names(data_std) <- inventory$key

  # adjust time periods
  freq_max <- max(sapply(data, frequency))

  data_adj <- lapply(inventory$key, function(ix){

    tim <- time(data_std[[ix]]) + (freq_max/frequency(data_std[[ix]]) - 1)/freq_max
    out <- as.data.frame(cbind(round(tim,5), data_std[[ix]]))
    colnames(out) <- c("time",ix)

    return(out)

  }); names(data_adj) <- inventory$key

  # coerce to sparse matrix
  idx = data.frame("time" = round(seq(1900, 2100, 1/freq_max),5))

  dfs <- lapply(data_adj, function(fx){

    tab <- left_join(idx, fx, by = "time")
    tab[,-1]

  })

  dfs <- do.call(cbind,dfs)
  colnames(dfs) <- inventory$key
  dfts <- ts(dfs, start = 1900, frequency = max(sapply(data, function(x) frequency(x))))
  dfts <- na.trim(dfts, is.na = "all")

  dfts[which(is.na(dfts))] <- 0

  return(dfts)

}
