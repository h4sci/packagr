#' Create an inventory for time series
#'
#' \code{create_inventory} constructs an inventory. It calculates the mean,
#' standard deviation of the various series and outputs its frequency, type, and
#' name in a table.
#'
#' The function reads a list of variables that are of class time series. For
#' every list the function looks for the following things: \itemize{ \item
#' name() \item frequency() \item mean() \item sd() }
#'
#' @param flows A list of time series variable (flow variables).
#' @param stocks A list of time series variable (stock variables).
#'
#' @return A table
#' @import stats
#' @examples
#' inventory <- create_inventory(flows = data$flows, stocks = data$stocks)
#' @export
#'
create_inventory <- function(flows, stocks) {

  inventory <- rbind(data.frame("key" = as.character(names(flows)),
                                "type" = factor("flow", levels = c("stock","flow")),
                                "freq" = sapply(flows, frequency),
                                "mean" = sapply(flows, mean, na.rm=T),
                                "sd" = sapply(flows, sd, na.rm=T),
                                stringsAsFactors = F,
                                row.names = NULL),
                     data.frame("key" = names(stocks),
                                "type" = factor("stock", levels = c("stock","flow")),
                                "freq" = sapply(stocks, frequency),
                                "mean" = sapply(stocks, mean, na.rm=T),
                                "sd" = sapply(stocks, sd, na.rm=T),
                                stringsAsFactors = F,
                                row.names = NULL))

  # remove NULL entires
  if(length(which(inventory$key == "") > 0)) inventory <- inventory[-which(inventory$key == ""),]

  return(inventory)

}
