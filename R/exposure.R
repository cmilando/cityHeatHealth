#' Scripts to clean the exposure data
#'
#' things to be aware of
# - there may be missing data, so fill in with na.approx
# - you need to include the full year so you can do xbasis matrix
# - you may need to know the mapping so you can average up
#' @param data
#' @importFrom dlnm lag
#' @importFrom zoo na.approx
#' @importFrom lubridate month
#' @returns
#' @export
#'
#' @examples
make_exposure_matrix <- function(data) {

  # *************
  # THINGS TO FIX
  # - names -- need to have the mapping be for TOWN20, date, tmax_C
  # - dont hardcode which months are the summer months: 5-9
  # - don't use lubridate for month
  # - warning if the dates are not the whole year because then
  #   there will be NAs
  # - dont hardcode the number of lags
  # - don't hardcode the name of the lagcolumn
  # *************

  # need to fill in NA values
  exposure1_l <- split(data, f = exposure1$TOWN20)


  exposure2_l <- vector("list", length(exposure1_l))

  for(i in 1:length(exposure1_l)) {

    x <- exposure1_l[[i]]

    ev1 <- is.na(x$tmax_C)

    if(all(ev1)) next

    x$tmax_C  <- zoo::na.approx(x$tmax_C)
    if(any(is.na(x$tmax_C))) stop()

    # have to do this here so it doesn't bleed over
    # you can hard-code this so it doesn't look weird
    x$Templag1 <- dlnm::lag(x$tmax_C, 1)
    x$Templag2 <- dlnm::lag(x$tmax_C, 2)
    x$Templag3 <- dlnm::lag(x$tmax_C, 3)
    x$Templag4 <- dlnm::lag(x$tmax_C, 4)
    x$Templag5 <- dlnm::lag(x$tmax_C, 5)
    x$Templag6 <- dlnm::lag(x$tmax_C, 6)
    x$Templag7 <- dlnm::lag(x$tmax_C, 7)
    x$Templag8 <- dlnm::lag(x$tmax_C, 8)

    exposure2_l[[i]] <- x
  }

  exposure2 <- do.call(rbind, exposure2_l)

  head(exposure2)

  summer_exposure <- exposure2[lubridate::month(exposure2$date) %in% 5:9, ]

  class(summer_exposure) <- c(class(summer_exposure), "exposure")

  return(summer_exposure)

}
