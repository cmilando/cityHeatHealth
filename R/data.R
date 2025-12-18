#' Massachusetts Towns and Cities (2020 Census)
#'
#' A dataset containing geographic identifiers, administrative classifications,
#' land and water area, population counts, housing units, and related attributes
#' for all cities and towns in Massachusetts, primarily based on the 2020 U.S.
#' Census and TIGER/MAF products. See [link](https://www.mass.gov/info-details/massgis-data-2020-us-census-towns#downloads)
#'
#' @format A shapefile with one row per Massachusetts city or town and the
#' following variables:
#' \describe{
#'   \item{STATEFP}{Character. 2020 Census state FIPS code.}
#'   \item{COUNTYFP}{Character. 2020 Census county FIPS code.}
#'   \item{COUSUBFP}{Character. 2020 Census county subdivision FIPS code.}
#'   \item{COUSUBNS}{Character. 2020 Census ANSI feature code for the county subdivision.}
#'   \item{GEOID}{Character. 2020 Census county subdivision identifier (STATEFP + COUNTYFP + COUSUBFP).}
#'   \item{NAMELSAD}{Character. Name and translated legal/statistical area description for the county subdivision.}
#'   \item{LSAD}{Character. Legal/statistical area description code for the county subdivision.}
#'   \item{CLASSFP}{Character. 2020 Census FIPS class code.}
#'   \item{MTFCC}{Character. MAF/TIGER Feature Class Code (G4040).}
#'   \item{CNECTAFP}{Character. 2020 Census combined New England city and town area code.}
#'   \item{NECTAFP}{Character. 2020 Census New England city and town area code.}
#'   \item{NCTADVFP}{Character. 2020 Census New England city and town area division code.}
#'   \item{FUNCSTAT}{Character. 2020 Census functional status.}
#'   \item{ALAND}{Numeric. 2020 Census land area in square meters.}
#'   \item{AWATER}{Numeric. 2020 Census water area in square meters.}
#'   \item{INTPTLAT}{Character. Latitude of the internal point (2020 Census).}
#'   \item{INTPTLON}{Character. Longitude of the internal point (2020 Census).}
#'   \item{TOWN}{Character. Official name of the city or town.}
#'   \item{TOWN_ID}{Integer. Town ID assigned to the official town name.}
#'   \item{FIPS_STCO}{Integer. Federal Information Processing Standard (state/county) code.}
#'   \item{COUNTY}{Character. Name of the county in which the town is located.}
#'   \item{TYPE}{Character. Municipality type: C = City; T = Town; TC = Town with city form of government.}
#'   \item{FOURCOLOR}{Integer. Four-color shading code (1–4) such that no adjacent towns share the same code.}
#'   \item{AREA_ACRES}{Numeric. Area of the town in acres.}
#'   \item{SQ_MILES}{Numeric. Area of the town in square miles.}
#'   \item{Shape_Length}{Numeric. Perimeter of the town in meters.}
#'   \item{Shape_Area}{Numeric. Area of the town in square meters.}
#'   \item{POP1960}{Integer. U.S. Census population in 1960.}
#'   \item{POP1970}{Integer. U.S. Census population in 1970.}
#'   \item{POP1980}{Integer. U.S. Census population in 1980.}
#'   \item{POP1990}{Integer. U.S. Census population in 1990.}
#'   \item{POP2000}{Integer. U.S. Census population in 2000.}
#'   \item{POP2010}{Integer. U.S. Census population in 2010.}
#'   \item{POP2020}{Integer. 2020 Census population count (PL 94-171 redistricting data).}
#'   \item{POPCH10_20}{Integer. Change in population from 2010 to 2020.}
#'   \item{HOUSING20}{Integer. 2020 Census count of housing units (PL 94-171).}
#' }
#'
#' @source U.S. Census Bureau, 2020 Decennial Census; TIGER/MAF geographic products.
#'
#' @name ma_towns
#' @docType data
#' @keywords datasets
NULL


#' Daily Maximum Temperature Exposure by Town
#'
#' A dataset containing _fictional_ daily maximum temperature measurements linked to
#' Massachusetts towns. Each row represents a single day–town observation
#' with the corresponding maximum temperature in degrees Celsius.
#'
#' @format A data frame with one row per date and town and the following variables:
#' \describe{
#'   \item{date}{Date. Calendar date of the observation.}
#'   \item{tmax_C}{Numeric. Daily maximum air temperature in degrees Celsius.}
#'   \item{TOWN20}{Character. Name of the Massachusetts town associated with the temperature observation.}
#' }
#'
#' @details
#' This dataset is typically used for environmental exposure assessment and
#' can be joined to town-level geographic or demographic datasets using
#' \code{TOWN20}.
#'
#' @name ma_exposure
#' @docType data
#' @keywords datasets
NULL
