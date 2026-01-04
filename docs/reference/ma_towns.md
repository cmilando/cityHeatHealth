# Massachusetts Towns and Cities (2020 Census)

A dataset containing geographic identifiers, administrative
classifications, land and water area, population counts, housing units,
and related attributes for all cities and towns in Massachusetts,
primarily based on the 2020 U.S. Census and TIGER/MAF products. See
[link](https://www.mass.gov/info-details/massgis-data-2020-us-census-towns#downloads)

## Format

A shapefile with one row per Massachusetts city or town and the
following variables:

- STATEFP:

  Character. 2020 Census state FIPS code.

- COUNTYFP:

  Character. 2020 Census county FIPS code.

- COUSUBFP:

  Character. 2020 Census county subdivision FIPS code.

- COUSUBNS:

  Character. 2020 Census ANSI feature code for the county subdivision.

- GEOID:

  Character. 2020 Census county subdivision identifier (STATEFP +
  COUNTYFP + COUSUBFP).

- NAMELSAD:

  Character. Name and translated legal/statistical area description for
  the county subdivision.

- LSAD:

  Character. Legal/statistical area description code for the county
  subdivision.

- CLASSFP:

  Character. 2020 Census FIPS class code.

- MTFCC:

  Character. MAF/TIGER Feature Class Code (G4040).

- CNECTAFP:

  Character. 2020 Census combined New England city and town area code.

- NECTAFP:

  Character. 2020 Census New England city and town area code.

- NCTADVFP:

  Character. 2020 Census New England city and town area division code.

- FUNCSTAT:

  Character. 2020 Census functional status.

- ALAND:

  Numeric. 2020 Census land area in square meters.

- AWATER:

  Numeric. 2020 Census water area in square meters.

- INTPTLAT:

  Character. Latitude of the internal point (2020 Census).

- INTPTLON:

  Character. Longitude of the internal point (2020 Census).

- TOWN20:

  Character. Official name of the city or town.

- TOWN_ID:

  Integer. Town ID assigned to the official town name.

- FIPS_STCO:

  Integer. Federal Information Processing Standard (state/county) code.

- COUNTY20:

  Character. Name of the county in which the town is located.

- TYPE:

  Character. Municipality type: C = City; T = Town; TC = Town with city
  form of government.

- FOURCOLOR:

  Integer. Four-color shading code (1–4) such that no adjacent towns
  share the same code.

- AREA_ACRES:

  Numeric. Area of the town in acres.

- SQ_MILES:

  Numeric. Area of the town in square miles.

- Shape_Length:

  Numeric. Perimeter of the town in meters.

- Shape_Area:

  Numeric. Area of the town in square meters.

- POP1960:

  Integer. U.S. Census population in 1960.

- POP1970:

  Integer. U.S. Census population in 1970.

- POP1980:

  Integer. U.S. Census population in 1980.

- POP1990:

  Integer. U.S. Census population in 1990.

- POP2000:

  Integer. U.S. Census population in 2000.

- POP2010:

  Integer. U.S. Census population in 2010.

- POP2020:

  Integer. 2020 Census population count (PL 94-171 redistricting data).

- POPCH10_20:

  Integer. Change in population from 2010 to 2020.

- HOUSING20:

  Integer. 2020 Census count of housing units (PL 94-171).

## Source

U.S. Census Bureau, 2020 Decennial Census; TIGER/MAF geographic
products.

## Details

This dataset is typically used for environmental exposure assessment and
can be joined to town-level geographic or demographic datasets using
`TOWN20`.
