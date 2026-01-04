# Massachusetts Towns and Cities (2020 Census)

A dataset containing geographic identifiers, administrative
classifications, land and water area, population counts, housing units,
and related attributes for all counties in Massachusetts. See
[link](https://www.mass.gov/info-details/massgis-data-counties#downloads)

## Format

A shapefile with one row per Massachusetts city or town and the
following variables:

- FIPSID:

  Character

- COUNTY20:

  Character. 2020 Census county FIPS code.

- AREA_ACRES:

  Numeric.

- OBJECTID:

  Numeric.

- SHAPE_AREA:

  Numeric.

- SHAPE_LEN:

  Numeric.

- geometry:

  Numeric.

## Source

U.S. Census Bureau, 2020 Decennial Census; TIGER/MAF geographic
products.

## Details

This dataset is typically used for environmental exposure assessment and
can be joined to town-level geographic or demographic datasets using
`COUNTY20`.
