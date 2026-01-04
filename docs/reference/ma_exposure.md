# Daily Maximum Temperature Exposure by Town

A SIMULATED dataset containing *fictional* daily maximum temperature
measurements linked to Massachusetts towns. Each row represents a single
day–town observation with the corresponding maximum temperature in
degrees Celsius.

## Format

A data frame with one row per date and town and the following variables:

- date:

  Date. Calendar date of the observation.

- tmax_C:

  Numeric. Daily maximum air temperature in degrees Celsius.

- TOWN20:

  Character. Name of the Massachusetts town associated with the
  temperature observation.

- COUNTY20:

  Character. Name of the Massachusetts county associated with TOWN20.

## Details

This dataset is typically used for environmental exposure assessment and
can be joined to town-level geographic or demographic datasets using
`TOWN20`.
