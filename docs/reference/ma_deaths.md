# Daily Deaths by Town

A SIMULATED dataset containing *fictional* daily deaths linked to
Massachusetts towns. Each row represents a single day–town observation.

## Format

A data frame with one row per date and town and the following variables:

- date:

  Date. Calendar date of the observation.

- daily_deaths:

  Numeric. Daily deaths - SIMULATED!

- age_grp:

  Character. Either '0-17', '18-64' or '65+'

- sex:

  Character. either 'M' or 'F' in this example

- TOWN20:

  Character. Name of the Massachusetts town associated with the
  temperature observation.

- COUNTY20:

  Character. Name of the Massachusetts county associated with TOWN20.

## Details

This dataset is typically used for environmental exposure assessment and
can be joined to town-level geographic or demographic datasets using
`TOWN20`.
