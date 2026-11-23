# Heatwaves and violent offenses in US cities

R code for the main analyses of *Disproportionate burden of heatwave-related violent crime across US cities in past and future climates* (in preparation).

`main_analysis.R` runs the three steps of the paper in order: the case-crossover association between heatwaves and violent offenses by social vulnerability group, the burden attributable to anthropogenic warming in 2005–2022, and the projected burden under four SSP scenarios. Results are written to `output/`.

`example_data.rds` contains simulated data for four hypothetical cities so that the script runs as is (`Rscript main_analysis.R`, under a minute). These are not study data. The study data are public and their sources are listed in the Data availability statement of the paper. 

Requirements: R 4.6.0 with data.table, dlnm 2.4.10, fixest 0.14.2, Rssa 1.1 and MASS.

License: MIT.
