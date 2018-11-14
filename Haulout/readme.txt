

The following scripts were used to process, analyze, and plot haul-out data.  Several require the glmmLDTS package, which Jay has on github.

-get_haulout_records_db.R: Took data from database, did some minimal filtering based on dates, animals missing covariates, etc.
It produces a spatial points data frame with weather covariates attached to hourly haul-out records
("Haulout_SpPtsDF_16May2018.RData")

-Assemble_weather_covariates_grid.R:  This also pulls data from the database.  This script pulls covariates from the larger BOSS/CHESS grid
at 3 hour intervals and summarizes them for haulout predictions.  For the haul-out paper, the relevant data produced are "mean_covs_for_ho_yr_effects.Rdata"
which averages haul-out weather covariates over the Bering Sea for predicting mean haul-out behavior (this will be used in the "plot_haulout_species.R" script).
This script also outputs a data.frame of covariates needed for BOSS haulout predictions ("BOSS_weather_covs.RData"); the BOSS grid stored in the database is also output 
to facilitate matching of cell ID with spatial coordinates (not needed for the haul-out paper).

-process_haulout_data4.R:  This script takes the output from get_haulout_records_db.R and does some more filtering and data transforming before running
the analyses.  It also produces the map of haulout observations (Observations.pdf or Observations.jpeg) from the paper. Breakdowns of sample sizes in the paper are obtained
from the 'summaryBy' statements and also by the "HO_data_summary_hours.csv" and "HO_data_summary_nindiv.csv" files). Note that plot creation
could potentially be moved to an .Rmd file if the pts_sp object in "Haulout_SPPtsDF_post_filtering.RData" were loaded).  The primary task of this script is 
to run GLMPM models, main output is the
model results (in the various "test_" .RData files).  Note that the fixed effects and p-values, etc. are available in e.g. test.ribbon$typeIII.hypoth.  I *HAVEN'T*
updated appendices with these values from new model runs.

-plot_haulout_species.R:  This script takes model results and uses them to produce predictive plots of haul-out behavior.  Many of the plots are
not finalized until the end of the script so that e.g. results can be visualized across species.  This script also outputs a data.frame used
in the next script for stable stage modeling.

-compute_stable_stage.R: This script computes stable stage distributions, produces a plot of these, and then loads in analysis output and the
data frame from the 2nd and 3rd scripts to predict consequences of stable stage modeling, including cases where (a) age class is accounted for in
haul-out models and stable stage distributions are used to predict mean haul-out probabilities for the population, and (b) age class is ignored
in haul-out analysis and results are taken at face value.  It actually takes awhile to run because it needs to fit different haul-out models.  I'd advise caching.
