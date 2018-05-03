

The following scripts were used to process, analyze, and plot haul-out data:

-Haulout_Covariates_ForPaul20180301_SKH.R: Took data from database, did some minimal filtering based on dates, animals missing covariates, etc.
note that the database just changed on 5/2 after importing covariate values that increased the number of records somewhat.  I have not updated
analyses to reflect any new records.  Someone else is welcome to, although the table of sample sizes, etc. will need to be updated by hand, model
selection results changed, etc.  I'm not expecting a huge change so don't feel it is worth the effort to update the analysis following the
database change.  For this reason, I suggest using the output from this script as run in March 2018 for analysis input
("Haulout_SpPtsDF_29Mar2018.RData")

-process_haulout_data3.R:  This script takes the output from the previous step and does some more filtering and data transforming before running
the analyses.  It also produces the map of haulout observations (Observations.pdf or Observations.jpeg) from the paper.  Note that plot creation
could potentially be moved to an .Rmd file if the pts_sp object in "Haulout_SPPtsDF_post_filtering.RData" were loaded).  The main output is the
model results.

-plot_haulout_species.R:  This script takes model results and uses them to produce predictive plots of haul-out behavior.  Many of the plots are
not finalized until the end of the script so that e.g. results can be visualized across species.  This script also outputs a data.frame for used
in the next script for stable stage modeling.

-compute_stable_stage.R: This script computes stable stage distributions, produces a plot of these, and then loads in analysis output and the
data frame from the 2nd and 3rd scripts to predict consequences of stable stage modeling, including cases where (a) age class is accounted for in
haul-out models and stable stage distributions are used to predict mean haul-out probabilities for the population, and (b) age class is ignored
in haul-out analysis and results are taken at face value.
