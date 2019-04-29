# Haulout: Export haulout data for Paul Conn's analysis
# S. Hardy, 01MAR2018

# Create functions -----------------------------------------------
# Function to install packages needed
install_pkg <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

# Install libraries ----------------------------------------------
install_pkg("RPostgreSQL")
install_pkg("sf")

# Run code -------------------------------------------------------
# Extract data from DB ------------------------------------------------------------------
con <- RPostgreSQL::dbConnect(PostgreSQL(), 
                              dbname = Sys.getenv("pep_db"), 
                              host = Sys.getenv("pep_ip"), 
                              user = Sys.getenv("pep_user"), 
                              rstudioapi::askForPassword(paste("Enter your DB password for user account: ", Sys.getenv("pep_user"), sep = "")))
#dat.sf <- sf::st_read_db(con, 
#                         query = "SELECT * FROM telem.res_iceseal_haulout_cov", 
#                         geom_column = "geom")

grid.sf <- sf::st_read_db(con, 
                         query = "SELECT * FROM base.geo_analysis_grid", 
                         geom_column = "geom")
save(grid.sf,file="BOSS_grid_from_db.Rdata")

### figure out cells to get average weather data for the "Bering Sea"
Coordinates = sf::st_coordinates(sf::st_centroid(grid.sf))
Which.inc = which(grid.sf$cell<15000)
#Which.inc = which(Coordinates[,"X"]>-300000 & grid.sf$cell<15000)
Cell.inc = grid.sf$cell[Which.inc]
plot(grid.sf[Which.inc,])
#grid.wx.Bering <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
#                                   WHERE cell>11824 AND cell<30542")
Cell_char = paste(Cell.inc,collapse=",")
Query_char = paste0("SELECT * FROM base.tbl_analysis_grid_cov_wx WHERE cell IN (",Cell_char,")")
grid.wx.Bering <- RPostgreSQL::dbGetQuery(con, Query_char)
#Which.inc.df = which(grid.wx.Bering$cell %in% Cell.inc)
#grid.wx.Bering = grid.wx.Bering[Which.inc.df,]
grid.wx.Bering$fdatetime_range_start=lubridate::with_tz(grid.wx.Bering$fdatetime_range_start,tz="UTC")
Hour = lubridate::hour(grid.wx.Bering$fdatetime_range_start)
#base on Hour = 23 UTC = 15 AKDST.  Convert to solar time later.
Which.23 = which(Hour %in% c(22,23))  #need both - daylight savings issue?
grid.wx.Bering = grid.wx.Bering[Which.23,]
#order by date-time
grid.wx.Bering = grid.wx.Bering[order(grid.wx.Bering$fdatetime_range_start),]
#get rid of 2004
Year = lubridate::year(grid.wx.Bering$fdatetime_range_start)
grid.wx.Bering=grid.wx.Bering[-which(Year==2004),]
#get rid of february
Month = lubridate::month(grid.wx.Bering$fdatetime_range_start)
grid.wx.Bering=grid.wx.Bering[-which(Month==2),]
Unique.dt = unique(grid.wx.Bering$fdatetime_range_start)
n.unique=length(Unique.dt)
Day = lubridate::day(grid.wx.Bering$fdatetime_range_start)
Year = lubridate::year(grid.wx.Bering$fdatetime_range_start)
Mean_covs = data.frame(Year=rep(0,n.unique),Day=rep(0,n.unique),
                       rast_acpcp=rep(0,n.unique),rast_air2m=rep(0,n.unique),rast_airsfc=rep(0,n.unique),
                       rast_prmsl=rep(0,n.unique),rast_uwnd=rep(0,n.unique),rast_vwnd=rep(0,n.unique))
for(i in 1:n.unique){
  Which = which(grid.wx.Bering$fdatetime_range_start==Unique.dt[i])
  Mean_covs[i,"Day"]=Day[Which][1]
  Mean_covs[i,"Year"]=Year[Which][1]
  Mean_covs[i,3:8]=colMeans(grid.wx.Bering[Which,3:8])
}
#standardize to format used in GLMPMs
Modeled_covs = data.frame(year = Mean_covs$Year,day=Mean_covs$Day,pressure=(Mean_covs$rast_prmsl-100000)/10000,
                          precip = Mean_covs$rast_acpcp, temp2=(Mean_covs$rast_air2m-270)/27, 
                          wind = sqrt(Mean_covs$rast_uwnd^2+Mean_covs$rast_vwnd^2)/10)
save(Modeled_covs,file="mean_covs_for_ho_yr_effects_15Apr2019.Rdata")

### BOSS covariate query

#2012: 4/10-5/8  (go one day ahead since times in DB are in UTC; note 5/10 will convert to 5/10 0:00)
grid.wx <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx 
                                   WHERE fdatetime_range_start BETWEEN '2012-04-10 00:00:00' AND '2012-05-10 00:00:00'")
grid.wx.modeled = data.frame(cell=grid.wx$cell,dt_UTC = grid.wx$fdatetime_range_start,pressure=(grid.wx$rast_prmsl-100000)/10000,
                          precip = grid.wx$rast_acpcp, temp2=(grid.wx$rast_air2m-270)/27, 
                          wind = sqrt(grid.wx$rast_uwnd^2+grid.wx$rast_vwnd^2)/10)
save(grid.wx.modeled,file="BOSS_weather_covs_15Apr2019.RData")
#grid.wx2 <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_wx LIMIT 10000")
#grid.seaice <- RPostgreSQL::dbGetQuery(con, "SELECT * FROM base.tbl_analysis_grid_cov_seaice --LIMIT 100")

# Once all are imported into memory, you should be able to join the covariates to the grid using the cell field

# In the above code, I added "--LIMIT 100". If you remove the --, this part of the query will no longer be commented out. 
#     The LIMIT ## is helpful for looking at a subset of the data before trying to get it all...

# The date/time field for the weather covariates is fdatetime_range_start...this field contains the first date/time of the 3-hour time range the covariates cover
# The date/time field for the sea ice covariate is fdate...there is a single record for each date