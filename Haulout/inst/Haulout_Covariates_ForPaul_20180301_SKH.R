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
dat.sf <- sf::st_read_db(con, 
                         query = "SELECT * FROM telem.res_iceseal_haulout_cov", 
                         geom_column = "geom")

#remove records between July and February
dat.df = data.frame(dat.sf)
Month = months(dat.df[,"haulout_dt"])
Which = which(Month %in% c("March","April","May","June"))
dat.sf = dat.sf[Which,]
dat.df = data.frame(dat.sf)

#remove records without covariate information
dat.sf = dat.sf[-which(is.na(dat.df[,"rast_vwnd"])),]
dat.df = data.frame(dat.sf)

bering_base <- nPacMaps::bering()

#take a look at  ringed seal data
ringed.sf = dat.sf[which(dat.df[,"species"]=="Ph"),]
ringed.df = data.frame(ringed.sf)
sex.age = apply(ringed.df[,c("sex","age")],1,'paste',collapse='.')
Years = factor(as.POSIXlt(ringed.df[,"haulout_dt"])$year)
tabulate(Years)
my.plot = ggplot2::ggplot()+ggplot2::geom_sf(data=ringed.sf)
my.plot = my.plot + ggplot2::geom_sf(data = bering_base,fill = "grey60", size = 0.2)


#take out ringed seals for main haulout paper
dat.sf = dat.sf[-which(dat.df[,"species"]=="Ph"),]
dat.df = data.frame(dat.sf)

#quick check to make sure points in reasonable places
Which.plot = c(1:13414)*10
my.plot = ggplot2::ggplot()+ggplot2::geom_sf(data=dat.sf[Which.plot,])
my.plot = my.plot + ggplot2::geom_sf(data = bering_base,fill = "grey60", size = 0.2)
my.plot

# change to spatial points for prior compatibility
pts_sp <- dat.sf %>% as("Spatial")
save(pts_sp,file = "Haulout_SpPtsDF_29Mar2018.RData")

