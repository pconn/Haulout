#process and run a glmmLDTS model on haul-out data; data input as SpPtsDF (see 'Haulout_Covariates_For_Paul....R')

library(lubridate)  #some date time functionality
library(nPacMaps)
library(glmmLDTS)
library(solaR)  #for getting hour of day in solar time
library(splines)
library(sf)

load("Haulout_SpPtsDF_1Nov2018.RData")

pts_sp <- pts_sp[order(pts_sp@data[,"speno"],pts_sp@data[,"haulout_dt"]),]

#Date formatting
pts_sp@data$hour = hour(pts_sp[["haulout_dt"]])
pts_sp@data$month = month(pts_sp[["haulout_dt"]])
pts_sp@data$year = year(pts_sp[["haulout_dt"]])
pts_sp@data$day = day(pts_sp[["haulout_dt"]])
pts_sp@data$jday = yday(pts_sp[["haulout_dt"]])

#format and increment ages if necessary - need to do this before filtering so that pups tagged in fall don't
#show up as pups in spring
get_tag_year <- function(ID){
  #different formatting constructs based on who did the tagging
  if(nchar(ID)==14)yr=2000+as.numeric(substr(ID,3,4))
  else{
    if(substr(ID,1,4)=="KOTE")yr=2000+as.numeric(substr(ID,10,11))
    else yr = as.numeric(substr(ID,3,6))
  }
  yr
}
#format ages into a few unique categories
uniqueIDs <- unique(pts_sp[["speno"]])
n.unique=length(uniqueIDs)
pts_sp@data[which(pts_sp@data$age=="PUP"),"age"]="YOY"
pts_sp@data[which(pts_sp@data$age=="YRL"),"age"]="SUB"
for(iid in 1:n.unique){
  CurDat = pts_sp@data[which(pts_sp[["speno"]]==uniqueIDs[iid]),]
  Years = unique(CurDat[,"year"])
  if(CurDat[1,"age"]=="YOY" & length(Years)>1){
    Which.later = which(CurDat$year>Years[1])
    CurDat[Which.later,"age"]="SUB"
  }
  #still some animals tagged as YOY previous year that first enter the dataset as subadults
  tag.year = get_tag_year(uniqueIDs[iid])
  if(CurDat[1,"age"]=="YOY" & min(Years)>tag.year){
    CurDat[,"age"]="SUB"
  }
  pts_sp@data[which(pts_sp[["speno"]]==uniqueIDs[iid]),]=CurDat
}

#remove points on land
pts_sf = st_as_sf(pts_sp)
library(nPacMaps)
us_arctic_base = nPacMaps::bering()
arctic_union = st_combine(us_arctic_base)
Intersects = st_intersects(pts_sf,arctic_union)
Which.intersects = which(unlist(lapply(Intersects, any))>0)
pts_sf = pts_sf[-Which.intersects,]  #takes out another 2160
pts_sp = pts_sf %>% as("Spatial")

###end filtering
#save SpPtsDF in current state for use in reconstructing
save(pts_sp,file="Haulout_SpPtsDF_post_filtering.RData")


### plot observations map
library(nPacMaps)
us_arctic_base = nPacMaps::bering()
pts_sf = st_as_sf(pts_sp)
b <- sf::st_bbox(pts_sf)
b[c(1,2)]=b[c(1,2)]-100000
b[c(3,4)]=b[c(3,4)]+100000
pts_sf$Species=pts_sf$species
pts_sf$Species[which(pts_sf$Species=="Eb")]="Bearded"
pts_sf$Species[which(pts_sf$Species=="Hf")]="Ribbon"
pts_sf$Species[which(pts_sf$Species=="Pl")]="Spotted"
pts_plot = pts_sf[c(1:13205)*10,]  #just plot every 10 observations
library(ggplot2)
crap = ggplot() + geom_sf(data=us_arctic_base,fill="grey60",size=0.2)
crap = crap + geom_sf(data=pts_plot,size=0.2,shape=20,alpha=.1,aes(col=Species))
crap = crap + coord_sf(xlim = c(b["xmin"], b["xmax"]), ylim = c(b["ymin"], b["ymax"]))
crap = crap + scale_color_manual(values=c("brown", "blue", "darkgreen"))
crap = crap + guides(col = guide_legend(override.aes = list(shape = 20,size=1)))
pdf("Observations.pdf")
crap
dev.off()

jpeg("Observations.jpeg")
crap
dev.off()

#pts_sp = pts_sp[which(Coords[,])]

#formulate "age.sex" covariate - 'yoy', 'sub', 'adult f', 'adult m'
#also, we have trouble estimating age.sex * DOY interaction for YOY so set YOY = ADULT.F for estimation of interactions
pts_sp@data[which(pts_sp[["sex"]]=="female"),"sex"]="F"
pts_sp@data[which(pts_sp[["sex"]]=="male"),"sex"]="M"

pts_sp@data[,"age.sex"]='YOY'
pts_sp@data[which(pts_sp[["age"]]=="SUB"),"age.sex"]="SUB"
pts_sp@data[which(pts_sp[["age"]]=="ADT" & pts_sp[["sex"]]=="M"),"age.sex"]="ADULT.M"
pts_sp@data[which(pts_sp[["age"]]=="ADT" & pts_sp[["sex"]]=="F"),"age.sex"]="ADULT.F"
pts_sp@data[,"age.sex.inter"]=pts_sp[["age.sex"]]
pts_sp@data[which(pts_sp[["age.sex.inter"]]=="YOY"),"age.sex.inter"]="ADULT.F"
pts_sp[["age.sex"]]=as.factor(pts_sp[["age.sex"]])
pts_sp[["age.sex.inter"]]=as.factor(pts_sp[["age.sex.inter"]])

uniqueIDs <- unique(pts_sp[["speno"]])
n.unique=length(uniqueIDs)

StartEnd=data.frame(Start=rep(pts_sp[["haulout_dt"]][1],n.unique))  #determine starting and ending date for each deployment (see if age needs to be altered)
StartEnd$End=StartEnd$Start
for(i in 1:n.unique){
  CurDat = pts_sp@data[which(pts_sp[["speno"]]==uniqueIDs[i]),]
  StartEnd$Start[i]=CurDat[1,"haulout_dt"]
  StartEnd$End[i]=CurDat[nrow(CurDat),"haulout_dt"]
}
StartEnd$End-StartEnd$Start
cat(paste("Maximum # days in deployment: ",max(StartEnd$End-StartEnd$Start)/24,"\n")) #almost 8 years?!  So need to adjust ages 

#pts_sf = pts_sf[-which(pts_sf$speno=="HF2005_5898" & pts_sf$year>2006),]   #this long deployment likely a tag on a beach somewhere
#pts_sp = pts_sp[-which(pts_sp[["speno"]]=="HF2005_5898" & pts_sp[["year"]]>2006),]

#take a look at breakdown of age - sex - species
DF.unique = pts_sp@data[1,]
for(iunique in 2:n.unique)DF.unique[iunique,]=pts_sp@data[which(pts_sp[["speno"]]==uniqueIDs[iunique])[1],]
library(doBy)
summaryBy(year~species+sex+age,data=DF.unique,FUN=length)
summaryBy(year~species+sex+age,data=pts_sp@data,FUN=length)
Summary=summaryBy(percent_dry~species+year+month+age.sex,data=pts_sp@data,FUN=length)

#summaryBy(sex~year+species,data=DF.unique,FUN=length)

Year = as.character(c(2005:2017))
Month = as.character(c(3:6))
Species = c("Eb","Hf","Pl")
Age.sex = c("YOY","SUB","ADULT.F","ADULT.M")
#plot number of seals and hourly records by species, year, and month
Sum.table = matrix(0,13*5,12)
for(iyear in 1:13){
  for(imonth in 1:4){
    for(isp in 1:3){
      for(iage in 1:4){
        cur.records = 0
        cur.which = which(Summary[,"species"]==Species[isp] & Summary[,"month"]==Month[imonth] & Summary[,"year"]==Year[iyear] & Summary[,"age.sex"]==Age.sex[iage])
        if(length(cur.which)>0)Sum.table[5*(iyear-1)+1+imonth,(isp-1)*4+iage]=Summary[cur.which,"percent_dry.length"]
      }
    }
  }
  Sum.table[5*(iyear-1)+1,]=colSums(Sum.table[5*(iyear-1)+c(2:5),])
}
Tots = matrix(0,5,12)
for(irow in 1:5)Tots[irow,]=colSums(Sum.table[(c(1:12)-1)*5+irow,])
write.csv(rbind(Sum.table,Tots), file="HO_data_summary_hours.csv")

pts_sp@data$IDyear = paste(pts_sp[["speno"]],pts_sp[["year"]])


#total number of seals and seal-years
NoDupe=pts_sp@data
if(length(which(duplicated(pts_sp@data[,"IDyear"])==0))>0)NoDupe = pts_sp@data[which(duplicated(pts_sp@data[,"IDyear"])==0),]
Summary=summaryBy(IDyear~species+year+age.sex,data=NoDupe,FUN=length)
Sum.table = matrix(0,13,12)
for(iyear in 1:13){
  for(isp in 1:3){
    for(iage in 1:4){
      cur.records = 0
      cur.which = which(Summary[,"species"]==Species[isp] & Summary[,"year"]==Year[iyear] & Summary[,"age.sex"]==Age.sex[iage])
      if(length(cur.which)>0)Sum.table[iyear,(isp-1)*4+iage]=Summary[cur.which,"IDyear.length"]
    }
  }
}
write.csv(Sum.table,file="HO_data_summary_nindiv.csv")

Unique = pts_sp@data[which(duplicated(pts_sp@data[,"speno"])==0),]
summaryBy(speno~species+age.sex,data=Unique,FUN=length)
       
#HO@data[,"datadatetime"]=with_tz(HO@data[,"datadatetime"],tzone="GMT")


#check to see if there are any duplicated entries for a given individual (looks okay)
ID.time=paste(pts_sp@data[,"speno"],pts_sp@data[,"haulout_dt"])
Dupes=duplicated(ID.time)
if(sum(Dupes)>0)print("Warning: duplicated entries")



#enumerate individuals
ID=pts_sp[["IDyear"]]
unique.ID=unique(ID)
n.indiv=length(unique.ID)
#give a new ID, going from 1 to # unique "individuals"
ID.num=ID
for(i in 1:n.indiv){
  ID.num[which(ID.num==unique.ID[i])]=i
}
ID.num=as.integer(ID.num)
pts_sp@data[,"ID"]=ID.num

#change response to binary
pts_sp@data[,"Dry"]=round(as.numeric(pts_sp@data[,"percent_dry"])/100)


#produce new ID variable for AR1 component (time gaps lead to new IDs even within individual)
AR1.ID=rep(0,length(ID.num))
HourDay = pts_sp[["jday"]]*24+pts_sp[["hour"]]  #discontinuities in HourDay will trigger new ID
AR1.ID[1]=1
for(irec in 2:length(AR1.ID)){
  if(pts_sp@data[irec,"ID"]!=pts_sp@data[irec-1,"ID"] | HourDay[irec]!=(HourDay[irec-1]+1)){
    AR1.ID[irec]=AR1.ID[irec-1]+1
  }
  else AR1.ID[irec] = AR1.ID[irec-1]
}
pts_sp@data$AR1.ID = AR1.ID
n.AR1=max(AR1.ID)


#more formatting for glmmLDTS
#add time vector for AR(1) modeling
pts_sp@data[,"time"]=rep(0,nrow(pts_sp@data))
n.id=max(AR1.ID)
Bout.length = rep(0,n.id)
for(iID in 1:n.id){
  Cur.which=which(AR1.ID==iID)
  Bout.length[iID]=length(Cur.which)
  Cur.time=c(0:(Bout.length[iID]-1))
  pts_sp@data[Cur.which,"time"]=Cur.time
}



#add indicator for whether observation occurs N of Bering Strait
pts_sp[["Chukchi"]]=0
pts_sp@data[which(pts_sp[["Northing"]]<0.85),"Chukchi"]=1


#get hour in solar time
Longs = coordinates(spTransform(pts_sp, CRS("+proj=longlat +datum=WGS84")))[,1]
SolarT = local2Solar(pts_sp@data$haulout_dt,lon=Longs)
pts_sp[["hour"]]=hour(SolarT)
#fourier hour terms
pts_sp@data$sin1= sin(pi*pts_sp[["hour"]]/12) 
pts_sp@data$cos1= cos(pi*pts_sp[["hour"]]/12)
pts_sp@data$sin2= sin(pi*pts_sp[["hour"]]/6)
pts_sp@data$cos2= cos(pi*pts_sp[["hour"]]/6)
pts_sp@data$sin3= sin(pi*pts_sp[["hour"]]/4)
pts_sp@data$cos3= cos(pi*pts_sp[["hour"]]/4)

pts_sp@data[,"day"]=(pts_sp@data[,"jday"]-120)/10
pts_sp@data[,"day2"]=pts_sp@data[,"day"]^2
pts_sp@data[,"day3"]=pts_sp@data[,"day"]^3
pts_sp@data[,"ID"]=as.factor(pts_sp@data[,"ID"])
#pts_sp@data[,"Easting"]=pts_sp@coords[,1]
pts_sp@data[,"Northing"]=pts_sp@coords[,2]/mean(pts_sp@coords[,2])
#pts_sp@data[,"iceconc2"]=(pts_sp@data[,"iceconc"])^2
pts_sp@data[,"species"]=as.factor(pts_sp@data[,"species"])
#pts_sp@data[,"age"]=as.factor(pts_sp@data[,"age"])
#pts_sp@data[,"yoy"]=(pts_sp@data$age=="YOY")  #in case we want to model pup/non-pup as binary
#pts_sp@data[,"sex"]=as.factor(pts_sp@data[,"sex"])
pts_sp@data[,"pressure"]=(pts_sp@data[,"rast_prmsl"]-100000)/10000
pts_sp@data[,"precip"]=pts_sp@data[,"rast_acpcp"]
pts_sp@data[,"temp0"]=(pts_sp@data[,"rast_airsfc"]-270)/27
pts_sp@data[,"temp2"]=(pts_sp@data[,"rast_air2m"]-270)/27
#pts_sp@data[,"yoy"]=(pts_sp@data$age=="YOY")
pts_sp@data[,"wind"]=sqrt(pts_sp@data$rast_uwnd^2+pts_sp@data$rast_vwnd^2)/10 #meters per second
pts_sp@data$year=as.factor(pts_sp@data$year)
pts_sp@data[,"hour"]=factor(pts_sp@data[,"hour"])


#set sea ice concentration for values that fall on land


#sample analysis, removing all zero histories
#undebug(glmmLDTS1)

#save data used in modeling
HO_df = pts_sp@data
save(HO_df,file="HO_df.RData")

# #sample analysis
# test.mod <- glm(Dry ~ species + age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + year + temp2 + wind + Northing + pressure + wind*temp2 +
#                   species*sin1*day + species*cos1*day + species*sin2*day + species*cos2*day + species*sin3*day + species*cos3*day +
#                   species*sin1*day2 + species*cos1*day2 + species*sin2*day2 + species*cos2*day2 + species*sin3*day2 + species*cos3*day2 +  
#                   species*sin1*day3 + species*cos1*day3 + species*sin2*day3 + species*cos2*day3 + species*sin3*day3 + species*cos3*day3 +
#                   age.sex*day + age.sex*day2 + age.sex*day3 + year*day + year*day2 +
#                   temp2*day + temp2*day2 + temp2*day3 +
#                   Northing*day + Northing*day2,family='binomial',data=HO_df)
# ###looks like year effects are poorly estimated.  get rid of them for now.
# test.mod <- glm(Dry ~ species + age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + Northing + pressure + wind*temp2 +
#                   species*sin1*day + species*cos1*day + species*sin2*day + species*cos2*day + species*sin3*day + species*cos3*day +
#                   species*sin1*day2 + species*cos1*day2 + species*sin2*day2 + species*cos2*day2 + species*sin3*day2 + species*cos3*day2 +  
#                   species*sin1*day3 + species*cos1*day3 + species*sin2*day3 + species*cos2*day3 + species*sin3*day3 + species*cos3*day3 +
#                   age.sex*day + age.sex*day2 + age.sex*day3 +
#                   temp2*day + temp2*day2 + temp2*day3 +
#                   Northing*day + Northing*day2,family='binomial',data=HO_df)
# 
# 
# ###age.sex = YOY has issues w/ DOY interaction; model age.sex interactions = for adult females and YOY
# test.mod <- glm(Dry ~ species + age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + Northing + pressure + wind*temp2 +
#                   species*sin1*day + species*cos1*day + species*sin2*day + species*cos2*day + species*sin3*day + species*cos3*day +
#                   species*sin1*day2 + species*cos1*day2 + species*sin2*day2 + species*cos2*day2 + species*sin3*day2 + species*cos3*day2 +  
#                   species*sin1*day3 + species*cos1*day3 + species*sin2*day3 + species*cos2*day3 + species*sin3*day3 + species*cos3*day3 +
#                   age.sex.inter:day + age.sex.inter:day2 + age.sex.inter:day3 +
#                   temp2*day + temp2*day2 + temp2*day3 +
#                   Northing*day + Northing*day2,family='binomial',data=HO_df)
# 
# 
# 
# Pred.dat=HO_df[1:24,]
# Pred.dat[,"temp2"]=mean(HO_df[,"temp2"])
# Pred.dat[,"wind"]=mean(HO_df[,"wind"])
# Pred.dat[,"pressure"]=mean(HO_df[,"pressure"])
# Time = c(0:23)
# Pred.dat[,"sin1"]=sin(pi*Time/12)
# Pred.dat[,"cos1"]=cos(pi*Time/12)
# Pred.dat[,"sin2"]=sin(pi*Time/6)
# Pred.dat[,"cos2"]=cos(pi*Time/6)
# Pred.dat[,"sin3"]=sin(pi*Time/4)
# Pred.dat[,"cos3"]=cos(pi*Time/4)
# 
# plot(predict(test.mod,newdata=Pred.dat))


# 
# test.mod <- glmmLDTS(fixed.formula = Dry ~ species + age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + Northing + pressure + wind*temp2 + #year +
#                        species*sin1*day + species*cos1*day + species*sin2*day + species*cos2*day + species*sin3*day + species*cos3*day +
#                        species*sin1*day2 + species*cos1*day2 + species*sin2*day2 + species*cos2*day2 + species*sin3*day2 + species*cos3*day2 +  
#                        #species*sin1*day3 + species*cos1*day3 + species*sin2*day3 + species*cos2*day3 + species*sin3*day3 + species*cos3*day3 +
#                        age.sex.inter:day + age.sex.inter:day2 + age.sex.inter:day3 + #year*day + year*day2 +
#                        temp2*day + temp2*day2 + temp2*day3 +
#                        Northing*day + Northing*day2,
#                        random.formula = Dry ~ ID,
#                        data = HO_df,
#                        EstMeth="REML",
#                        timecol = "time", 
#                        #ridge.reg = "global",
#                        #lambda = 0.5,
#                        group.vec = "AR1.ID")

#save(test.mod,file="test_big.RData")

HO_ribbon = HO_df[HO_df$species=="Hf",]
test.ribbon <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + pressure + precip + wind*temp2 +
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                          age.sex:day + age.sex:day2 + age.sex:day3,
                        random.formula = Dry ~ speno,
                        data = HO_ribbon,
                        EstMeth="REML",
                        timecol = "time", 
                        #ridge.reg = "global",
                        #lambda = 0.5,
                        group.vec = "AR1.ID")

save(test.ribbon,file="test_ribbon.RData")

HO_spotted = HO_df[HO_df$species=="Pl",]
HO_spotted = HO_spotted[-which(HO_spotted$year=="2012"),]  #get rid of 2012 - only 1 seal and it causes numerical problems if included
test.spotted <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + pressure + precip + wind*temp2 +
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 + 
                          age.sex:day + age.sex:day2 + age.sex:day3,
                        random.formula = Dry ~ speno,
                        data = HO_spotted,
                        EstMeth="REML",
                        timecol = "time", 
                        #ridge.reg = "global",
                        #lambda = 0.5,
                        group.vec = "AR1.ID")
save(test.spotted,file="test_spotted.RData")

HO_bearded = HO_df[HO_df$species=="Eb",]
HO_bearded$age.sex = factor(HO_bearded$age.sex) #get rid of YOY
test.bearded <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + Northing + temp2 + wind + pressure + precip + wind*temp2 +
                           sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                           sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2+  
                           Northing*day + Northing*day2,
                         random.formula = Dry ~ speno,
                         data = HO_bearded,
                         EstMeth="REML", #"ML",  #REML was default
                         timecol = "time", 
                         #ridge.reg = "global",
                         #lambda = 0.5,
                         group.vec = "AR1.ID")
save(test.bearded,file="test_bearded.RData")
#crap = inla(Dry~day+day2+day3+f(IDyear,model="iid"),family="binomial",data=HO_bearded)

HO_ribbon$year = factor(HO_ribbon$year)  #no 2013 in data
test.ribbon.year <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + pressure + 
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                          age.sex:day + age.sex:day2 + age.sex:day3 + year:day + year:day2,
                        random.formula = Dry ~ speno,
                        data = HO_ribbon,
                        EstMeth="REML",
                        timecol = "time", 
                        #ridge.reg = "global",
                        #lambda = 0.5,
                        group.vec = "AR1.ID")

save(test.ribbon.year,file="test_ribbon_year.RData")


HO_spotted$year = factor(HO_spotted$year)  #get rid of years w/o data
test.spotted.year <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + pressure + 
                               sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                               sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                               age.sex:day + age.sex:day2 + age.sex:day3 + year:day + year:day2,
                             random.formula = Dry ~ speno,
                             data = HO_spotted,
                             EstMeth="REML",
                             timecol = "time", 
                             #ridge.reg = "global",
                             #lambda = 0.5,
                             group.vec = "AR1.ID")

save(test.spotted.year,file="test_spotted_year.RData")

#note year effects are going to be confounded with age-class since e.g. YOY are the only age class in 2005.
# HO_bearded$year = factor(HO_bearded$year)  
# test.bearded.year <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + temp2 + wind + pressure + 
#                                sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
#                                sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
#                                 year:day + year:day2,
#                              random.formula = Dry ~ speno,
#                              data = HO_bearded,
#                              EstMeth="REML",
#                              timecol = "time", 
#                              #ridge.reg = "global",
#                              #lambda = 0.5,
#                              group.vec = "AR1.ID")
# 
# save(test.bearded.year,file="test_bearded_year.RData")


#run models without weather or spatial effects for Okhotsk haul-out predictions

test.ribbon.nowx <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                          age.sex:day + age.sex:day2 + age.sex:day3,
                        random.formula = Dry ~ speno,
                        data = HO_ribbon,
                        EstMeth="REML",
                        timecol = "time", 
                        #ridge.reg = "global",
                        #lambda = 0.5,
                        group.vec = "AR1.ID")

save(test.ribbon.nowx,file="test_ribbon_nowx.RData")

test.spotted.nowx <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 +
                           sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                           sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 + 
                           age.sex:day + age.sex:day2 + age.sex:day3,
                         random.formula = Dry ~ speno,
                         data = HO_spotted,
                         EstMeth="REML",
                         timecol = "time", 
                         #ridge.reg = "global",
                         #lambda = 0.5,
                         group.vec = "AR1.ID")
save(test.spotted.nowx,file="test_spotted_nowx.RData")

test.bearded.nowx <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                           sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                           sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2+
                           Northing + Northing*day + Northing*day2,
                         random.formula = Dry ~ speno,
                         data = HO_bearded,
                         EstMeth="REML", #"ML",  #REML was default
                         timecol = "time", 
                         #ridge.reg = "global",
                         #lambda = 0.5,
                         group.vec = "AR1.ID")
save(test.bearded.nowx,file="test_bearded_nowx.RData")


#JUST use day-of-year, time-of-day
bearded.simple <- glmmLDTS(fixed.formula = Dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                                sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                                sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                              random.formula = Dry ~ speno,
                              data = HO_bearded,
                              EstMeth="REML", #"ML",  #REML was default
                              timecol = "time", 
                              #ridge.reg = "global",
                              #lambda = 0.5,
                              group.vec = "AR1.ID")

spotted.simple <- glmmLDTS(fixed.formula = Dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                        random.formula = Dry ~ speno,
                        data = HO_spotted,
                        EstMeth="REML", #"ML",  #REML was default
                        timecol = "time", 
                        #ridge.reg = "global",
                        #lambda = 0.5,
                        group.vec = "AR1.ID")

ribbon.simple <- glmmLDTS(fixed.formula = Dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                        random.formula = Dry ~ speno,
                        data = HO_ribbon,
                        EstMeth="REML", #"ML",  #REML was default
                        timecol = "time", 
                        #ridge.reg = "global",
                        #lambda = 0.5,
                        group.vec = "AR1.ID")


load('c:/users/paul.conn/git/BOSS/JayPowerR/glmmLDTS.bearded.fit.RDa')
max(glmmLDTS.bearded.fit$fit.table$mu)

load('c:/users/paul.conn/git/BOSS/JayPowerR/glmmLDTS.spotted.fit.RDa')
max(glmmLDTS.spotted.fit$fit.table$mu)

load('c:/users/paul.conn/git/BOSS/JayPowerR/glmmLDTS.ribbon.fit.RDa')
max(glmmLDTS.ribbon.fit$fit.table$mu)



