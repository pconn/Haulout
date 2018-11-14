# plot effects from species-specific haul-out model
library(mgcv)  #for predictions of Northing variable over time
expit<-function(x)1/(1+exp(-x))


###############################
#  Ribbon seals
###############################

load("test_ribbon.RData")  #model output list object - includes data set as one of its objects
data = test.ribbon$dataset  
FE = test.ribbon$fixed.effects
#FE = test.mod$coefficients

#possible plots: 
# wind effect by species
# species specific bivariate hour * day effects
# year * day of year effects - also plot total ice?
# species-specific Northing effect
# age.sex * day of year (x-axis = day of year effect)
# Northing by day-of-year

######################################
# Plot 1
# HO probability for each species by day of year and time of day
######################################
day.start = 60
day.end = 181
Day = (c(day.start:day.end)-120)/10
Day2 = Day^2
Day3 = Day^3
Hr = c(0:24)
Sin1 = sin(pi*Hr/12)
Cos1 = cos(pi*Hr/12)
Sin2 = sin(pi*Hr/6)
Cos2 = cos(pi*Hr/6)
Sin3 = sin(pi*Hr/4)
Cos3 = cos(pi*Hr/4)
#mean.yr.eff = mean(FE[FE[,"effect"]=="year","estimate"])
#mean.day.yr.eff = mean(FE[FE[,"effect"]=="day:year","estimate"])
#mean.day2.yr.eff = mean(FE[FE[,"effect"]=="day2:year","estimate"])
#mean.day3.yr.eff = mean(FE[FE[,"effect"]=="day3:year","estimate"])
#mean.temp = mean(data$temp2)
#mean.wind = mean(data$wind)
mean.Northing = mean(data$Northing)
Tmp.data = data[1:length(Day),]
Tmp.data$day = Day
Tmp.data$day2 = Day2
Tmp.data$day3 = Day3
Tmp.data2 = data[1:(length(Day)*24),]
Tmp.data2$day = rep(Day,24)
Tmp.data2$day2 = (rep(Day2,24))
Tmp.data2$day3 = (rep(Day3,24))
Tmp.data2$hour = rep(c(0:23),each=length(Day))
Tmp.data2$day.index = rep(c(1:length(Day)),24)

library(mgcv)
#Crap=data[1:24,]
#Crap$hour=c(0:23)
#gam.wind.hr = gam(wind~s(as.numeric(hour)),data=data)
#Wind.hr.pred = predict(gam.wind.hr,newdata=Crap)
#plot(Wind.hr.pred)  #not much of an effect of wind by hr of day

gam.baro <- gam(pressure~s(day),data=data)
gam.northing <- gam(Northing~s(day),data=data)
gam.temp <- gam(temp2~s(day)+s(as.numeric(hour)),data=data)
gam.wind <- gam(wind~s(day),data=data)
gam.precip <- gam(precip~s(day),data=data)
Temp.pred=predict(gam.temp,newdata=Tmp.data2)
Temp.mat = matrix(0,length(Day),24)
for(i in 1:nrow(Tmp.data2))Temp.mat[Tmp.data2[i,"day.index"],Tmp.data2[i,"hour"]+1]=Temp.pred[i]
Northing.pred = predict(gam.northing,newdata=Tmp.data)
Wind.pred = predict(gam.wind,newdata=Tmp.data)
Baro.pred = predict(gam.baro,newdata=Tmp.data)
Precip.pred = predict(gam.precip,newdata=Tmp.data)
#mean.pressure = mean(data$pressure)

n.days = day.end-day.start+1
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
#ASinter.vec = c("SUBADULT","ADULT.F","ADULT.M","ADULT.F")  #YOY get adult F interaction w day of year
#const.eff = FE[FE[,"effect"]=="intercept","estimate"]+mean.yr.eff+FE[FE[,"effect"]=="temp2","estimate"]*mean.temp+
#  FE[FE[,"effect"]=="wind","estimate"]*mean.wind+FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing+
#  FE[FE[,"effect"]=="pressure","estimate"]*mean.pressure+FE[FE[,"effect"]=="temp2:wind","estimate"]*mean.temp*mean.wind
const.eff = FE[FE[,"effect"]=="intercept","estimate"] #+FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing #FE[FE[,"effect"]=="temp2","estimate"]*mean.temp+
  #FE[FE[,"effect"]=="wind","estimate"]*mean.wind+ #FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing+
  #FE[FE[,"effect"]=="pressure","estimate"]*mean.pressure+FE[FE[,"effect"]=="temp2:wind","estimate"]*mean.temp*mean.wind
DHS = array(const.eff,dim=c(n.days,24,4))

for(iday in 1:n.days){
  DHS[iday,,]=DHS[iday,,]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    #FE[FE[,"effect"]=="Northing","estimate"]*Northing.pred[iday]+
    FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iday]+
    FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iday]+
    FE[FE[,"effect"]=="precip","estimate"]*Precip.pred[iday]
    #FE[FE[,"effect"]=="day:Northing","estimate"]*mean.Northing*Day[iday] #Northing.pred[iday]*Day[iday]+     
    #FE[FE[,"effect"]=="day2:Northing","estimate"]*mean.Northing*Day2[iday] #Northing.pred[iday]*Day2[iday]
}

for(ihr in 1:24){
    DHS[,ihr,] = DHS[,ihr,] + 
      FE[FE[,"effect"]=="sin1","estimate"]*Sin1[ihr] +
      FE[FE[,"effect"]=="cos1","estimate"]*Cos1[ihr] +
      FE[FE[,"effect"]=="sin2","estimate"]*Sin2[ihr] +
      FE[FE[,"effect"]=="cos2","estimate"]*Cos2[ihr] +
      FE[FE[,"effect"]=="sin3","estimate"]*Sin3[ihr] +
      FE[FE[,"effect"]=="cos3","estimate"]*Cos3[ihr]
}

for(iday in 1:n.days){
  for(ihr in 1:24){
    DHS[iday,ihr,] = DHS[iday,ihr,] + 
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin1:day3","estimate"]*Sin1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin2:day3","estimate"]*Sin2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin3:day3","estimate"]*Sin3[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos1:day3","estimate"]*Cos1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos2:day3","estimate"]*Cos2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos3:day3","estimate"]*Cos3[ihr]*Day3[iday]+     
      #FE[FE[,"effect"]=="day:temp2","estimate"]*mean.temp*Day[iday]+     
      #FE[FE[,"effect"]=="day2:temp2","estimate"]*mean.temp*Day2[iday]+     
      FE[FE[,"effect"]=="temp2","estimate"]*Temp.mat[iday,ihr]+   
      FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]
  }
}
   
for(iage in 1:4){
  DHS[,,iage]=DHS[,,iage]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

for(iday in 1:n.days){
  for(iage in 1:4){
    DHS[iday,,iage]=DHS[iday,,iage]+
        FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
        FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
        FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  }
}

library(reshape2)

Plot.df=melt(DHS,varnames=c("Day","Hour","Age.sex"))

#Plot.df = melt(DHS,varnames=c("Day","Hour","Species","Age.sex"))
Plot.df[,"value"]=1/(1+exp(-Plot.df[,"value"]))
#Plot.df[which(Plot.df[,"Species"]==1),"Species"]="Bearded"
#Plot.df[which(Plot.df[,"Species"]=='2'),"Species"]="Ribbon"
#Plot.df[which(Plot.df[,"Species"]=='3'),"Species"]="Spotted"
Plot.df[which(Plot.df[,"Age.sex"]==1),"Age.sex"]="Subadult"
Plot.df[which(Plot.df[,"Age.sex"]=='2'),"Age.sex"]="Adult.female"
Plot.df[which(Plot.df[,"Age.sex"]=='3'),"Age.sex"]="Adult.male"
Plot.df[which(Plot.df[,"Age.sex"]=='4'),"Age.sex"]="YOY"

#put an NA in place of modeled value for days of year w/ no observations before that day
for(iday in 1:(day.end-day.start+1)){
  if(sum((data$jday-59)==iday & data$age.sex=="SUB")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Subadult","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="ADULT.F")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Adult.female","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="ADULT.M")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Adult.male","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="YOY")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="YOY","value"]=NA
}

library(ggplot2)
library(RColorBrewer)
myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
DHS.plot = ggplot(Plot.df,aes(x=Day,y=Hour))+geom_raster(aes(fill=value))+facet_grid(~Age.sex)+scale_fill_gradientn(colours=myPalette(100),name='Pr(HO)')
DHS.plot

DH.ribbon = Plot.df

############################################## 
#plot 2: Haulout probability by wind
##############################################
Wind = c(0:20)
n=length(Wind)
Xtemp = model.matrix(test.ribbon$fixed.formula,data)
Wind.col = which(colnames(Xtemp)=="wind")
Wind.temp.col = which(colnames(Xtemp)=="temp2*wind")
Temp.col = which(colnames(Xtemp)=="temp2")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Wind.col]=Wind[i]/10
  Xtemp[,Wind.temp.col]=Xtemp[,Temp.col]*Wind[i]/10
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Wind=Wind,HOresp=HOresp)
wind.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Wind,y=HOresp))+xlab("Wind speed (m/s)")+ylab("Haul-out probability")
wind.plot = wind.plot + geom_density(data=data,fill="gray",aes(x=wind*10))
wind.plot

Wind.ribbon = Plot.df

##############################################
#  Plot 3: HO by temperature
##############################################
Temp = c(-10:10)
n=length(Temp)
Xtemp = model.matrix(test.ribbon$fixed.formula,data)
Wind.col = which(colnames(Xtemp)=="wind")
Wind.temp.col = which(colnames(Xtemp)=="temp2*wind")
Temp.col = which(colnames(Xtemp)=="temp2")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Temp.col]=Temp[i]/10
  Xtemp[,Wind.temp.col]=Xtemp[,Wind.col]*Temp[i]/10
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Temp=Temp*2.7-3.15,HOresp=HOresp)
temp.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Temp,y=HOresp))+xlab(expression(paste("Temperature (",~degree~C,")")))+ylab("Haul-out probability")
#temp.plot = temp.plot + geom_density(data=data,fill="gray",aes(x=temp*10))
temp.plot

Temp.ribbon = Plot.df

##############################################
#  Plot 4: HO by pressure
##############################################
Pressure = c(-50:47)
n=length(Pressure)
Xtemp = model.matrix(test.ribbon$fixed.formula,data)
P.col = which(colnames(Xtemp)=="pressure")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,P.col]=Pressure[i]/100
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Pressure=(Pressure*100+100000)/1000,HOresp=HOresp) # response in kPa
pressure.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Pressure,y=HOresp))+xlab("Pressure at mean sea level (kPa)")+ylab("Haul-out probability")+theme(text=element_text(size=14))
#temp.plot = temp.plot + geom_density(data=data,fill="gray",aes(x=temp*10))
pressure.plot

Pressure_ribbon = Plot.df

# pdf("ribbon_pressure.pdf")
#  pressure.plot
# dev.off()
# 
# jpeg("ribbon_pressure.jpg")
#  pressure.plot
# dev.off()

###############################################
#  Plot 5: Year effect models
##############################################
load("test_ribbon_year.RData")  #model output list object - includes data set as one of its objects
data = test.ribbon.year$dataset  
FE = test.ribbon.year$fixed.effects

Yrs = unique(test.ribbon.year$dataset$year)
n_yrs = length(Yrs)
n_days = day.end-day.start+1
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
iage = 2  #base predictions on adult females

const.eff = FE[FE[,"effect"]=="intercept","estimate"] + FE[FE[,"levels"]==AS.vec[iage],"estimate"]+
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1[12] +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1[12] +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2[12] +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2[12] +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3[12] +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3[12]

YrDay = array(const.eff,dim=c(n_yrs,n_days))

#load("mean_covs_for_ho_yr_effects.Rdata") #load space-averaged environmental covariates 
Pred.data = data.frame(day=Day,day2=Day2,day3=Day3,hour=rep(12,n_days))
Temp.pred = Wind.pred = Baro.pred = matrix(0,n_yrs,n_days)
gam.baro <- gam(pressure~s(day),data=data)
gam.temp <- gam(temp2~s(day)+s(as.numeric(hour)),data=data)
gam.wind <- gam(wind~s(day),data=data)
Temp.pred[1,]=predict(gam.temp,newdata=Pred.data)
Baro.pred[1,]=predict(gam.baro,newdata=Pred.data)
Wind.pred[1,]=predict(gam.wind,newdata=Pred.data)
for(iy in 2:n_yrs){
  Temp.pred[iy,]=Temp.pred[iy-1,]
  Baro.pred[iy,]=Baro.pred[iy-1,]
  Wind.pred[iy,]=Wind.pred[iy-1,]
  #Which.rows = which(Modeled_covs$year==as.numeric(as.character(Yrs[iy])))
  #Temp.pred[iy,]=c(Modeled_covs[Which.rows,"temp2"],Modeled_covs[Which.rows[length(Which.rows)],"temp2"])  #1 day short on covariates
  #Baro.pred[iy,]=c(Modeled_covs[Which.rows,"pressure"],Modeled_covs[Which.rows[length(Which.rows)],"pressure"])  #1 day short on covariates
  #Wind.pred[iy,]=c(Modeled_covs[Which.rows,"wind"],Modeled_covs[Which.rows[length(Which.rows)],"wind"])  #1 day short on covariates
}

ihr=12
for(iyr in 1:n_yrs){
  for(iday in 1:n.days){
  YrDay[iyr,iday]=YrDay[iyr,iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iyr,iday]+
    FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iyr,iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
    #FE[FE[,"effect"]=="sin1:day3","estimate"]*Sin1[ihr]*Day3[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
    #FE[FE[,"effect"]=="sin2:day3","estimate"]*Sin2[ihr]*Day3[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
    #FE[FE[,"effect"]=="sin3:day3","estimate"]*Sin3[ihr]*Day3[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
    #FE[FE[,"effect"]=="cos1:day3","estimate"]*Cos1[ihr]*Day3[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
    #FE[FE[,"effect"]=="cos2:day3","estimate"]*Cos2[ihr]*Day3[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
    #FE[FE[,"effect"]=="cos3:day3","estimate"]*Cos3[ihr]*Day3[iday]+     
    #FE[FE[,"effect"]=="day:temp2","estimate"]*mean.temp*Day[iday]+     
    #FE[FE[,"effect"]=="day2:temp2","estimate"]*mean.temp*Day2[iday]+     
    FE[FE[,"effect"]=="temp2","estimate"]*Temp.pred[iyr,iday]+   
    #FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]+
    FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0("ADULT.F",','),"estimate"]*Day[iday]+
    FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0("ADULT.F",','),"estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0("ADULT.F",','),"estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="day:year" & FE[,"levels"]==paste0(",",Yrs[iyr]),"estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2:year" & FE[,"levels"]==paste0(",",Yrs[iyr]),"estimate"]*Day2[iday]
  }
}

library(reshape2)

#Yrs = Yrs[-11]
Plot.df=data.frame("Year"=rep(Yrs,n_days),"Day"=rep(c(day.start:day.end),each=n_yrs),"HO"=plogis(as.vector(YrDay)))

#determine julian day maximums
HO_max_ribbon = rep(0,n_yrs)
for(iyr in 1:n_yrs){
  Cur_dat = Plot.df[which(Plot.df$Year==Yrs[iyr]),]
  HO_max_ribbon[iyr] = Cur_dat[which(Cur_dat$HO==max(Cur_dat$HO)),"Day"]
}
Sea_ice_Apr1 = c(438125,451250,558125,771875,736875,737500,487500,802500,643125,442500,480000,453125)
#par(mfrow(c(1,1)))
#plot(Sea_ice_Apr1[Yrs],HO_max_ribbon)
summary(lm(HO_max_ribbon~Sea_ice_Apr1[Yrs]))


#replace days w missing records with NA
for(iday in 1:(day.end-day.start+1)){
  for(iyr in 1:n_yrs){
    if(sum((data$jday-59)==iday & data$year==Yrs[iyr])==0)Plot.df[Plot.df$Day==(iday+59) & Plot.df$Year==Yrs[iyr],"HO"]=NA
  }
}

library(ggplot2)

YrDay.plot = ggplot(Plot.df)+geom_line(aes(x=Day,y=HO,colour=Year),size=1.3) + xlab("Julian day") + ylab("Haul-out probability")+ theme(text=element_text(size=14))
YrDay.plot

YrDay.ribbon = YrDay

pdf("HO_ribbon_YrDay.pdf")
YrDay.plot
dev.off()

jpeg("HO_ribbon_YrDay.jpg")
YrDay.plot
dev.off()

#take a look at 2011 

crap = crap=data[data$year=="2011",]
U.S = unique(crap$speno)
par(mfrow=c(5,1))
plot(crap[crap$speno==U.S[1],"jday"],crap[crap$speno==U.S[1],"percent_dry"],xlim=c(min(crap$jday),max(crap$jday)),xlab="Julian day",ylab="% Dry")
for(iind in 2:length(U.S)){
  plot(crap[crap$speno==U.S[iind],"jday"],crap[crap$speno==U.S[iind],"percent_dry"],xlim=c(min(crap$jday),max(crap$jday)),xlab="Julian day",ylab="% Dry")
}







# ###################### 
# #Northing by DOY plot    NOT DOING NORTHING FOR RIBBON AND SPOTTED
# ######################
# northing.sb=1.1  #southern Bering value
# northing.nb=0.91  #northern Bering value
# northing.ch=0.75  #chukchi value
# Northings = c(northing.sb,northing.nb,northing.ch)
# ihr = 12  #predictions all for solar noon
# iage=1 #plot for subadults
# 
# const.eff = FE[FE[,"effect"]=="intercept","estimate"]+
#   FE[FE[,"effect"]=="sin1","estimate"]*Sin1[ihr] +
#   FE[FE[,"effect"]=="cos1","estimate"]*Cos1[ihr] +
#   FE[FE[,"effect"]=="sin2","estimate"]*Sin2[ihr] +
#   FE[FE[,"effect"]=="cos2","estimate"]*Cos2[ihr] +
#   FE[FE[,"effect"]=="sin3","estimate"]*Sin3[ihr] +
#   FE[FE[,"effect"]=="cos3","estimate"]*Cos3[ihr] +
#   FE[FE[,"levels"]==AS.vec[iage],"estimate"]
# 
# 
# DN = matrix(const.eff,n.days,3)
# 
# for(iN in 1:3)DN[,iN]=DN[,iN]+FE[FE[,"effect"]=="Northing","estimate"]*Northings[iN]
# for(iday in 1:n.days){
#   DN[iday,]=DN[iday,]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
#     FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
#     FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
#     FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iday]+
#     FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iday]+
#     FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
#     FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
#     FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
#     FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
#     FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
#     FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
#     FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
#     FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
#     FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
#     FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
#     FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
#     FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
#     FE[FE[,"effect"]=="temp2","estimate"]*Temp.mat[iday,ihr]+
#     FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]+
#     FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
#     FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
#     FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
#   for(iN in 1:3){
#     DN[iday,iN]=DN[iday,iN]+FE[FE[,"effect"]=="day:Northing","estimate"]*Day[iday]*Northings[iN]+
#       FE[FE[,"effect"]=="day2:Northing","estimate"]*Day2[iday]*Northings[iN]
#   }
# }
# Plot.df=melt(DN,varnames=c("Day","Location"))
# 
# Plot.df[,"value"]=1/(1+exp(-Plot.df[,"value"]))
# Plot.df[which(Plot.df[,"Location"]==1),"Location"]="Southern.Bering"
# Plot.df[which(Plot.df[,"Location"]=='2'),"Location"]="Northern.Bering"
# Plot.df[which(Plot.df[,"Location"]=='3'),"Location"]="Chukchi"
# 
# #put an NA in place of modeled value for days of year w/ no observations before that day
# for(iday in 1:(day.end-day.start+1)){
#   if(sum((data$jday-59)==iday & data$Northing<0.85)==0)Plot.df[Plot.df$Day==iday & Plot.df$Location=="Chukchi","value"]=NA
#   if(sum((data$jday-59)==iday & (data$Northing<1 & data$Northing>0.85))==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Northern.Bering","value"]=NA
#   if(sum((data$jday-59)==iday & data$Northing>1)==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Southern.Bering","value"]=NA
# }
# 
# DN.plot = ggplot(Plot.df)+geom_line(aes(x=Day,y=value,colour=Location))
# DN.plot
# 
# DN.ribbon = Plot.df

##########################################
##############################################
#Spotted seals
##############################################
##########################################

load("test_spotted.RData")  #model output list object - includes data set as one of its objects
data = test.spotted$dataset  
FE = test.spotted$fixed.effects
#FE = test.mod$coefficients

gam.baro <- gam(pressure~s(day),data=data)
gam.northing <- gam(Northing~s(day),data=data)
gam.temp <- gam(temp2~s(day)+s(as.numeric(hour)),data=data)
gam.wind <- gam(wind~s(day),data=data)
gam.precip <- gam(precip~s(day),data=data)
Temp.pred=predict(gam.temp,newdata=Tmp.data2)
Temp.mat = matrix(0,length(Day),24)
for(i in 1:nrow(Tmp.data2))Temp.mat[Tmp.data2[i,"day.index"],Tmp.data2[i,"hour"]+1]=Temp.pred[i]
Northing.pred = predict(gam.northing,newdata=Tmp.data)
Wind.pred = predict(gam.wind,newdata=Tmp.data)
Baro.pred = predict(gam.baro,newdata=Tmp.data)
Precip.pred = predict(gam.precip,newdata=Tmp.data)

#mean.pressure = mean(data$pressure)

n.days = day.end-day.start+1
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
#ASinter.vec = c("SUBADULT","ADULT.F","ADULT.M","ADULT.F")  #YOY get adult F interaction w day of year
#const.eff = FE[FE[,"effect"]=="intercept","estimate"]+mean.yr.eff+FE[FE[,"effect"]=="temp2","estimate"]*mean.temp+
#  FE[FE[,"effect"]=="wind","estimate"]*mean.wind+FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing+
#  FE[FE[,"effect"]=="pressure","estimate"]*mean.pressure+FE[FE[,"effect"]=="temp2:wind","estimate"]*mean.temp*mean.wind
const.eff = FE[FE[,"effect"]=="intercept","estimate"] #+FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing #FE[FE[,"effect"]=="temp2","estimate"]*mean.temp+
#FE[FE[,"effect"]=="wind","estimate"]*mean.wind+ #FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing+
#FE[FE[,"effect"]=="pressure","estimate"]*mean.pressure+FE[FE[,"effect"]=="temp2:wind","estimate"]*mean.temp*mean.wind
DHS = array(const.eff,dim=c(n.days,24,4))

for(iday in 1:n.days){
  DHS[iday,,]=DHS[iday,,]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    #FE[FE[,"effect"]=="Northing","estimate"]*Northing.pred[iday]+
    FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iday]+
    FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iday]+
    FE[FE[,"effect"]=="precip","estimate"]*Precip.pred[iday]
  #FE[FE[,"effect"]=="day:Northing","estimate"]*mean.Northing*Day[iday] #Northing.pred[iday]*Day[iday]+     
  #FE[FE[,"effect"]=="day2:Northing","estimate"]*mean.Northing*Day2[iday] #Northing.pred[iday]*Day2[iday]
}

for(ihr in 1:24){
  DHS[,ihr,] = DHS[,ihr,] + 
    FE[FE[,"effect"]=="sin1","estimate"]*Sin1[ihr] +
    FE[FE[,"effect"]=="cos1","estimate"]*Cos1[ihr] +
    FE[FE[,"effect"]=="sin2","estimate"]*Sin2[ihr] +
    FE[FE[,"effect"]=="cos2","estimate"]*Cos2[ihr] +
    FE[FE[,"effect"]=="sin3","estimate"]*Sin3[ihr] +
    FE[FE[,"effect"]=="cos3","estimate"]*Cos3[ihr]
}

for(iday in 1:n.days){
  for(ihr in 1:24){
    DHS[iday,ihr,] = DHS[iday,ihr,] + 
      FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin1:day3","estimate"]*Sin1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin2:day3","estimate"]*Sin2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin3:day3","estimate"]*Sin3[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos1:day3","estimate"]*Cos1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos2:day3","estimate"]*Cos2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos3:day3","estimate"]*Cos3[ihr]*Day3[iday]+     
      #FE[FE[,"effect"]=="day:temp2","estimate"]*mean.temp*Day[iday]+     
      #FE[FE[,"effect"]=="day2:temp2","estimate"]*mean.temp*Day2[iday]+     
      FE[FE[,"effect"]=="temp2","estimate"]*Temp.mat[iday,ihr]+   
      FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]
  }
}

for(iage in 1:4){
  DHS[,,iage]=DHS[,,iage]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

for(iday in 1:n.days){
  for(iage in 1:4){
    DHS[iday,,iage]=DHS[iday,,iage]+
      FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  }
}

library(reshape2)

Plot.df=melt(DHS,varnames=c("Day","Hour","Age.sex"))

#Plot.df = melt(DHS,varnames=c("Day","Hour","Species","Age.sex"))
Plot.df[,"value"]=1/(1+exp(-Plot.df[,"value"]))
#Plot.df[which(Plot.df[,"Species"]==1),"Species"]="Bearded"
#Plot.df[which(Plot.df[,"Species"]=='2'),"Species"]="Ribbon"
#Plot.df[which(Plot.df[,"Species"]=='3'),"Species"]="Spotted"
Plot.df[which(Plot.df[,"Age.sex"]==1),"Age.sex"]="Subadult"
Plot.df[which(Plot.df[,"Age.sex"]=='2'),"Age.sex"]="Adult.female"
Plot.df[which(Plot.df[,"Age.sex"]=='3'),"Age.sex"]="Adult.male"
Plot.df[which(Plot.df[,"Age.sex"]=='4'),"Age.sex"]="YOY"


#put an NA in place of modeled value for days of year w/ no observations before that day
for(iday in 1:(day.end-day.start+1)){
  if(sum((data$jday-59)==iday & data$age.sex=="SUB")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Subadult","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="ADULT.F")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Adult.female","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="ADULT.M")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Adult.male","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="YOY")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="YOY","value"]=NA
}

library(ggplot2)
library(RColorBrewer)
myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
DHS.plot = ggplot(Plot.df,aes(x=Day,y=Hour))+geom_raster(aes(fill=value))+facet_grid(~Age.sex)+scale_fill_gradientn(colours=myPalette(100),name='Pr(HO)')
DHS.plot

DH.spotted = Plot.df

############################################## 
#plot 2: Haulout probability by wind
##############################################
Wind = c(0:20)
n=length(Wind)
Xtemp = model.matrix(test.spotted$fixed.formula,data)
Wind.col = which(colnames(Xtemp)=="wind")
Wind.temp.col = which(colnames(Xtemp)=="temp2*wind")
Temp.col = which(colnames(Xtemp)=="temp2")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Wind.col]=Wind[i]/10
  Xtemp[,Wind.temp.col]=Xtemp[,Temp.col]*Wind[i]/10
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Wind=Wind,HOresp=HOresp)
wind.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Wind,y=HOresp))+xlab("Wind speed (m/s)")+ylab("Haul-out probability")
wind.plot = wind.plot + geom_density(data=data,fill="gray",aes(x=wind*10))
wind.plot

Wind.spotted = Plot.df

##############################################
#  Plot 3: HO by temperature
##############################################
Temp = c(-10:10)
n=length(Temp)
Xtemp = model.matrix(test.spotted$fixed.formula,data)
Wind.col = which(colnames(Xtemp)=="wind")
Wind.temp.col = which(colnames(Xtemp)=="temp2*wind")
Temp.col = which(colnames(Xtemp)=="temp2")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Temp.col]=Temp[i]/10
  Xtemp[,Wind.temp.col]=Xtemp[,Wind.col]*Temp[i]/10
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Temp=Temp*2.7-3.15,HOresp=HOresp)
temp.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Temp,y=HOresp))+xlab(expression(paste("Temperature (",~degree~C,")")))+ylab("Haul-out probability")
#temp.plot = temp.plot + geom_density(data=data,fill="gray",aes(x=temp*10))
temp.plot

Temp.spotted = Plot.df

################ 
# precipitation plot
########
Precip = c(0:15)/10
n=length(Precip)
Xtemp = model.matrix(test.spotted$fixed.formula,data)
Precip.col = which(colnames(Xtemp)=="precip")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Precip.col]=Precip[i]
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Precip=Precip,HOresp=HOresp)
precip.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Precip,y=HOresp))+xlab(bquote('Convective precipitation (kg/'*m^2*')'))+ylab("Haul-out probability")
#temp.plot = temp.plot + geom_density(data=data,fill="gray",aes(x=temp*10))
precip.plot

jpeg("spotted_precip.jpg")
  precip.plot
dev.off()

pdf("spotted_precip.pdf")
  precip.plot
dev.off()

Precip.spotted = Plot.df


###############################################
#  Plot 4: Year effect models
##############################################
load("test_spotted_year.RData")  #model output list object - includes data set as one of its objects
data = test.spotted.year$dataset  
FE = test.spotted.year$fixed.effects

#Yrs = unique(test.spotted.year$dataset$year)  #out of order
Yrs = factor(c(as.character(2006:2011),"2014","2016","2017"))
n_yrs = length(Yrs)
n_days = day.end-day.start+1
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
iage = 2  #base predictions on adult females

const.eff = FE[FE[,"effect"]=="intercept","estimate"] + FE[FE[,"levels"]==AS.vec[iage],"estimate"]+
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1[12] +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1[12] +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2[12] +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2[12] +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3[12] +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3[12]

YrDay = array(const.eff,dim=c(n_yrs,n_days))

#spatially averaged environmental covariates 
#temporarily substitute w gams
Pred.data = data.frame(day=Day,day2=Day2,day3=Day3,hour=rep(12,n_days))
Temp.pred = Wind.pred = Baro.pred = matrix(0,n_yrs,n_days)
gam.baro <- gam(pressure~s(day),data=data)
gam.temp <- gam(temp2~s(day)+s(as.numeric(hour)),data=data)
gam.wind <- gam(wind~s(day),data=data)
Temp.pred[1,]=predict(gam.temp,newdata=Pred.data)
Baro.pred[1,]=predict(gam.baro,newdata=Pred.data)
Wind.pred[1,]=predict(gam.wind,newdata=Pred.data)
for(iy in 2:n_yrs){
  Temp.pred[iy,]=Temp.pred[iy-1,]
  Baro.pred[iy,]=Baro.pred[iy-1,]
  Wind.pred[iy,]=Wind.pred[iy-1,]
  #Which.rows = which(Modeled_covs$year==as.numeric(as.character(Yrs[iy])))
  #Temp.pred[iy,]=c(Modeled_covs[Which.rows,"temp2"],Modeled_covs[Which.rows[length(Which.rows)],"temp2"])  #1 day short on covariates
  #Baro.pred[iy,]=c(Modeled_covs[Which.rows,"pressure"],Modeled_covs[Which.rows[length(Which.rows)],"pressure"])  #1 day short on covariates
  #Wind.pred[iy,]=c(Modeled_covs[Which.rows,"wind"],Modeled_covs[Which.rows[length(Which.rows)],"wind"])  #1 day short on covariates
}


ihr=12
for(iyr in 1:n_yrs){
  for(iday in 1:n.days){
    YrDay[iyr,iday]=YrDay[iyr,iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
      FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
      FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iyr,iday]+
      FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iyr,iday]+
      FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin1:day3","estimate"]*Sin1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin2:day3","estimate"]*Sin2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin3:day3","estimate"]*Sin3[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos1:day3","estimate"]*Cos1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos2:day3","estimate"]*Cos2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos3:day3","estimate"]*Cos3[ihr]*Day3[iday]+     
      #FE[FE[,"effect"]=="day:temp2","estimate"]*mean.temp*Day[iday]+     
      #FE[FE[,"effect"]=="day2:temp2","estimate"]*mean.temp*Day2[iday]+     
      FE[FE[,"effect"]=="temp2","estimate"]*Temp.pred[iyr,iday]+   
      #FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]+
      FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0("ADULT.F",','),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0("ADULT.F",','),"estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0("ADULT.F",','),"estimate"]*Day3[iday]+
      FE[FE[,"effect"]=="day:year" & FE[,"levels"]==paste0(",",Yrs[iyr]),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="day2:year" & FE[,"levels"]==paste0(",",Yrs[iyr]),"estimate"]*Day2[iday]
  }
}

library(reshape2)

Plot.df=data.frame("Year"=rep(Yrs,n_days),"Day"=rep(c(day.start:day.end),each=n_yrs),"HO"=plogis(as.vector(YrDay)))

for(iday in 1:(day.end-day.start+1)){
  for(iyr in 1:n_yrs){
    if(sum((data$jday-59)==iday & data$year==Yrs[iyr])==0)Plot.df[Plot.df$Day==(iday+59) & Plot.df$Year==Yrs[iyr],"HO"]=NA
  }
}
HO_max_spotted = rep(0,n_yrs)
for(iyr in 1:n_yrs){
  Cur_dat = Plot.df[which(Plot.df$Year==Yrs[iyr]),]
  HO_max_spotted[iyr] = Cur_dat[which(Cur_dat$HO==max(Cur_dat$HO,na.rm=TRUE)),"Day"]
}
#par(mfrow(c(1,1)))
#plot(Sea_ice_Apr1[Yrs],HO_max_spotted)
summary(lm(HO_max_spotted~Sea_ice_Apr1[Yrs]))


YrDay.plot = ggplot(Plot.df)+geom_line(aes(x=Day,y=HO,colour=Year),size=1.3) + xlab("Julian day") + ylab("Haul-out probability")+ theme(text=element_text(size=14))
YrDay.plot

YrDay.spotted = YrDay

pdf("HO_spotted_YrDay.pdf")
YrDay.plot
dev.off()

jpeg("HO_spotted_YrDay.jpg")
YrDay.plot
dev.off()

#take a look at 2014, 2016 
crap = crap=data[data$year=="2014",]
U.S = unique(crap$speno)
par(mfrow=c(length(U.S),1))
plot(crap[crap$speno==U.S[1],"jday"],crap[crap$speno==U.S[1],"percent_dry"],xlim=c(min(crap$jday),max(crap$jday)),xlab="Julian day",ylab="% Dry")
for(iind in 2:length(U.S)){
  plot(crap[crap$speno==U.S[iind],"jday"],crap[crap$speno==U.S[iind],"percent_dry"],xlim=c(min(crap$jday),max(crap$jday)),xlab="Julian day",ylab="% Dry")
}






##########################################
##############################################
#Bearded seals
##############################################
##########################################

load("test_bearded.RData")  #model output list object - includes data set as one of its objects
data = test.bearded$dataset  
FE = test.bearded$fixed.effects
#FE = test.mod$coefficients

gam.baro <- gam(pressure~s(day),data=data)
gam.northing <- gam(Northing~s(day),data=data)
gam.temp <- gam(temp2~s(day)+s(as.numeric(hour)),data=data)
gam.wind <- gam(wind~s(day),data=data)
gam.precip <- gam(wind~s(precip),data=data)
Temp.pred=predict(gam.temp,newdata=Tmp.data2)
Temp.mat = matrix(0,length(Day),24)
for(i in 1:nrow(Tmp.data2))Temp.mat[Tmp.data2[i,"day.index"],Tmp.data2[i,"hour"]+1]=Temp.pred[i]
Northing.pred = predict(gam.northing,newdata=Tmp.data)
Wind.pred = predict(gam.wind,newdata=Tmp.data)
Baro.pred = predict(gam.baro,newdata=Tmp.data)
Precip.pred = predict(gam.precip,newdata=Tmp.data)
mean.Northing = mean(data$Northing)

n.days = day.end-day.start+1
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
#ASinter.vec = c("SUBADULT","ADULT.F","ADULT.M","ADULT.F")  #YOY get adult F interaction w day of year
#const.eff = FE[FE[,"effect"]=="intercept","estimate"]+mean.yr.eff+FE[FE[,"effect"]=="temp2","estimate"]*mean.temp+
#  FE[FE[,"effect"]=="wind","estimate"]*mean.wind+FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing+
#  FE[FE[,"effect"]=="pressure","estimate"]*mean.pressure+FE[FE[,"effect"]=="temp2:wind","estimate"]*mean.temp*mean.wind
const.eff = FE[FE[,"effect"]=="intercept","estimate"]#+FE[FE[,"effect"]=="Northing","estimate"]*mean.Northing#+ #FE[FE[,"effect"]=="temp2","estimate"]*mean.temp+
  #FE[FE[,"effect"]=="wind","estimate"]*mean.wind+
  #FE[FE[,"effect"]=="pressure","estimate"]*mean.pressure+FE[FE[,"effect"]=="temp2:wind","estimate"]*mean.temp*mean.wind
DHS = array(const.eff,dim=c(n.days,24,4))

for(iday in 1:n.days){
  DHS[iday,,]=DHS[iday,,]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="Northing","estimate"]*Northing.pred[iday]+
    FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iday]+
    FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iday]+
    FE[FE[,"effect"]=="precip","estimate"]*Precip.pred[iday]+
    FE[FE[,"effect"]=="day:Northing","estimate"]*Northing.pred[iday]*Day[iday]+     
    FE[FE[,"effect"]=="day2:Northing","estimate"]*Northing.pred[iday]*Day2[iday]
}

for(ihr in 1:24){
  DHS[,ihr,] = DHS[,ihr,] +
    FE[FE[,"effect"]=="sin1","estimate"]*Sin1[ihr] +
    FE[FE[,"effect"]=="cos1","estimate"]*Cos1[ihr] +
    FE[FE[,"effect"]=="sin2","estimate"]*Sin2[ihr] +
    FE[FE[,"effect"]=="cos2","estimate"]*Cos2[ihr] +
    FE[FE[,"effect"]=="sin3","estimate"]*Sin3[ihr] +
    FE[FE[,"effect"]=="cos3","estimate"]*Cos3[ihr]
}

for(iday in 1:n.days){
  for(ihr in 1:24){
    DHS[iday,ihr,] = DHS[iday,ihr,] +
      FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin1:day3","estimate"]*Sin1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin2:day3","estimate"]*Sin2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="sin3:day3","estimate"]*Sin3[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos1:day3","estimate"]*Cos1[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos2:day3","estimate"]*Cos2[ihr]*Day3[iday]+
      FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
      FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
      #FE[FE[,"effect"]=="cos3:day3","estimate"]*Cos3[ihr]*Day3[iday]+
      #FE[FE[,"effect"]=="day:temp2","estimate"]*mean.temp*Day[iday]+
      #FE[FE[,"effect"]=="day2:temp2","estimate"]*mean.temp*Day2[iday]+
      FE[FE[,"effect"]=="temp2","estimate"]*Temp.mat[iday,ihr]+
      FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]
  }
}

for(iage in 1:3){
  DHS[,,iage]=DHS[,,iage]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

# for(iday in 1:n.days){  #no age.sex * doy interactions
#   for(iage in 1:3){
#     DHS[iday,,iage]=DHS[iday,,iage]+
#       FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
#       FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
#       FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
#   }
# }

library(reshape2)

Plot.df=melt(DHS,varnames=c("Day","Hour","Age.sex"))

#Plot.df = melt(DHS,varnames=c("Day","Hour","Species","Age.sex"))
Plot.df[,"value"]=1/(1+exp(-Plot.df[,"value"]))
#Plot.df[which(Plot.df[,"Species"]==1),"Species"]="Bearded"
#Plot.df[which(Plot.df[,"Species"]=='2'),"Species"]="Ribbon"
#Plot.df[which(Plot.df[,"Species"]=='3'),"Species"]="Spotted"
Plot.df[which(Plot.df[,"Age.sex"]==1),"Age.sex"]="Subadult"
Plot.df[which(Plot.df[,"Age.sex"]=='2'),"Age.sex"]="Adult.female"
Plot.df[which(Plot.df[,"Age.sex"]=='3'),"Age.sex"]="Adult.male"
Plot.df[which(Plot.df[,"Age.sex"]=='4'),"Age.sex"]="YOY"

#put an NA in place of modeled value for days of year w/ no observations before that day
for(iday in 1:(day.end-day.start+1)){
  if(sum((data$jday-59)==iday & data$age.sex=="SUB")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Subadult","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="ADULT.F")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Adult.female","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="ADULT.M")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Adult.male","value"]=NA
  if(sum((data$jday-59)==iday & data$age.sex=="YOY")==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="YOY","value"]=NA
}

library(ggplot2)
library(RColorBrewer)
myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
DHS.plot = ggplot(Plot.df,aes(x=Day,y=Hour))+geom_raster(aes(fill=value))+facet_grid(~Age.sex)+scale_fill_gradientn(colours=myPalette(100),name='Pr(HO)')
DHS.plot

DH.bearded = Plot.df

############################################## 
#plot 2: Haulout probability by wind
##############################################
Wind = c(0:20)
n=length(Wind)
Xtemp = model.matrix(test.bearded$fixed.formula,data)
Wind.col = which(colnames(Xtemp)=="wind")
Wind.temp.col = which(colnames(Xtemp)=="temp2*wind")
Temp.col = which(colnames(Xtemp)=="temp2")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Wind.col]=Wind[i]/10
  Xtemp[,Wind.temp.col]=Xtemp[,Temp.col]*Wind[i]/10
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Wind=Wind,HOresp=HOresp)
wind.plot = ggplot()+geom_density(data=data,fill="gray",aes(x=wind*10))+xlab("Wind speed (m/s)")+ylab("Haul-out probability")
wind.plot = wind.plot + geom_line(size=1.3,data=Plot.df,aes(x=Wind,y=HOresp))
wind.plot

Wind.bearded = Plot.df

##############################################
#  Plot 3: HO by temperature
##############################################
Temp = c(-10:10)
n=length(Temp)
Xtemp = model.matrix(test.bearded$fixed.formula,data)
Wind.col = which(colnames(Xtemp)=="wind")
Wind.temp.col = which(colnames(Xtemp)=="temp2*wind")
Temp.col = which(colnames(Xtemp)=="temp2")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,Temp.col]=Temp[i]/10
  Xtemp[,Wind.temp.col]=Xtemp[,Wind.col]*Temp[i]/10
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Temp=Temp*2.7-3.15,HOresp=HOresp)
temp.plot = ggplot()+geom_line(size=1.3,data=Plot.df,aes(x=Temp,y=HOresp))+xlab(expression(paste("Temperature (",~degree~C,")")))+ylab("Haul-out probability")
#temp.plot = temp.plot + geom_density(data=data,fill="gray",aes(x=temp*10))
temp.plot

Temp.bearded = Plot.df

######
# pressure
######
Pressure = c(-50:47)
n=length(Pressure)
Xtemp = model.matrix(test.bearded$fixed.formula,data)
P.col = which(colnames(Xtemp)=="pressure")
HOresp = rep(0,n)
Ests = FE$estimate
Ests = Ests[-which(Ests==0)]
for(i in 1:n){
  Xtemp[,P.col]=Pressure[i]/100
  HOresp[i]=mean(expit(Xtemp%*%Ests))
}
Plot.df = data.frame(Pressure=(Pressure*100+100000)/1000,HOresp=HOresp) # response in kPa
pressure_bearded = Plot.df

###################### 
#Northing by DOY plot
######################
northing.sb=1.1  #southern Bering value
northing.nb=0.91  #northern Bering value
northing.ch=0.75  #chukchi value
Northings = c(northing.sb,northing.nb,northing.ch)
ihr = 12  #predictions all for solar noon
iage=1 #plot for subadults

const.eff = FE[FE[,"effect"]=="intercept","estimate"]+
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1[ihr] +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1[ihr] +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2[ihr] +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2[ihr] +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3[ihr] +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3[ihr] +
  FE[FE[,"levels"]==AS.vec[iage],"estimate"]


DN = matrix(const.eff,n.days,3)

for(iN in 1:3)DN[,iN]=DN[,iN]+FE[FE[,"effect"]=="Northing","estimate"]*Northings[iN]
for(iday in 1:n.days){
  DN[iday,]=DN[iday,]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="wind","estimate"]*Wind.pred[iday]+
    FE[FE[,"effect"]=="pressure","estimate"]*Baro.pred[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1[ihr]*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1[ihr]*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2[ihr]*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2[ihr]*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3[ihr]*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3[ihr]*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1[ihr]*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1[ihr]*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2[ihr]*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2[ihr]*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3[ihr]*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3[ihr]*Day2[iday]+
    FE[FE[,"effect"]=="temp2","estimate"]*Temp.mat[iday,ihr]+
    FE[FE[,"effect"]=="temp2:wind","estimate"]*Temp.mat[iday,ihr]*Wind.pred[iday]
    #FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
    #FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
    #FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  for(iN in 1:3){
    DN[iday,iN]=DN[iday,iN]+FE[FE[,"effect"]=="day:Northing","estimate"]*Day[iday]*Northings[iN]+
      FE[FE[,"effect"]=="day2:Northing","estimate"]*Day2[iday]*Northings[iN]
  }
}
Plot.df=melt(DN,varnames=c("Day","Location"))

Plot.df[,"value"]=1/(1+exp(-Plot.df[,"value"]))
Plot.df[which(Plot.df[,"Location"]==1),"Location"]="Southern.Bering"
Plot.df[which(Plot.df[,"Location"]=='2'),"Location"]="Northern.Bering"
Plot.df[which(Plot.df[,"Location"]=='3'),"Location"]="Chukchi"

#put an NA in place of modeled value for days of year w/ no observations before that day
for(iday in 1:(day.end-day.start+1)){
  if(sum((data$jday-59)==iday & data$Northing<0.85)==0)Plot.df[Plot.df$Day==iday & Plot.df$Location=="Chukchi","value"]=NA
  if(sum((data$jday-59)==iday & (data$Northing<1 & data$Northing>0.85))==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Northern.Bering","value"]=NA
  if(sum((data$jday-59)==iday & data$Northing>1)==0)Plot.df[Plot.df$Day==iday & Plot.df$Age.sex=="Southern.Bering","value"]=NA
}
Plot.df$Day = Plot.df$Day+59
DN.plot = ggplot(Plot.df)+geom_line(aes(x=Day,y=value,colour=Location),size=1.2)+ylab('Haul-out probability')+xlab('Julian day')+theme(text=element_text(size=14))
DN.plot
pdf("bearded_HO_by_location.pdf")
DN.plot
dev.off()
jpeg("bearded_HO_by_location.jpg")
DN.plot
dev.off()
DN.bearded = Plot.df


####### plot all together

DH.all = rbind(DH.bearded,DH.ribbon,DH.spotted)
DH.all$Species = c(rep("Bearded",nrow(DH.bearded)),rep("Ribbon",nrow(DH.ribbon)),rep("Spotted",nrow(DH.spotted)))

DH.all$Day = DH.all$Day+60  #convert to julian day
DH.plot = ggplot(DH.all,aes(x=Day,y=Hour))+geom_raster(aes(fill=value))+facet_grid(Species~Age.sex)+scale_fill_gradientn(colours=myPalette(100),name='Pr(HO)')+xlab('Julian day')+ylab('Solar hour')+theme(text=element_text(size=14))+scale_x_continuous(breaks=c(75,125,175))
DH.plot

pdf("HO_doy_tod_species_age.pdf")
DH.plot
dev.off()

jpeg("HO_doy_tod_species_age.jpeg")
DH.plot
dev.off()

Wind.all = rbind(Wind.bearded,Wind.ribbon,Wind.spotted)
Wind.all$Species=c(rep("Bearded",nrow(Wind.bearded)),rep("Ribbon",nrow(Wind.ribbon)),rep("Spotted",nrow(Wind.spotted)))
Dens.all = data.frame(Wind = c(test.bearded$dataset$wind,test.ribbon$dataset$wind,test.spotted$dataset$wind),Species=c(rep("Bearded",nrow(test.bearded$dataset)),rep("Ribbon",nrow(test.ribbon$dataset)),rep("Spotted",nrow(test.spotted$dataset))))

wind.plot = ggplot()+geom_density(data=Dens.all,fill="gray",aes(x=Wind*10))+geom_line(size=1.3,data=Wind.all,aes(x=Wind,y=HOresp))+xlab("Wind speed (m/s)")+ylab("Haul-out probability")+ facet_grid(~Species)+scale_x_continuous(limits=c(0,20))+theme(text=element_text(size=14))
wind.plot
pdf("HO_wind.pdf")
wind.plot
dev.off()

jpeg("HO_wind.jpeg")
wind.plot
dev.off()

Temp.all = rbind(Temp.bearded,Temp.ribbon,Temp.spotted)
Temp.all$Species=c(rep("Bearded",nrow(Temp.bearded)),rep("Ribbon",nrow(Temp.ribbon)),rep("Spotted",nrow(Temp.spotted)))
#Dens.all = data.frame(Temp = c(test.bearded$dataset$temp2,test.ribbon$dataset$temp2,test.spotted$dataset$temp2),Species=c(rep("Bearded",nrow(test.bearded$dataset)),rep("Ribbon",nrow(test.ribbon$dataset)),rep("Spotted",nrow(test.spotted$dataset))))

Temp.plot = ggplot()+geom_line(size=1.3,data=Temp.all,aes(x=Temp,y=HOresp,linetype=Species))+xlab(expression(paste("Temperature (",degree*C,")")))+ylab("Haul-out probability")+ theme(text=element_text(size=14),legend.key.width=unit(3,"line"))
Temp.plot
pdf("HO_temp.pdf")
Temp.plot
dev.off()

jpeg("HO_temp.jpeg")
Temp.plot
dev.off()

Pressure.all = rbind(pressure_bearded,Pressure_ribbon)
Pressure.all$Species=c(rep("Bearded",nrow(pressure_bearded)),rep("Ribbon",nrow(Pressure_ribbon)))
Pressure.plot = ggplot()+geom_line(size=1.3,data=Pressure.all,aes(x=Pressure,y=HOresp,linetype=Species))+xlab("Pressure at mean sea level (kPa)")+ylab("Haul-out probability")+theme(text=element_text(size=14),legend.key.width=unit(3,"line"))
Pressure.plot
pdf("HO_Pressure.pdf")
Pressure.plot
dev.off()

jpeg("HO_Pressure.jpeg")
Pressure.plot
dev.off()





