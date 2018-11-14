#check to see if there are still any missing data from VerHoef et al.'s haul-out
# analysis circa BOSS power calculations (2012?) and the data I'm no analyzing

load('c:/users/paul.conn/git/BOSS/JayPowerR/glmmLDTS.bearded.fit.RDa')
Speno_Jay = unique(glmmLDTS.bearded.fit$dataset$Speno)

load('c:/users/paul.conn/git/haulout/test_bearded.RData')
Speno_new = unique(test.bearded$dataset$speno)

Which.missing = which((Speno_Jay %in% Speno_new)==FALSE)

#see if any of these data in March - June
Data.missing = glmmLDTS.bearded.fit$dataset[which(glmmLDTS.bearded.fit$dataset$Speno %in% Speno_Jay[Which.missing]),]
Data.missing2 = Data.missing[which(Data.missing$month %in% c(3:6)),]

#okay, not very much missing - may have just been removed as out of deployment data or something

#compare old fitted values for May to new fitted values for May
mean(glmmLDTS.bearded.fit$dataset$Drytime[glmmLDTS.bearded.fit$dataset$month==6])
#0.66
mean(test.bearded$dataset$percent_dry[test.bearded$dataset$month==6])
#33
mean(test.bearded$dataset$Dry[test.bearded$dataset$month==6])
#29

#limit new data to speno's which were analyzed in common
All_speno = unique(c(Speno_Jay,Speno_new))
Common_speno = All_speno[which(All_speno %in% Speno_Jay & All_speno %in% Speno_new)]

mean(glmmLDTS.bearded.fit$dataset$Drytime[glmmLDTS.bearded.fit$dataset$month==6 & glmmLDTS.bearded.fit$dataset$Speno %in% Common_speno])
mean(test.bearded$dataset$Dry[test.bearded$dataset$month==6 & test.bearded$dataset$speno %in% Common_speno])
mean(test.bearded$dataset$percent_dry[test.bearded$dataset$month==6 & test.bearded$dataset$speno %in% Common_speno])

glmmLDTS.bearded.fit$dataset$year = as.numeric(as.character(glmmLDTS.bearded.fit$dataset$year))
test.bearded$dataset$year = as.numeric(as.character(test.bearded$dataset$year))
test.bearded$dataset$hour = as.numeric(as.character(test.bearded$dataset$hour))

Jay_common1 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[1],]
Me_common1 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[1],]
Jay_common1 = Jay_common1[Jay_common1$month>=3,]
plot(Jay_common1[1:1192,"Drytime"],type="l")
lines(Me_common1$percent_dry/100,col="blue")
#off by some hours - prob Jays is in UTC, mine is in Solar
# but raw data is almost the same ~21 percent hauled out

Jay_common2 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[2],]
Me_common2 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[2],]
Jay_common2 = Jay_common2[Jay_common2$month>=3,]
plot(Jay_common2[,"Drytime"],type="l")
lines(Me_common2$percent_dry/100,col="blue")
mean(Jay_common2$Drytime)
mean(Me_common2$percent_dry)
# 24%

Jay_common3 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[3],]
Me_common3 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[3],]
Jay_common3 = Jay_common3[Jay_common3$month>=3,]
plot(Jay_common3[,"Drytime"],type="l")
lines(Me_common3$percent_dry/100,col="blue")
mean(Jay_common3$Drytime)
mean(Me_common3$percent_dry)
# both about 21%

Jay_common4 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[4],]
Me_common4 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[4],]
Jay_common4 = Jay_common4[Jay_common4$month>=3,]
plot(Jay_common4[,"Drytime"],type="l")
lines(Me_common4$percent_dry/100,col="blue")
mean(Jay_common4$Drytime)
mean(Me_common4$percent_dry)
# both about 28%

Jay_common6 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[6],]
Me_common6 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[6],]
Jay_common6 = Jay_common6[Jay_common6$month>=3,]
Me_day = 366*(Me_common6$year-min(Me_common6$year))+Me_common6$jday+Me_common6$hour/24
Jay_day = 366*(Jay_common6$year-min(Me_common6$year))+Jay_common6$yrday+Jay_common6$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common6$Drytime,type="l")
plot(Me_day,Me_common6$percent_dry/100,col="blue",type="l")
mean(Jay_common6$Drytime)
mean(Me_common6$percent_dry)
# Jay - 26%, me- 19%, but more data for Jay

Jay_common6 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[6],]
Me_common6 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[6],]
Jay_common6 = Jay_common6[Jay_common6$month>=3,]
Me_day = 365*(Me_common6$year-min(Me_common6$year))+Me_common6$jday+Me_common6$hour/24
Jay_day = 365*(Jay_common6$year-min(Me_common6$year))+Jay_common6$yrday+Jay_common6$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common6$Drytime,type="l")
plot(Me_day,Me_common6$percent_dry/100,col="blue",type="l")
mean(Jay_common6$Drytime)
mean(Me_common6$percent_dry)
# both about 19.5%

Jay_common7 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[7],]
Me_common7 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[7],]
Jay_common7 = Jay_common7[Jay_common7$month>=3,]
Me_day = 365*(Me_common7$year-min(Me_common7$year))+Me_common7$jday+Me_common7$hour/24
Jay_day = 365*(Jay_common7$year-min(Me_common7$year))+Jay_common7$yrday+Jay_common7$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common7$Drytime,type="l")
plot(Me_day,Me_common7$percent_dry/100,col="blue",type="l")
mean(Jay_common7$Drytime)
mean(Me_common7$percent_dry)
# both about 20.4%

Jay_common8 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[8],]
Me_common8 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[8],]
Jay_common8 = Jay_common8[Jay_common8$month>=3,]
Me_day = 365*(Me_common8$year-min(Me_common8$year))+Me_common8$jday+Me_common8$hour/24
Jay_day = 365*(Jay_common8$year-min(Me_common8$year))+Jay_common8$yrday+Jay_common8$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common8$Drytime,type="l")
plot(Me_day,Me_common8$percent_dry/100,col="blue",type="l")
mean(Jay_common8$Drytime)
mean(Me_common8$percent_dry)
# 22.0 vs 21.6

Jay_common9 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[9],]
Me_common9 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[9],]
Jay_common9 = Jay_common9[Jay_common9$month>=3,]
Me_day = 365*(Me_common9$year-min(Me_common9$year))+Me_common9$jday+Me_common9$hour/24
Jay_day = 365*(Jay_common9$year-min(Me_common9$year))+Jay_common9$yrday+Jay_common9$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common9$Drytime,type="l")
plot(Me_day,Me_common9$percent_dry/100,col="blue",type="l")
mean(Jay_common9$Drytime)
mean(Me_common9$percent_dry)
# 21.8 vs 21.6

Jay_common10 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[10],]
Me_common10 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[10],]
Jay_common10 = Jay_common10[Jay_common10$month>=3,]
Me_day = 365*(Me_common10$year-min(Me_common10$year))+Me_common10$jday+Me_common10$hour/24
Jay_day = 365*(Jay_common10$year-min(Me_common10$year))+Jay_common10$yrday+Jay_common10$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common10$Drytime,type="p")
plot(Me_day,Me_common10$percent_dry/100,col="blue",type="p")
mean(Jay_common10$Drytime)
mean(Me_common10$percent_dry)
#Jay 43, me 7  - looks like Jay's includes a period of 100% haulout at end of time series

Jay_common11 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[11],]
Me_common11 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[11],]
Jay_common11 = Jay_common11[Jay_common11$month>=3,]
Me_day = 365*(Me_common11$year-min(Me_common11$year))+Me_common11$jday+Me_common11$hour/24
Jay_day = 365*(Jay_common11$year-min(Me_common11$year))+Jay_common11$yrday+Jay_common11$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common11$Drytime,type="p")
plot(Me_day,Me_common11$percent_dry/100,col="blue",type="p")
mean(Jay_common11$Drytime)
mean(Me_common11$percent_dry)
# me - 9.8%, Jay - 42%
#my dataset only includes about 8 days of observations.  Jays includes far more
#EB2009_3001
#come of these in the second year (2010) are a long string of 1's

Jay_common12 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[12],]
Me_common12 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[12],]
Jay_common12 = Jay_common12[Jay_common12$month>=3,]
Me_day = 365*(Me_common12$year-min(Me_common12$year))+Me_common12$jday+Me_common12$hour/24
Jay_day = 365*(Jay_common12$year-min(Me_common12$year))+Jay_common12$yrday+Jay_common12$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common12$Drytime,type="p")
plot(Me_day,Me_common12$percent_dry/100,col="blue",type="p")
mean(Jay_common12$Drytime)
mean(Me_common12$percent_dry)
# jay - 49 to me - 10%.  Again, more data and long strings of 1 in Jay's that are presumably erroneous

Jay_common13 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[13],]
Me_common13 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[13],]
Jay_common13 = Jay_common13[Jay_common13$month>=3,]
Me_day = 365*(Me_common13$year-min(Me_common13$year))+Me_common13$jday+Me_common13$hour/24
Jay_day = 365*(Jay_common13$year-min(Me_common13$year))+Jay_common13$yrday+Jay_common13$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common13$Drytime,type="p")
plot(Me_day,Me_common13$percent_dry/100,col="blue",type="p")
mean(Jay_common13$Drytime)
mean(Me_common13$percent_dry)
# jay - 11 to me - 8%.  More data for Jay (prob July)

Jay_common14 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[14],]
Me_common14 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[14],]
Jay_common14 = Jay_common14[Jay_common14$month>=3,]
Me_day = 365*(Me_common14$year-min(Me_common14$year))+Me_common14$jday+Me_common14$hour/24
Jay_day = 365*(Jay_common14$year-min(Me_common14$year))+Jay_common14$yrday+Jay_common14$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common14$Drytime,type="p")
plot(Me_day,Me_common14$percent_dry/100,col="blue",type="p")
mean(Jay_common14$Drytime)
mean(Me_common14$percent_dry)
# jay - 43 to me - 9.6%.  More data for Jay (prob going into July) as well as a subsequent year w a big string of 1's

Jay_common15 = glmmLDTS.bearded.fit$dataset[glmmLDTS.bearded.fit$dataset$Speno==Common_speno[15],]
Me_common15 = test.bearded$dataset[test.bearded$dataset$speno==Common_speno[15],]
Jay_common15 = Jay_common15[Jay_common15$month>=3,]
Me_day = 365*(Me_common15$year-min(Me_common15$year))+Me_common15$jday+Me_common15$hour/24
Jay_day = 365*(Jay_common15$year-min(Me_common15$year))+Jay_common15$yrday+Jay_common15$hour/24
par(mfrow=c(2,1))
plot(Jay_day,Jay_common15$Drytime,type="p")
plot(Me_day,Me_common15$percent_dry/100,col="blue",type="p")
mean(Jay_common15$Drytime)
mean(Me_common15$percent_dry)
# jay - 7 to me - 9%.  More data for Jay (prob going into July) 

