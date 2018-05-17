### compute stable stage distributions for haul-out paper

#load natural mortality estimates from meta-analysis
Ests <- read.csv('Surv_Haz_4_spp.csv')

Surv = matrix(0,4,40)
Surv[1,]=Ests[which(Ests[,"Species"]=="bearded" & Ests[,"Type"]=="Survival"),"Probability"]
Surv[2,]=Ests[which(Ests[,"Species"]=="ribbon" & Ests[,"Type"]=="Survival"),"Probability"]
Surv[3,]=Ests[which(Ests[,"Species"]=="ringed" & Ests[,"Type"]=="Survival"),"Probability"]
Surv[4,]=Ests[which(Ests[,"Species"]=="spotted" & Ests[,"Type"]=="Survival"),"Probability"]
Prob = Surv
for(iage in 2:40){
  Prob[,iage]=Surv[,iage]/Surv[,iage-1]
}
write.csv(t(Prob),"Surv_ests.csv")

Ages = c(0:39)
Surv_df = data.frame(Age = rep(Ages,4),Species=factor(c(rep("Bearded",40),rep("Ribbon",40),rep("Ringed",40),rep("Spotted",40))),Survival.prob=as.vector(t(Prob)))
library(ggplot2)
pdf("Survival.pdf")
ggplot(Surv_df)+geom_line(aes(colour=Species,x=Age,y=Survival.prob),size=1.5)+theme(text = element_text(size=16))+scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40))+scale_color_manual(values=c("red","blue","darkgreen","orange"))+ylab("Annual survival probability")
dev.off()

A = array(0,dim=c(4,40,40))  
for(isp in 1:4){
  for(iage in 1:39){
    A[isp,iage+1,iage]=Prob[isp,iage+1]  #assume post-breeding census so first survival is for yearlings; pup survival to be incorporated in fecundity term
  }
}



#ribbon fecundity from Fedoseev table 24
library(mgcv)
N_fec_ribbon = c(269,273,141+79+48,114+88+36,201,414+304+160)
y_fec_ribbon = round(c(0,0.013*83+0.183*82,0.589*141+0.62*79+0.563*48,0.886*114+0.977*88+0.944*36,
                       0.95*101+0.962*53+0.872*47,
                       0.936*78+26+0.974*310+0.974*304+0.963*160))
age = c(0:5)  #references ages 1 to 6+
Resp = cbind(y_fec_ribbon,N_fec_ribbon-y_fec_ribbon)
gam_ribbon = gam(Resp~s(age,k=3),family=binomial) 
fec_ribbon=0.5*predict(gam_ribbon,type="response")

#spotted seals
#age categories rported as 1, 2, ..., 8, 9+ ; using all data from Bering Sea from Table 38 in Fedoseev (2000)
N_fec_spotted = c(217,186,128,155,126,88,64,52,345)
y_fec_spotted = round(c(0,0,.05*21+0.06*17,0.18*22+0.12*.17+0.3*19+0.19*33+0.03*28,
                  0.59*22+0.18*11+0.65*37+0.44*18+0.47*17+0.09*21,
                  0.96*24+0.6*10+0.81*16+10+0.83*16+0.42*12,
                  8+0.63*9+0.96*22+6+14*.93+5,
                  4+8+17+5+0.92*13+5,
                  50+0.94*35+0.97*76+45+0.99*119+0.95*20))
age = c(0:8)
Resp = cbind(y_fec_spotted,N_fec_spotted-y_fec_spotted)
gam_spotted = gam(Resp~s(age,k=3),family=binomial) 
fec_spotted=0.5*predict(gam_spotted,type="response")

# ringed seals  from Table 11 of Fedoseev (2000)  - putting in a sample size of 1 for ages 0-3 since no observations in Bering
N_fec_ringed = c(1,1,1,14,11,8,14,11,33)
y_fec_ringed = round(c(0,0,0,0,0,2,8,.727*11,0.88*33))
Resp = cbind(y_fec_ringed,N_fec_ringed-y_fec_ringed)
gam_ringed = gam(Resp~s(age,k=3),family=binomial) 
fec_ringed=0.5*predict(gam_ringed,type="response")

# bearded seals - from table 47
age = c(0:8) #(indexes ages 1:9+)
N_fec_bearded = c(54,40,29,36,27,25,36,32,254)
y_fec_bearded = round(c(0,0,0,0,.191*23,.583*21+1,.667*19+.882*17,.947*18+.929*14,.935*107+.98*147))
Resp = cbind(y_fec_bearded,N_fec_bearded-y_fec_bearded)
gam_bearded = gam(Resp~s(age,k=3),family=binomial) 
fec_bearded=0.5*predict(gam_bearded,type="response")

Fec = matrix(0,9,4)
Fec[,1]=2*fec_bearded
Fec[1:6,2]=2*fec_ribbon
Fec[,3]=2*fec_ringed
Fec[,4]=2*fec_spotted
Fec[7:9,2]=Fec[6,2]
write.csv(Fec,file="Fec_table.csv")  #multiply by 2 so we can present parturition
#Ages = as.character(1:9)
#Ages[9]="9+"
Ages=c(1:9)
Fec_df = data.frame(Age = rep(Ages,4),Species=factor(c(rep("Bearded",9),rep("Ribbon",9),rep("Ringed",9),rep("Spotted",9))),Parturition=as.vector(Fec))
library(ggplot2)
pdf("Parturition.pdf")
ggplot(Fec_df)+geom_line(aes(colour=Species,x=Age,y=Parturition),size=1.5)+theme(text = element_text(size=16))+scale_x_continuous(breaks=c(1,3,5,7,9))+scale_color_manual(values=c("red","blue","darkgreen","orange"))
dev.off()


#fill in leslie matrices
A[1,1,2:10]=fec_bearded*Prob[1,1]
A[1,1,11:40]=A[1,1,10]
A[2,1,2:7]=fec_ribbon*Prob[2,1]
A[2,1,8:40]=A[2,1,7]
A[3,1,2:10]=fec_ringed*Prob[3,1]
A[3,1,11:40]=A[3,1,10]
A[4,1,2:10]=fec_spotted*Prob[4,1]
A[4,1,11:40]=A[4,1,10]

#Lambda values
cat("Lambda values: ")
eigen(A[1,,])$values[1]
eigen(A[2,,])$values[1]
eigen(A[3,,])$values[1]
eigen(A[4,,])$values[1]


#stable age distributions
Stable = matrix(0,4,40)
Stable[1,] = eigen(A[1,,])$vectors[,1]
Stable[1,] = Stable[1,]/sum(Stable[1,]) 
Stable[2,] = eigen(A[2,,])$vectors[,1]
Stable[2,] = Stable[2,]/sum(Stable[2,])
Stable[3,] = eigen(A[3,,])$vectors[,1]
Stable[3,] = Stable[3,]/sum(Stable[3,]) 
Stable[4,] = eigen(A[4,,])$vectors[,1]
Stable[4,] = Stable[4,]/sum(Stable[4,])

Stable= matrix(as.numeric(Stable),4,40)


#stable stage calculations
Stage = array(0,dim=c(4,2,3))  #species, sex (M/F), stage class

#ribbon
#males - based on proportions mature as reported by Tikhomirov (1966) via Fedoseev 2000, pg 126 
Stage[2,1,1]=Stable[2,1]
Stage[2,1,2]=Stable[2,2]*0.85 + Stable[2,3]*.16
Stage[2,1,3]=Stable[2,2]*.15 + Stable[2,3]*.84 + sum(Stable[2,4:40])

#females - based on percentage of immagure seals at age from Table 24, Fedoseev 2000
N_mat_ribbon = c(269,273,141+79+48,114+88+36,201,104,310+304+160)
y_mat_ribbon = c(0.109*46+0.062*64,.74*108+.88*83+.585*82,.979*141+.987*79+.938*48)
prop_mat_ribbon = rep(1,40)
prop_mat_ribbon[1]=0
prop_mat_ribbon[2:4]=y_mat_ribbon/N_mat_ribbon[1:3]
Stage[2,2,1]=Stable[2,1]
Stage[2,2,2]=sum((1-prop_mat_ribbon[2:40])*Stable[2,2:40])
Stage[2,2,3]=1-sum(Stage[2,2,1:2])

#ringed seals - assume male maturity similar to females (no data provided in Fedoseev 2000 for males)
#use table 11 for Bering Strait
prop_mat_ringed = rep(1,39)
prop_mat_ringed[1:5]=0
prop_mat_ringed[6:8]=c(0.55,0.75,0.84)
Stage[3,,1]=Stable[3,1]
Stage[3,1,2]=Stage[3,2,2]=sum((1-prop_mat_ringed)*Stable[3,2:40])
Stage[3,1,3]=Stage[3,2,3]=sum(prop_mat_ringed*Stable[3,2:40])

#spotted seals - use table 38 for females; use note that "males attain maturity at 5-6" (pg 159) to get male maturity 
prop_mat_spotted = rep(1,39)
prop_mat_spotted[1]=0
y_mat_spotted = c(0,.04*22+.05*20,.24*21+.43*14+.12*33+.12*17+.35*22+.1*21,.73*22+12+57*.84+.65*17+.7*19+.46*28,.96*22+.91*11+.97*37+18+.94*17+.71*21,60+.96*16+.92*12)
prop_mat_spotted[1:6]=y_mat_spotted/N_fec_spotted[1:6]
Stage[4,,1]=Stable[4,1]
Stage[4,2,2]=sum((1-prop_mat_spotted)*Stable[4,2:40])
Stage[4,2,3]=sum(prop_mat_spotted*Stable[4,2:40])

prop_mat_spotted_male = rep(1,39)
prop_mat_spotted_male[1:4]=0
prop_mat_spotted_male[5]=1/3
prop_mat_spotted_male[6]=2/3
Stage[4,1,2]=sum((1-prop_mat_spotted_male)*Stable[4,2:40])
Stage[4,1,3]=sum(prop_mat_spotted_male*Stable[4,2:40])

# bearded: use table 47 for females; use note in text (pg 180) that "males attain maturity at the age of 5-6 years"
#males
Stage[1,,1] = Stable[1,1]
Stage[1,1,2]=sum((1-prop_mat_spotted_male)*Stable[1,2:40])
Stage[1,1,3]=sum(prop_mat_spotted_male*Stable[1,2:40])
#females
prop_mat_bearded = rep(1,39)
y_mat_bearded = c(0,0.143*7,0.067*28+1,0.346*25+0.546*11,0.619*23+4,.833*21+.75*4)
prop_mat_bearded[1:6]=y_mat_bearded/N_fec_bearded[1:6]
Stage[1,2,2]=sum((1-prop_mat_bearded)*Stable[1,2:40])
Stage[1,2,3]=sum(prop_mat_bearded*Stable[1,2:40])

Stage = Stage*0.5 #assumes 50/50 sex ratio

#plot stage classes using bar chart 
Plot_df = data.frame(Proportion=as.numeric(Stage))
Plot_df$Species = rep(c("Bearded","Ribbon","Ringed","Spotted"),6)
Plot_df$Sex = rep(c(rep("Male",4),rep("Female",4)),3)
Plot_df$Stage = c(rep("YOY",8),rep("Sub-adult",8),rep("Adult",8))
Plot_df$Stage = factor(Plot_df$Stage,levels=c("YOY","Sub-adult","Adult")) #change so order of x-axis isn't alphabetical
library(ggplot2)
SS_plot = ggplot(data=Plot_df,aes(x=Stage,y=Proportion,fill=Sex))+
  geom_bar(stat="identity",position = 'dodge') + 
  geom_text(aes(label=sprintf("%0.2f", round(Proportion, digits = 2))), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_wrap(~Species)+
  scale_fill_manual(values=c("black","gray"))+
  ylim(0,0.45)

jpeg('stable_stage.jpg')
  SS_plot
dev.off()

pdf('stable_stage.pdf')
  SS_plot
dev.off()


####################
#  fit GLMPMs to spotted and ribbon seal data w/ (1) just DOY and TOD and (2) DOY, TOD, and age-sex class.
#  then produce predictions at solar noon and adjust with stable stage distributions to compare
####################

load("HO_df.RData")  #from process_haulout_data4.R
library(glmmLDTS)

HO_ribbon = HO_df[HO_df$species=="Hf",]
AS.ribbon <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                          sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                          sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                          age.sex:day + age.sex:day2 + age.sex:day3,
                        random.formula = Dry ~ speno,
                        data = HO_ribbon,
                        EstMeth="REML",
                        timecol = "time", 
                        group.vec = "AR1.ID")
Const.ribbon <- glmmLDTS(fixed.formula = Dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                        sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                        sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                      random.formula = Dry ~ speno,
                      data = HO_ribbon,
                      EstMeth="REML",
                      timecol = "time", 
                      group.vec = "AR1.ID")


HO_spotted = HO_df[HO_df$species=="Pl",]
AS.spotted <- glmmLDTS(fixed.formula = Dry ~ age.sex + sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                        sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                        sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2 +   
                        age.sex:day + age.sex:day2 + age.sex:day3,
                      random.formula = Dry ~ speno,
                      data = HO_spotted,
                      EstMeth="REML",
                      timecol = "time", 
                      group.vec = "AR1.ID")
Const.spotted <- glmmLDTS(fixed.formula = Dry ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + day + day2+ day3 + 
                           sin1*day + cos1*day + sin2*day + cos2*day + sin3*day + cos3*day +
                           sin1*day2 + cos1*day2 + sin2*day2 + cos2*day2 + sin3*day2 + cos3*day2,
                         random.formula = Dry ~ speno,
                         data = HO_spotted,
                         EstMeth="REML",
                         timecol = "time", 
                         group.vec = "AR1.ID")

## ribbon seal predictions at solar noon
#  1) constant
data = Const.ribbon$dataset  
FE = Const.ribbon$fixed.effects
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
day.start = day.start.ribbon = min(data$jday)
day.end = day.end.ribbon = max(data$jday)
n.days=day.end-day.start+1
n.days.ribbon=n.days
Day = (c(day.start:day.end)-120)/10
Day2 = Day^2
Day3 = Day^3
hr = 12
Sin1 = sin(pi*hr/12)
Cos1 = cos(pi*hr/12)
Sin2 = sin(pi*hr/6)
Cos2 = cos(pi*hr/6)
Sin3 = sin(pi*hr/4)
Cos3 = cos(pi*hr/4)


const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = rep(const.eff,n.days)

for(iday in 1:n.days){
  Pred[iday]=Pred[iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}
Pred.ribbon.const = plogis(Pred)
plot(Pred.ribbon.const)

# 1b) ribbon - age specific
FE = AS.ribbon$fixed.effects

const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = matrix(const.eff,4,n.days)

for(iday in 1:n.days){
  Pred[,iday]=Pred[,iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}

for(iage in 1:4){
  Pred[iage,]=Pred[iage,]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

for(iday in 1:n.days){
  for(iage in 1:4){
    Pred[iage,iday]=Pred[iage,iday]+
      FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  }
}

Pred = plogis(Pred)
Pred.ribbon.AS = Pred[1,]*sum(Stage[2,,2]) + Pred[2,]*Stage[2,2,3] + Pred[3,]*Stage[2,1,3] + Pred[4,]*sum(Stage[2,,1])

lines(Pred.ribbon.AS)

#### spotted seals

## spotted seal predictions at solar noon
#  1) constant
data = Const.spotted$dataset  
FE = Const.spotted$fixed.effects
AS.vec = c("SUB","ADULT.F","ADULT.M","YOY")
day.start = min(data$jday)
day.end = max(data$jday)
n.days=day.end-day.start+1
Day = (c(day.start:day.end)-120)/10
Day2 = Day^2
Day3 = Day^3
hr = 12
Sin1 = sin(pi*hr/12)
Cos1 = cos(pi*hr/12)
Sin2 = sin(pi*hr/6)
Cos2 = cos(pi*hr/6)
Sin3 = sin(pi*hr/4)
Cos3 = cos(pi*hr/4)


const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = rep(const.eff,n.days)

for(iday in 1:n.days){
  Pred[iday]=Pred[iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}
Pred.spotted.const = plogis(Pred)
plot(Pred.spotted.const)

# 1b) spotted - age specific
FE = AS.spotted$fixed.effects

const.eff = FE[FE[,"effect"]=="intercept","estimate"] +     
  FE[FE[,"effect"]=="sin1","estimate"]*Sin1 +
  FE[FE[,"effect"]=="cos1","estimate"]*Cos1 +
  FE[FE[,"effect"]=="sin2","estimate"]*Sin2 +
  FE[FE[,"effect"]=="cos2","estimate"]*Cos2 +
  FE[FE[,"effect"]=="sin3","estimate"]*Sin3 +
  FE[FE[,"effect"]=="cos3","estimate"]*Cos3
Pred = matrix(const.eff,4,n.days)

for(iday in 1:n.days){
  Pred[,iday]=Pred[,iday]+FE[FE[,"effect"]=="day","estimate"]*Day[iday]+
    FE[FE[,"effect"]=="day2","estimate"]*Day2[iday]+
    FE[FE[,"effect"]=="day3","estimate"]*Day3[iday]+
    FE[FE[,"effect"]=="sin1:day","estimate"]*Sin1*Day[iday]+
    FE[FE[,"effect"]=="sin1:day2","estimate"]*Sin1*Day2[iday]+
    FE[FE[,"effect"]=="sin2:day","estimate"]*Sin2*Day[iday]+
    FE[FE[,"effect"]=="sin2:day2","estimate"]*Sin2*Day2[iday]+
    FE[FE[,"effect"]=="sin3:day","estimate"]*Sin3*Day[iday]+
    FE[FE[,"effect"]=="sin3:day2","estimate"]*Sin3*Day2[iday]+
    FE[FE[,"effect"]=="cos1:day","estimate"]*Cos1*Day[iday]+
    FE[FE[,"effect"]=="cos1:day2","estimate"]*Cos1*Day2[iday]+
    FE[FE[,"effect"]=="cos2:day","estimate"]*Cos2*Day[iday]+
    FE[FE[,"effect"]=="cos2:day2","estimate"]*Cos2*Day2[iday]+
    FE[FE[,"effect"]=="cos3:day","estimate"]*Cos3*Day[iday]+
    FE[FE[,"effect"]=="cos3:day2","estimate"]*Cos3*Day2[iday]
}

for(iage in 1:4){
  Pred[iage,]=Pred[iage,]+FE[FE[,"levels"]==AS.vec[iage],"estimate"]
}

for(iday in 1:n.days){
  for(iage in 1:4){
    Pred[iage,iday]=Pred[iage,iday]+
      FE[FE[,"effect"]=="age.sex:day" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day[iday]+
      FE[FE[,"effect"]=="age.sex:day2" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day2[iday]+
      FE[FE[,"effect"]=="age.sex:day3" & FE[,"levels"]==paste0(AS.vec[iage],','),"estimate"]*Day3[iday]
  }
}

Pred = plogis(Pred)
Pred.spotted.AS = Pred[1,]*sum(Stage[4,,2]) + Pred[2,]*Stage[4,2,3] + Pred[3,]*Stage[4,1,3] + Pred[4,]*sum(Stage[4,,1])

lines(Pred.spotted.AS)

#plot AS-specific vs. constant AS mean haulout behavior
library(ggplot2)
Plot.df = data.frame(Day=c(rep(c(day.start.ribbon:day.end.ribbon),2),rep(c(day.start:day.end),2)),
                     HO=c(Pred.ribbon.const,Pred.ribbon.AS,Pred.spotted.const,Pred.spotted.AS),
                     Species = c(rep("Ribbon",n.days.ribbon*2),rep("Spotted",n.days*2)),
                     Type = c(rep("Constant",n.days.ribbon),rep("Stable stage",n.days.ribbon),rep("Constant",n.days),rep("Stable stage",n.days)))
SS.plot = ggplot(Plot.df) + geom_line(size=1.2,aes(x=Day,y=HO,colour=Species,linetype=Type)) + xlab("Julian day")+ylab("Haul-out probability")+ theme(text=element_text(size=14))

jpeg("HO_stable_stage.jpg")
  SS.plot
dev.off()

pdf("HO_stable_stage.pdf")
  SS.plot
dev.off()

