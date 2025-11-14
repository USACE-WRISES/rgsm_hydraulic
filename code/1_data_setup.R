# Data set-up for models
# Uses Yackulic code to format data for model runs
# Code before OOS data is unaltered from Yackulic et al. 2022
# Code for OOS data is new because of compiling additional years of data

# Calculation of the larval carrying capacity index is in subsequent script

####################
#### STEP 1: Read in all data
#####################



t2asir<-read.csv("data/yackulic2022_data/ASIRQuery4Combined.csv")
tdry<-read.csv("data/yackulic2022_data/rivereyes_v2.csv")
angQ<-read.csv("data/yackulic2022_data/abq_gage_08330000.csv")
sanaQ<-read.csv("data/yackulic2022_data/SanAcacia_gage_08354900.csv")
t2vie<-read.csv("data/yackulic2022_data/releases 8_26_2019.csv")
converter<-read.csv("data/yackulic2022_data/fws2br_rmconverter.csv")
tmeso<-read.csv("data/yackulic2022_data/mesohab_sum.csv")
trescue<-read.csv("data/yackulic2022_data/Fish Rescue 2009-2020.csv")
phiR<-read.csv("data/yackulic2022_data/rescue_surv.csv")
ee1<-read.csv("data/yackulic2022_data/sum_ee1.csv")
ee2<-read.csv("data/yackulic2022_data/sum_ee2.csv")
ee3<-read.csv("data/yackulic2022_data/sum_ee3.csv")
ee4<-read.csv("data/yackulic2022_data/sum_ee4.csv")
ee5<-read.csv("data/yackulic2022_data/sum_ee5.csv")
aQ<-read.csv("data/yackulic2022_data/aQ.csv",header=T)[,c(3:19)]
sQ<-read.csv("data/yackulic2022_data/sQ.csv",header=T)[,c(3:19)]



######################
#### STEP 2: Define variables and reformat data
######################



# define two grids
grid<-c(55.4,116,170,210.1) # start and end points for 3 main river segments
vgrid<-seq(55.4,210.2,0.2) # represents all 200 meter segments in study area
crossgrid<-findInterval(vgrid[-length(vgrid)]+.05,grid)
gridst_end<-rbind(range(which(crossgrid==1)),range(which(crossgrid==2)),range(which(crossgrid==3)))
years<-c(2002:2018) ##years used to fit model
Nyears<-length(years)
# reformat drying data
tdry$ustrata<-findInterval(tdry$URM,vgrid)
tdry$dstrata<-findInterval(tdry$LRM,vgrid)
tdry$jul<-NA
for (i in 1:length(tdry$jul)){
  tdry$jul[i]<-julian(as.Date(paste(tdry$Year[i],tdry$month[i],tdry$dom[i],sep="-")),as.Date(paste(tdry$Year[i],"03","31",sep="-")))[[1]]}
fdryday<-matrix(NA,nrow=(Nyears),ncol=(length(vgrid)-1))
ldryday<-matrix(NA,nrow=(Nyears),ncol=(length(vgrid)-1))
for (k in 1:(length(vgrid)-1)){
  for (j in 1:Nyears){
    temp2<-subset(tdry$jul,k<=tdry$ustrata&k>=tdry$dstrata&tdry$Year==(2001+j))
    fdryday[j,k]<-ifelse(length(temp2)==0,NA,min(temp2))
    ldryday[j,k]<-ifelse(length(temp2)==0,NA,max(temp2))
  }}
# reformat augmented fish release data
vkeep<-match(c("Date","Month","Year","C","L","S","Number","RM_up","RM_down"),names(t2vie))
tvie<-t2vie[,vkeep]
tvie$strata_up<-findInterval(converter[match(tvie$RM_up,converter[,1]),2],vgrid)
tvie$strata_down<-findInterval(converter[match(tvie$RM_down,converter[,1]),2],vgrid)
tlen<-length(tvie$strata_down)
vie<-tvie
for (i in 1:tlen){
  temp<-tvie$strata_up[i]-tvie$strata_down[i]
  vie$Number[i]<-tvie$Number[i]/(1+temp)
  if (temp>0){
    for (j in 1:temp){
      temp2<-tvie[i,]
      temp2$strata_down<-temp2$strata_down+j
      temp2$Number<-tvie$Number[i]/(1+temp)
      vie<-rbind(vie,temp2)
    }}}
tvie<-vie
tvie$C<-ifelse(tvie$Year>2007|(tvie$Year==2007&tvie$Month>9),paste(tvie$C),"")
tvie$L<-ifelse(tvie$Year>2007|(tvie$Year==2007&tvie$Month>9),paste(tvie$L),"")
tvie$S<-ifelse(tvie$Year>2007|(tvie$Year==2007&tvie$Month>9),paste(tvie$S),"")
tvie$VIE<-paste(tvie$C,tvie$L,sep="")
Vcodes<-sort(unique(tvie$VIE))
tvie$Vc<-match(tvie$VIE,Vcodes)
tvie$period<-ifelse(tvie$Month<4,8*(tvie$Year-2002)+1,ifelse(tvie$Month>9,8*(tvie$Year-2002)+9,8*(tvie$Year-2002)+tvie$Month-3))
tvie$winmon<-ifelse(tvie$Month<4,4-tvie$Month,ifelse(tvie$Month>9,16-tvie$Month,0))
vkeep<-match(c("Number","strata_down","Vc","period","winmon"),names(tvie))
vie<-tvie[,vkeep]
vie<-subset(vie,vie$period!=130&vie$period!=132) ###complicate things too much and not that many fish
rm(tvie,vkeep)
####
Nstrata<-3
Nvstrata<-length(vgrid)-1
nw_um<-matrix(0,nrow=30,ncol=Nstrata) #summarized at coarse grid
w_um<-matrix(0,nrow=6,ncol=Nstrata)
nmons<-matrix(0,nrow=6,ncol=Nstrata)
t1<-subset(vie,vie$Vc==1&8*(vie$period/8-floor(vie$period/8))==1)
t2<-subset(vie,vie$Vc==1&8*(vie$period/8-floor(vie$period/8))!=1)
t3<-subset(vie,vie$Vc>1)
tw_m<-array(0,dim=c((Nyears-5),Nvstrata,2)) #summarized at fine grid
w_m<-array(0,dim=c((Nyears-5),Nstrata,2)) #summarized at coarse grid
for (i in 1:6){
  for (j in 1:Nstrata){
    temp<-subset(t1,((t1$period-1)/8+1)==i&crossgrid[t1$strata_down]==j)
    if (dim(temp)[1]>0){
      w_um[i,j]<-sum(temp[,1])
      nmons[i,j]<-round(sum(temp[,1]*temp[,5])/sum(temp[,1]))
    }}}
for (i in 1:6){
  for (t in 1:5){
    for (j in 1:Nstrata){
      temp<-subset(t2,(floor(t2$period/8)+1)==i&(8*(t2$period/8-floor(t2$period/8))-1)==t&crossgrid[t2$strata_down]==j)
      if (dim(temp)[1]>0){
        nw_um[((i-1)*5+t),j]<-sum(temp[,1])
      }}}}
for (i in 1:(Nyears-5)){
  for (j in 1:Nvstrata){
    temp<-subset(t3,((t3$period-1)/8-5)==i&t3$strata_down==j)
    if (dim(temp)[1]>0){
      tw_m[i,j,1]<-sum(subset(temp[,1],temp$winmon<3))
      tw_m[i,j,2]<-sum(subset(temp[,1],temp$winmon>3))
    }}
  for (j in 1:Nstrata){
    w_m[i,j,1]<-sum(tw_m[i,gridst_end[j,1]:gridst_end[j,2],1])
    w_m[i,j,2]<-sum(tw_m[i,gridst_end[j,1]:gridst_end[j,2],2])
  }}
w_um[6,3]<-sum(nw_um[25:26,3])
##### done reformating vie release data	

#simplify and reformat monitoring database
keep<-match(c("ProjectName","DateSampled","SiteID","year", "month","RMStart","HabitatNumber","Habitat","SamplingEffort","RepeatedSamplingNumber",
              "Species","DepletionNumber","Gear","NumberCaptured","AgeClass","LengthSL","LengthMinSL","LengthMaxSL","VIEColor","VIELocation"),names(t2asir))
tasir<-t2asir[,keep]
tasir<-subset(tasir,tasir$Habitat!=""&tasir$Gear!="larval"&tasir$year>2001)
ghab<-c("","BW","MCPLPO","MCPO","MCSHPLPO","MCSHPO","PO","SCPLPO","SCPO","SCSHPLPO","SCSHPO","SHPO","MCED","MCSHED","SCED","SCSHED") # sampled pool habitats
ghabAv<-c("BW NS","PO NS","SHPO NS") #not sampled pool habitats - only measured during population estimation
bhabAv<-c("RU NS","SHRU NS") #not sampled run/riffle habitat.
tasir$cHab<-ifelse(is.na(match(tasir$Habitat,ghab))==FALSE,1,ifelse(is.na(match(tasir$Habitat,ghabAv))==FALSE,2,
                                                                    ifelse(is.na(match(tasir$Habitat,bhabAv))==FALSE,4,3)))
tasir$strata<-findInterval(tasir$RMStart,vgrid)
tasir$period<-8*(tasir$year-2002)+tasir$month-3
tasir<-subset(tasir,tasir$period>0)
extractdom<-function(x){
  lx<-nchar(x)
  tx<-numeric()
  for (i in 1:lx){
    tx[i]<-substr(x,i,i)}
  wx<-which(tx=="/")
  as.numeric(substr(x,wx[1]+1,wx[2]-1))}
tasir$dom<-NA
tasir$jul<-NA
tasir$dry<-NA
for (j in 1:length(tasir[,1])){
  tasir$dom[j]<-extractdom(paste(tasir$DateSampled[j]))
  tasir$jul[j]<-julian(as.Date(paste(tasir$year[j],tasir$month[j],tasir$dom[j],sep="-")),as.Date(paste(tasir$year[j],"03","31",sep="-")))[[1]]
  tf1<-fdryday[(tasir$year[j]-2001),tasir$strata[j]]
  tl1<-ldryday[(tasir$year[j]-2001),tasir$strata[j]]
  tasir$dry[j]<-ifelse(is.na(tf1)==TRUE|(tf1>tasir$jul[j])|(tl1<tasir$jul[j]),0,1)}
# subset different asir data for more reformating
###start with regular monitoring data
t2mon<-subset(tasir,tasir$ProjectName=="Hybognathus Amarus Population Monitoring"&tasir$month>3&tasir$month<11&tasir$cHab!=4&tasir$SamplingEffort!=""&tasir$SamplingEffort!="#N/A"&tasir$dry==0)
t2mon$uni_id<-paste(t2mon$year,t2mon$jul,t2mon$RMStart,t2mon$HabitatNumber)
t2mon$C<-ifelse(t2mon$year>2007,substr(paste(t2mon$VIEColor),1,1),"")
t2mon$L<-ifelse(t2mon$year>2007,substr(paste(t2mon$VIELocation),1,1),"")
t2mon$VC<-match(paste(t2mon$C,t2mon$L,sep=""),Vcodes)
t2mon$VC[which(t2mon$AgeClass==0&t2mon$VC>1)]<-1
# summarize data to seine haul
tmon<-data.frame(uni=sort(unique(t2mon$uni_id)))
tmon$RMStart<-t2mon$RMStart[match(tmon[,1],t2mon$uni_id)]
tmon$period<-t2mon$period[match(tmon[,1],t2mon$uni_id)]
tmon$jul<-t2mon$jul[match(tmon[,1],t2mon$uni_id)]
tmon$year<-t2mon$year[match(tmon[,1],t2mon$uni_id)]
tmon$month<-t2mon$month[match(tmon[,1],t2mon$uni_id)]
tmon$date<-t2mon$DateSampled[match(tmon[,1],t2mon$uni_id)]
tmon$angQ<-angQ$cfs[match(tmon$date,angQ$Date)]
tmon$sanaQ<-sanaQ$cfs[match(tmon$date,sanaQ$Date)]
tmon$HaulNo<-t2mon$HabitatNumber[match(tmon[,1],t2mon$uni_id)]
tmon$cHab<-t2mon$cHab[match(tmon[,1],t2mon$uni_id)]
tmon$effort<-as.numeric(paste(t2mon$SamplingEffort[match(tmon[,1],t2mon$uni_id)]))
tmon$strata<-t2mon$strata[match(tmon[,1],t2mon$uni_id)]
tmon$cQ<-ifelse(tmon$strata>573,tmon$angQ,tmon$sanaQ)
tmon$hybama0<-0
tmon$hybama1<-0
tmon$hybama2<-0
tmon$hybama3<-0
tmon$hybama4<-0
tmon$hybama5<-0
tmon$hybama6<-0
tmon$hybama7<-0
tmon$hybama8<-0
tmon$hybama9<-0
tmon$hybama10<-0
tmon$hybama11<-0
tmon$hybama12<-0
tmon$hybama13<-0
tmon$hybama14<-0
tmon$hybama01<-0
tstart<-which(names(tmon)=="hybama0")
t3mon<-subset(t2mon,t2mon$Species=="HYBAMA")
for (i in 1:length(t3mon[,1])){
  tx<-match(t3mon$uni_id[i],tmon$uni)
  ty<-ifelse(is.na(t3mon$AgeClass[i])==TRUE,1+tstart+length(Vcodes),ifelse(
    t3mon$AgeClass[i]==0,tstart,t3mon$VC[i]+tstart))
  tz<-t3mon$NumberCaptured[i]
  tmon[tx,ty]<-tmon[tx,ty]+tz
}
tmon$type<-ifelse(tmon$cHab==1,1,2)
tmon$SPt_id<-paste(tmon$year,tmon$jul,tmon$strata,tmon$type)
tmon<-subset(tmon,tmon$cQ<1000) #remove seine hauls were discharge was greater than 1000 cfs
####subset data to summarize all hauls on the same day, in the same habitat and same river segment
mon<-data.frame(spt=sort(unique(tmon$SPt_id)))
mon$period<-tmon$period[match(mon[,1],tmon$SPt_id)]
mon$jul<-tmon$jul[match(mon[,1],tmon$SPt_id)]
mon$year<-tmon$year[match(mon[,1],tmon$SPt_id)]
mon$month<-tmon$month[match(mon[,1],tmon$SPt_id)]
mon$strata<-tmon$strata[match(mon[,1],tmon$SPt_id)]
mon$type<-tmon$type[match(mon[,1],tmon$SPt_id)]
mon$effort<-0
mon$cQ<-0
mon$hybama0<-0
mon$hybama1<-0
mon$hybama2<-0
mon$hybama3<-0
mon$hybama4<-0
mon$hybama5<-0
mon$hybama6<-0
mon$hybama7<-0
mon$hybama8<-0
mon$hybama9<-0
mon$hybama10<-0
mon$hybama11<-0
mon$hybama12<-0
mon$hybama13<-0
mon$hybama14<-0
mon$hybama01<-0
mean_names<-c("cQ")
sum_names<-c("effort","hybama0","hybama1","hybama2","hybama3","hybama4","hybama5",
             "hybama6","hybama7","hybama8","hybama9","hybama10","hybama11","hybama12",
             "hybama13","hybama14","hybama01")
tsm<-match(sum_names,names(tmon))
sm<-match(sum_names,names(mon))
tmn<-match(mean_names,names(tmon))
smn<-match(mean_names,names(mon))
for (i in 1:length(mon[,1])){
  temp<-subset(tmon[,tsm],tmon$SPt_id==mon$spt[i])
  mon[i,sm]<-colSums(temp)
  temp<-subset(tmon[,tmn],tmon$SPt_id==mon$spt[i])
  mon[i,smn]<-mean(temp)
}
vie_names<-c("hybama2","hybama3","hybama4","hybama5","hybama6","hybama7","hybama8",
             "hybama9","hybama10","hybama11","hybama12","hybama13","hybama14")
mon$hybama_vie<-rowSums(mon[,vie_names])
mon$c_id<-paste(mon$year,mon$jul,crossgrid[mon$strata],mon$type)

## fit vie dispersal analysis described in appendix S1 to determine weights 
wm<-tw_m[,,1]+tw_m[,,2]
wm2<-tapply(wm[1,],crossgrid,sum)
for (j in 2:12){
  wm2<-rbind(wm2,tapply(wm[j,],crossgrid,sum))}
tt<-subset(mon,mon$year>2007)
tt$temp<-NA
for (i in 1:length(tt$temp)){
  tt$temp[i]<-wm2[(tt$year[i]-2007),crossgrid[tt$strata[i]]]}
tt2<-subset(tt,tt$temp!=0&tt$jul<60)
fitcauchy<-function(par){ # in units of 200 m sites
  a<-exp(par[2:3])
  predV<-numeric()	
  for (i in 1:401){
    t2<-ifelse(crossgrid==crossgrid[tt2$strata[i]],1,0)
    t3<-wm[(tt2$year[i]-2007),]*t2
    t4<-which(t3!=0)
    t4b<-which(t2!=0)
    t5<-numeric()
    for (j in 1:length(t4)){t5[j]<-sum(dcauchy(t4b,t4[j]+par[5],exp(par[1])))}
    t6<-sum(t3[t4]*dcauchy(tt2$strata[i],t4+par[5],exp(par[1]))/t5)
    predV[i]<-t6*tt2$effort[i]*a[tt2$type[i]]/400000}
  -1*sum(dnbinom(tt2$hybama_vie,mu=predV,size=exp(par[4]),log=TRUE))}
m<-optim(c(3,2,0,-2,0),fitcauchy,method="BFGS",hessian=TRUE)
mweights<-matrix(NA,nrow=12,ncol=Nvstrata)
for (i in 1:12){
  for (s in 1:3){
    t3<-wm[i,]*ifelse(crossgrid==s,1,0)
    t4<-which(t3!=0)
    if (length(t4)>0){
      t5<-numeric()
      for (j in 1:length(t4)){t5[j]<-sum(dcauchy(c(gridst_end[s,1]:gridst_end[s,2]),t4[j]+m$par[5],exp(m$par[1])))}
      t6<-numeric()
      for (k in gridst_end[s,1]:gridst_end[s,2]){t6[(k+1-gridst_end[s,1])]<-sum(t3[t4]*dcauchy(k,t4+m$par[5],exp(m$par[1]))/t5)}
      mweights[i,gridst_end[s,1]:gridst_end[s,2]]<-t6/sum(t6)
    }}}
mon$vieW<-0
for (i in 1:length(mon[,1])){
  mon$vieW[i]<-ifelse(mon$year[i]<2008,NA,mweights[(mon$year[i]-2007),mon$strata[i]])}
mon$Cstrata<-crossgrid[mon$strata]
monC<-data.frame(cid=sort(unique(mon$c_id)))
monC$jul<-mon$jul[match(monC[,1],mon$c_id)]
monC$year<-mon$year[match(monC[,1],mon$c_id)]
monC$month<-mon$month[match(monC[,1],mon$c_id)]
monC$Cstrata<-crossgrid[mon$strata[match(monC[,1],mon$c_id)]]
monC$type<-mon$type[match(monC[,1],mon$c_id)]
monC$effort<-0
monC$cQ<-0
monC$hybama0<-0
monC$hybama1<-0
monC$hybama01<-0
monC$hybama_vie<-0
monC$vieW<-0
monC$period<-mon$period[match(monC[,1],mon$c_id)]

mean_names<-c("cQ")
sum_names<-c("effort","hybama0","hybama1","hybama_vie","hybama01","vieW")
tsm<-match(sum_names,names(mon))
sm<-match(sum_names,names(monC))
tmn<-match(mean_names,names(mon))
smn<-match(mean_names,names(monC))
for (i in 1:length(monC[,1])){
  temp<-subset(mon[,tsm],mon$c_id==monC$cid[i])
  monC[i,sm]<-colSums(temp)
  temp<-subset(mon[,tmn],mon$c_id==monC$cid[i])
  monC[i,smn]<-mean(temp)
}
#### Done reformating monthly April to October monitoring data
#### reformat November data
t2monre<-subset(tasir,tasir$ProjectName=="Hybognathus Amarus Population Monitoring Repeated"&tasir$VIEColor==""&tasir$VIELocation==""&tasir$dry==0)
t2monre$uni_id<-paste(t2monre$year,t2monre$jul,t2monre$RMStart,t2monre$HabitatNumber)
tmonre<-data.frame(uni_YSH=sort(unique(t2monre$uni_id)))
tmonre$Date<-t2monre$Date[match(tmonre[,1],t2monre$uni_id)]
tmonre$jul<-t2monre$jul[match(tmonre[,1],t2monre$uni_id)]
tmonre$angQ<-angQ$cfs[match(tmonre$Date,angQ$Date)]
tmonre$sanaQ<-sanaQ$cfs[match(tmonre$Date,sanaQ$Date)]
tmonre$RMStart<-t2monre$RMStart[match(tmonre[,1],t2monre$uni_id)]
tmonre$year<-t2monre$year[match(tmonre[,1],t2monre$uni_id)]
tmonre$HabNo<-t2monre$HabitatNumber[match(tmonre[,1],t2monre$uni_id)]
tmonre$cHab<-t2monre$cHab[match(tmonre[,1],t2monre$uni_id)]
tmonre$effort<-as.numeric(paste(t2monre$SamplingEffort[match(tmonre[,1],t2monre$uni_id)]))
tmonre$Cstrata<-crossgrid[t2monre$strata[match(tmonre[,1],t2monre$uni_id)]]
tmonre$cQ<-ifelse(tmonre$Cstrata==3,tmonre$angQ,tmonre$sanaQ)
tmonre$hybama_1<-0
tmonre$hybama_2<-0
tmonre$hybama_3<-0
tmonre$hybama_4<-0

tstart<-which(names(tmonre)=="hybama_1")
t3monre<-subset(t2monre,t2monre$Species=="HYBAMA")
for (i in 1:length(t3monre[,1])){
  tx<-match(t3monre$uni_id[i],tmonre$uni_YSH)
  ty<-t3monre$RepeatedSamplingNumber[i]+tstart-1
  tz<-t3monre$NumberCaptured[i]
  tmonre[tx,ty]<-tmonre[tx,ty]+tz
}
tmonre$type<-ifelse(tmonre$cHab==1,1,2)
tmonre$period<-t2monre$period[match(tmonre[,1],t2monre$uni_id)]
tmonre$SPt_id<-paste(tmonre$year,tmonre$jul,tmonre$Cstrata,tmonre$type)
## summarize by day, river segement and habitat type
monre<-data.frame(spt=sort(unique(tmonre$SPt_id)))
monre$period<-tmonre$period[match(monre[,1],tmonre$SPt_id)]
monre$year<-tmonre$year[match(monre[,1],tmonre$SPt_id)]
monre$month<-tmonre$month[match(monre[,1],tmonre$SPt_id)]
monre$Cstrata<-tmonre$Cstrata[match(monre[,1],tmonre$SPt_id)]
monre$type<-tmonre$type[match(monre[,1],tmonre$SPt_id)]
monre$jul<-tmonre$jul[match(monre[,1],tmonre$SPt_id)]
monre$effort<-0
monre$cQ<-0
monre$hybama_1<-0
mean_names<-c("cQ")
sum_names<-c("effort","hybama_1")
tsm<-match(sum_names,names(tmonre))
sm<-match(sum_names,names(monre))
tmn<-match(mean_names,names(tmonre))
smn<-match(mean_names,names(monre))
for (i in 1:length(monre[,1])){
  temp<-subset(tmonre[,tsm],tmonre$SPt_id==monre$spt[i])
  monre[i,sm]<-colSums(temp)
  temp<-subset(tmonre[,tmn],tmonre$SPt_id==monre$spt[i])
  monre[i,smn]<-mean(temp)
}
monre<-subset(monre,monre$cQ<1000)
#### only looking at these data for habitat availabilityestimates
t2est<-subset(tasir,tasir$ProjectName=="Hybognathus Amarus Population Estimation")
t2est$uni_YSH<-paste(t2est$year,t2est$RMStart,t2est$HabitatNumber)
test<-data.frame(uni_YSH=unique(t2est$uni_YSH))
test$RMStart<-t2est$RMStart[match(test[,1],t2est$uni_YSH)]
test$year<-t2est$year[match(test[,1],t2est$uni_YSH)]
test$Date<-t2est$Date[match(test[,1],t2est$uni_YSH)]
test$angQ<-angQ$cfs[match(test$Date,angQ$Date)]
test$sanaQ<-sanaQ$cfs[match(test$Date,sanaQ$Date)]
test$HabNo<-t2est$HabitatNumber[match(test[,1],t2est$uni_YSH)]
test$cHab<-t2est$cHab[match(test[,1],t2est$uni_YSH)]
test$effort<-as.numeric(paste(t2est$SamplingEffort[match(test[,1],t2est$uni_YSH)]))
test$Cstrata<-crossgrid[t2est$strata[match(test[,1],t2est$uni_YSH)]]
test$cQ<-ifelse(test$Cstrata==3,test$angQ,test$sanaQ)
#
test$uni_YS<-paste(test$year,test$RMStart)
av<-data.frame(uni_YS=unique(test$uni_YS))
av$RMStart<-test$RMStart[match(av[,1],test$uni_YS)]
av$Date<-test$Date[match(av[,1],test$uni_YS)]
av$angQ<-angQ$cfs[match(av$Date,angQ$Date)]
av$sanaQ<-sanaQ$cfs[match(av$Date,sanaQ$Date)]
av$year<-test$year[match(av[,1],test$uni_YS)]
av$Cstrata<-test$Cstrata[match(av[,1],test$uni_YS)]
av$cQ<-ifelse(av$Cstrata==3,av$angQ,av$sanaQ)
av$gHab<-0
av$bHab<-0
for (i in 1:length(av[,1])){
  temp<-subset(test,test$uni_YS==av[i,1])
  av$gHab[i]<-sum(subset(temp$effort,temp$cHab<3))
  av$bHab[i]<-sum(subset(temp$effort,temp$cHab>2))
}
av$TotHab<-av$gHab+av$bHab
av$pro<-av$gHab/av$TotHab
av<-subset(av,av$year>2008) ### only use 2009-2011 because 2008 data is incomplete
#### reformat data from Braun et al. 2015
tmeso$Chab<-ifelse(tmeso$mesohab_cl==1|tmeso$mesohab_cl==7|tmeso$mesohab_cl==8|tmeso$mesohab_cl==14|tmeso$mesohab_cl==15,1,3)
tmeso$Cstrata<-findInterval(tmeso$RM,grid)
tmeso$angQ<-angQ$cfs[match(tmeso$Date,angQ$Date)]
tmeso$sanaQ<-sanaQ$cfs[match(tmeso$Date,sanaQ$Date)]
tmeso$cQ<-ifelse(tmeso$Cstrata==3,tmeso$angQ,tmeso$sanaQ)
tmeso$uni_RD<-paste(tmeso$RM,tmeso$Date)
meso<-data.frame(uni_RD=unique(tmeso$uni_RD))
meso$RM<-tmeso$RM[match(meso[,1],tmeso$uni_RD)]
meso$site_abv<-tmeso$site_abv[match(meso[,1],tmeso$uni_RD)]
meso$Date<-tmeso$Date[match(meso[,1],tmeso$uni_RD)]
meso$Cstrata<-tmeso$Cstrata[match(meso[,1],tmeso$uni_RD)]
meso$cQ<-tmeso$cQ[match(meso[,1],tmeso$uni_RD)]
meso$gHab<-0
meso$bHab<-0
for (j in 1:length(meso[,1])){
  temp<-subset(tmeso,tmeso$uni_RD==meso$uni_RD[j])
  meso$gHab[j]<-sum(subset(temp$Shape_Area,temp$Chab==1))
  meso$bHab[j]<-sum(subset(temp$Shape_Area,temp$Chab==3))
}
meso$TotHab<-meso$gHab+meso$bHab
meso$pro<-meso$gHab/meso$TotHab
meso<-subset(meso,meso$Cstrata!=4)
mAV<-rbind(cbind(meso$Cstrata,meso$cQ/1000,meso$gHab,meso$gHab+meso$bHab),
           cbind(av$Cstrata,av$cQ/1000,av$gHab,av$gHab+av$bHab)) 
NAVsamps<-dim(mAV)[1]
#### done reformating and combining mesohabitat availability data

#### further summarize catch data 
# isolate data for age - 1+ fish only
tt<-subset(monC,(monC$hybama01==0&monC$month<10)|monC$month<7)
mon1<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,ifelse(tt$month<6,tt$hybama1+tt$hybama01,tt$hybama1))
mon1_effQ<-cbind(tt$effort,tt$cQ/1000)
colnames(mon1)<-c("year","julian","Cstrata","habitat","catch")
Nobs_mon1<-dim(mon1)[1]

# isolate data for augmented fish only
###lump different vie marks and start at period 49 for coarse strata 1 and 2 and 97 for coarse strata 3...
tt<-subset(mon,mon$Cstrata<3&mon$period>48&mon$jul<60)
monV<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama_vie)
monV_effQ<-cbind(tt$effort,tt$cQ/1000,tt$vieW)
tt<-subset(mon,mon$Cstrata>2&mon$period>96&mon$jul<60)
monV<-rbind(monV,cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama_vie))
monV_effQ<-rbind(monV_effQ,cbind(tt$effort,tt$cQ/1000,tt$vieW))
monV<-subset(monV,is.na(monV_effQ[,3])==FALSE)
monV_effQ<-subset(monV_effQ,is.na(monV_effQ[,3])==FALSE)
Nobs_monV<-dim(monV)[1]
colnames(monV)<-c("year","julian","Cstrata","habitat","catch")

# isolate data for age 0 fish only
tt<-subset(monC,monC$hybama01==0&monC$month>6&monC$month<10)
mon0<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama0)
colnames(mon0)<-c("year","julian","Cstrata","habitat","catch")
mon0_effQ<-cbind(tt$effort,tt$cQ/1000)
Nobs_mon0<-dim(mon0)[1]

#  isolate data where age -0  and  age - 1+ fish were not separated
tt<-subset(monC,monC$hybama01>0&monC$month>6&monC$month<10)
mon01<-cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama0+tt$hybama1+tt$hybama01)
mon01_effQ<-cbind(tt$effort,tt$cQ/1000)
#now do oct and november data -
tt<-subset(monC,monC$month==10)
mon01<-rbind(mon01,cbind((tt$year-2001),tt$jul,tt$Cstrata,tt$type,tt$hybama0+tt$hybama1+tt$hybama01))
mon01_effQ<-rbind(mon01_effQ,cbind(tt$effort,tt$cQ/1000))
# add november data
mon01<-rbind(mon01,cbind((monre$year-2001),215,monre$Cstrata,monre$type,monre$hybama_1))
mon01_effQ<-rbind(mon01_effQ,cbind(monre$effort,monre$cQ/1000))
colnames(mon01)<-c("year","julian","Cstrata","habitat","catch")
Nobs_mon01<-dim(mon01)[1]

##read in results of population estimation drawn from Dudley et al., 2012
Nz<-c(1108430,1387948,267272,122381)
lNz<-log(Nz)
cvz<-c(0.3,0.259,0.372,0.376)
Cz<-exp(1.96*sqrt(log(1+cvz^2)))
lCz<-log(Cz)/1.96

###read in, subset and reformat fish rescue data using drying data
trescue<-subset(trescue,trescue[,1]>2008&trescue[,1]<2019)
trescue$jul<-NA
for (i in 1:length(trescue$jul)){
  trescue$jul[i]<-julian(as.Date(paste(trescue$year[i],trescue$month[i],trescue$dom[i],sep="-")),as.Date(paste(trescue$year[i],"03","31",sep="-")))[[1]]}
trescue$strata<-findInterval(trescue$rm,vgrid)
tR0<-matrix(NA,nrow=(Nyears-7),ncol=Nvstrata)
tR1<-matrix(NA,nrow=(Nyears-7),ncol=Nvstrata)
for (t in 1:(Nyears-7)){
  for (j in 1:Nvstrata){
    if(length(subset(trescue$jul,trescue$year==(t+2008)&trescue$strata==j))>0){
      tday1<-max(c(min(subset(trescue$jul,trescue$year==(t+2008)&trescue$strata==j))+14,fdryday[(t+7),j]+4),na.rm=T)
      if (tday1>0){
        temp<-subset(trescue,trescue$year==(t+2008)&trescue$strata==j&trescue$jul<tday1)
        tR1[t,j]<-sum(as.numeric(paste(temp$adult.alive)),na.rm=T)
        tR0[t,j]<-sum(temp$yoy.alive,na.rm=T)
      }}}}
#
Ntotjul<-215
StrataLen_rm<-array(NA,dim=c(Nyears,Ntotjul,3))
prop_nd<-array(NA,dim=c(Nyears,Ntotjul,3))
cum_nd<-array(NA,dim=c(Nyears,Ntotjul,3))
R0<-matrix(NA,ncol=4,nrow=(Nyears-7)*Ntotjul*Nstrata)
R1<-matrix(NA,ncol=4,nrow=(Nyears-7)*Ntotjul*Nstrata)
R0[,1]<-rep(c(8:Nyears),each=Ntotjul*Nstrata)
R1[,1]<-rep(c(8:Nyears),each=Ntotjul*Nstrata)
R0[,2]<-rep(c(1:Ntotjul),(Nyears-7)*Nstrata)
R1[,2]<-rep(c(1:Ntotjul),(Nyears-7)*Nstrata)
R0[,3]<-rep(rep(c(1:Nstrata),each=Ntotjul),(Nyears-7))
R1[,3]<-rep(rep(c(1:Nstrata),each=Ntotjul),(Nyears-7))
maxStrataLen<-table(crossgrid)/5
cum_phiR<-array(NA,dim=c(Nyears,Ntotjul,3))
# calculate weighted average survival of rescued fish
phiR$jul<-NA	
for (j in 1:12){phiR$jul[j]<-julian(as.Date(paste(2000,phiR$month[j],phiR$day[j],sep="-")),as.Date("2000-03-31"))[[1]]}
tm<-glm(cbind(phiR[,4],phiR[,3]-phiR[,4])~phiR[,5],family="binomial")
pred_phiR<-plogis(coef(tm)[1]+coef(tm)[2]*c(1:Ntotjul))
for (t in 1:Nyears){
  for (d in 1:Ntotjul){
    StrataLen_rm[t,d,1]<-length(which(crossgrid==1&(is.na(fdryday[t,])==T|fdryday[t,]>d|ldryday[t,]<d)))/5
    StrataLen_rm[t,d,2]<-length(which(crossgrid==2&(is.na(fdryday[t,])==T|fdryday[t,]>d|ldryday[t,]<d)))/5
    StrataLen_rm[t,d,3]<-length(which(crossgrid==3&(is.na(fdryday[t,])==T|fdryday[t,]>d|ldryday[t,]<d)))/5
    t1<-which(crossgrid==1&fdryday[t,]==d)
    t2<-which(crossgrid==2&fdryday[t,]==d)
    t3<-which(crossgrid==3&fdryday[t,]==d)
    prop_nd[t,d,1]<-length(t1)/5/maxStrataLen[1]
    prop_nd[t,d,2]<-length(t2)/5/maxStrataLen[2]
    prop_nd[t,d,3]<-length(t3)/5/maxStrataLen[3]
    if (t>7&length(t1)>0){
      R0[which(R0[,1]==t&R0[,2]==d&R0[,3]==1),4]<-sum(tR0[(t-7),t1])
      R1[which(R1[,1]==t&R1[,2]==d&R1[,3]==1),4]<-sum(tR1[(t-7),t1])
    }
    if (t>7&length(t2)>0){
      R0[which(R0[,1]==t&R0[,2]==d&R0[,3]==2),4]<-sum(tR0[(t-7),t2])
      R1[which(R1[,1]==t&R1[,2]==d&R1[,3]==2),4]<-sum(tR1[(t-7),t2])
    }
    if (t>7&length(t3)>0){
      R0[which(R0[,1]==t&R0[,2]==d&R0[,3]==3),4]<-sum(tR0[(t-7),t3])
      R1[which(R1[,1]==t&R1[,2]==d&R1[,3]==3),4]<-sum(tR1[(t-7),t3])
    }
  }
  cum_nd[t,,1]<-cumsum(prop_nd[t,,1])
  cum_nd[t,,2]<-cumsum(prop_nd[t,,2])
  cum_nd[t,,3]<-cumsum(prop_nd[t,,3])
  cum_phiR[t,,1]<-cumsum(prop_nd[t,,1]*pred_phiR)
  cum_phiR[t,,2]<-cumsum(prop_nd[t,,2]*pred_phiR)
  cum_phiR[t,,3]<-cumsum(prop_nd[t,,3]*pred_phiR)
}
R0<-subset(R0,is.na(R0[,4])==FALSE&R0[,2]>91)
R1<-subset(R1,is.na(R1[,4])==FALSE)
NR0<-length(R0[,1])
NR1<-length(R1[,1])
# convert to rkm
StrataLen<-array(NA,dim=c(Nyears,Ntotjul,3))
StrataLen[,,1]<-92.1*StrataLen_rm[,,1]/maxStrataLen[1]
StrataLen[,,2]<-85.5*StrataLen_rm[,,2]/maxStrataLen[2]
StrataLen[,,3]<-65*StrataLen_rm[,,3]/maxStrataLen[3]

## reformat summarize expert ellicitation data
ee<-list(ee1=ee1,ee2=ee2,ee3=ee3,ee4=ee4,ee5=ee5)
# function to calculate standard error from quantile, upper, lower and mean info
calcsig<-function(up,down,mean,q,INT=c(0,1000)){
  sig<-function(x){abs(pnorm(up,mean,x)-pnorm(down,mean,x)-q)}
  optimize(sig,interval=INT)$minimum}
## 5a - create prior on movement out of reach
#extract and convert info to sd
emove<-matrix(NA,nrow=5,ncol=2)
for (i in 1:5){
  temp<-ee[[i]][ee[[i]][,1]=="5a",]
  emove[i,1]<-temp[1,6]
  emove[i,2]<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
}
## 6a,6b,6c - create priors on river width
#extract and convert info to sd
Nexperts<-5
ewidths<-array(NA,dim=c(5,3,2,8))
for (i in 1:5){
  temp<-ee[[i]][ee[[i]][,1]=="6a - river width",]
  temp2<-ee[[i]][ee[[i]][,1]=="6b - river width",]
  temp3<-ee[[i]][ee[[i]][,1]=="6c - river width",]
  for (j in 1:8){
    ewidths[i,1,1,j]<-temp[j,6]
    ewidths[i,2,1,j]<-temp2[j,6]
    ewidths[i,3,1,j]<-temp3[j,6]
    ewidths[i,1,2,j]<-calcsig(temp[j,4],temp[j,3],temp[j,6],temp[j,5]/100)
    ewidths[i,2,2,j]<-calcsig(temp2[j,4],temp2[j,3],temp2[j,6],temp2[j,5]/100)
    ewidths[i,3,2,j]<-calcsig(temp3[j,4],temp3[j,3],temp3[j,6],temp3[j,5]/100)
  }}
refQ<-c(5,50,100,150,200,250,500,1000)/1000

save(Nyears, Nstrata, NAVsamps, Nobs_mon1, Nobs_mon0, Nobs_mon01, Nobs_monV, ewidths,
                    Nexperts, emove, refQ, mon1, mon0, mon01, monV, StrataLen, lNz, lCz, R0, R1,
                    cum_nd, Ntotjul, cum_phiR, NR0, NR1, w_um, mAV, nmons, w_m, mon1_effQ, mon01_effQ,
                    mon0_effQ, monV_effQ, prop_nd, ee, sQ, aQ, sanaQ, angQ,
                    file = "output/input_data.RData")

# Save separately for use in 3_scenarios
saveRDS(ee, "output/ee_list.RData")

####################
#### Prepare out of sample data for model evaluation
#####################

# This code is new to this analysis although it uses some of the out-of-sample code from Yackulic et al. 2022

library(tidyverse)

## Original out of sample data
oosD_orig<-read.csv("data/yackulic2022_data/oos_rivereyes.csv")
oosR_orig<-read.csv("data/yackulic2022_data/oos_releases.csv")
oosc_orig <-read.csv("data/yackulic2022_data/oos_catch.csv") 

# Format for just necessary columns, to see what's needed
# Release data
vkeep2<-match(c("Date","Month","Year","C","L","S","Number","RM"),names(oosR_orig))
oosR_orig<-oosR_orig[,vkeep2] %>%
  mutate(Date = as.Date(Date, "%d-%b-%y"))

# Catch
keep<-match(c("Date.Collected","RM_Start","Haul","Habitat","Effort_m.2","year","month",
              "Species_Codes","Gear","Length..SL.mm.","AgeClass","SumOfSPEC",
              "Sampling_Period"),names(oosc_orig))
oosc_orig<-oosc_orig[,keep]

### New Data

# Release data
#Downloaded from:https://data.mendeley.com/datasets/nwc7k6rm47/6
new_oosd <- read_csv("data/new_oos/dry_eyes_download_reyes.csv") %>%
  mutate(date = as.Date(Date, "%m/%d/%Y"),
         Year = as.numeric(format(date, "%Y")),
         month = as.numeric(format(date, "%m")),
         dom = as.numeric(format(date, "%d"))) %>%
  rename(URM = `Upstream Dry River Mile`, LRM = `Downstream Dry River Mile`, Distance = `Dry Length (River Miles)`) %>%
  select(names(oosD_orig)) %>%
  filter(Year > 2018)

#Order of new_oosd data differs from old oosd data, but appears to be identical based on mergeR 
new_oosr <- read_csv("data/new_oos/oosR_raw.csv") %>%
  mutate(Date = as.Date(Date, "%m/%d/%y"),
         Month = as.numeric(format(Date, "%m")),
         RM = as.numeric(RM)) %>%
  filter(Date > "2018-09-30") %>%
  select(names(oosR_orig)) %>%
  filter(!is.na(C))

oosD <- new_oosd
oosR <- new_oosr

#### Out of sample predictions
### format out of sample data
oosD$ustrata<-findInterval(oosD$URM,vgrid)
oosD$dstrata<-findInterval(oosD$LRM,vgrid)
oosD$jul<-NA
for (i in 1:length(oosD$jul)){
  oosD$jul[i]<-julian(as.Date(paste(oosD$Year[i],oosD$month[i],oosD$dom[i],sep="-")),as.Date(paste(oosD$Year[i],"03","31",sep="-")))[[1]]}
fdryday_oos<-matrix(NA,nrow=6,ncol=(length(vgrid)-1))
ldryday_oos<-matrix(NA,nrow=6,ncol=(length(vgrid)-1))
for (k in 1:(length(vgrid)-1)){
  for (j in 1:6){
    temp2<-subset(oosD$jul,k<=oosD$ustrata&k>=oosD$dstrata&oosD$Year==(2018+j))
    fdryday_oos[j,k]<-ifelse(length(temp2)==0,NA,min(temp2))
    ldryday_oos[j,k]<-ifelse(length(temp2)==0,NA,max(temp2))
  }}
StrataLen_rm_oos<-array(NA,dim=c(6,Ntotjul,3))
prop_nd_oos<-array(NA,dim=c(6,Ntotjul,3))
cum_nd_oos<-array(NA,dim=c(6,Ntotjul,3))
cum_phiR_oos<-array(NA,dim=c(6,Ntotjul,3))
for (t in 1:6){
  for (d in 1:Ntotjul){
    StrataLen_rm_oos[t,d,1]<-length(which(crossgrid==1&(is.na(fdryday_oos[t,])==T|fdryday_oos[t,]>d|ldryday_oos[t,]<d)))/5
    StrataLen_rm_oos[t,d,2]<-length(which(crossgrid==2&(is.na(fdryday_oos[t,])==T|fdryday_oos[t,]>d|ldryday_oos[t,]<d)))/5
    StrataLen_rm_oos[t,d,3]<-length(which(crossgrid==3&(is.na(fdryday_oos[t,])==T|fdryday_oos[t,]>d|ldryday_oos[t,]<d)))/5
    t1<-which(crossgrid==1&fdryday_oos[t,]==d)
    t2<-which(crossgrid==2&fdryday_oos[t,]==d)
    t3<-which(crossgrid==3&fdryday_oos[t,]==d)
    prop_nd_oos[t,d,1]<-length(t1)/5/maxStrataLen[1]
    prop_nd_oos[t,d,2]<-length(t2)/5/maxStrataLen[2]
    prop_nd_oos[t,d,3]<-length(t3)/5/maxStrataLen[3]
  }
  cum_nd_oos[t,,1]<-cumsum(prop_nd_oos[t,,1])
  cum_nd_oos[t,,2]<-cumsum(prop_nd_oos[t,,2])
  cum_nd_oos[t,,3]<-cumsum(prop_nd_oos[t,,3])
  cum_phiR_oos[t,,1]<-cumsum(prop_nd_oos[t,,1]*pred_phiR)
  cum_phiR_oos[t,,2]<-cumsum(prop_nd_oos[t,,2]*pred_phiR)
  cum_phiR_oos[t,,3]<-cumsum(prop_nd_oos[t,,3]*pred_phiR)
}
StrataLen_oos<-array(NA,dim=c(6,Ntotjul,3))
StrataLen_oos[,,1]<-92.1*StrataLen_rm_oos[,,1]/maxStrataLen[1]
StrataLen_oos[,,2]<-85.5*StrataLen_rm_oos[,,2]/maxStrataLen[2]
StrataLen_oos[,,3]<-65*StrataLen_rm_oos[,,3]/maxStrataLen[3]
#
#Original code commented out because these steps have already been done
#vkeep2<-match(c("Date","Month","Year","C","L","S","Number","RM"),names(oosR))
#oosR<-oosR[,vkeep2]
oosR$strata<-findInterval(converter[match(oosR$RM,converter[,1]),2],vgrid)
oosR$period<-ifelse(oosR$Month<4,8*(oosR$Year-2002)+1,ifelse(oosR$Month>9,8*(oosR$Year-2002)+9,8*(oosR$Year-2002)+oosR$Month-3))
oosR$winmon<-ifelse(oosR$Month<4,4-oosR$Month,ifelse(oosR$Month>9,16-oosR$Month,0))
#
tw_m_oos<-array(0,dim=c(6,Nvstrata,2)) #summarized at fine grid
w_m_oos<-array(0,dim=c(6,Nstrata,2)) #summarized at coarse grid
#
for (i in 1:6){
  for (j in 1:Nvstrata){
    temp<-subset(oosR,((oosR$period-1)/8-16)==i&oosR$strata==j)
    if (dim(temp)[1]>0){
      tw_m_oos[i,j,1]<-sum(subset(temp$Number,temp$winmon<3))
      tw_m_oos[i,j,2]<-sum(subset(temp$Number,temp$winmon>3))
    }}
  for (j in 1:Nstrata){
    w_m_oos[i,j,1]<-sum(tw_m_oos[i,gridst_end[j,1]:gridst_end[j,2],1])
    w_m_oos[i,j,2]<-sum(tw_m_oos[i,gridst_end[j,1]:gridst_end[j,2],2])
  }}


# Need 2019 catch data from oosc original

oosc_2019 <- oosc_orig %>%
  filter(Sampling_Period == 201910) 

# New ASIR data
raw_oosc <- read_csv("data/new_oos/PopMon_MonthlyHaulUSBR.csv")

oosc <- raw_oosc %>%
  mutate(date = as.Date(`Date Collected`, "%m/%d/%y"),
         year = as.numeric(format(date, "%Y")),
         month = as.numeric(format(date, "%m")),
         day = as.numeric(format(date, "%d")),
         Date.Collected =  paste(month, day, year, sep = "/")) %>%
  filter(Sampling_Period %in% c(202010, 202110, 202210, 202310, 202410)) %>%
  select(Date.Collected, 
         RM_Start, 
         Haul,
         Habitat,
         Effort_m.2 =`Effort_m^2`, 
         year, month,
         Species_Codes, 
         Gear,
         Length..SL.mm. = `Length (SL,mm)`,
         AgeClass,
         SumOfSPEC) %>%
  bind_rows(oosc_2019) %>%
  select(-Sampling_Period) # This had to be added initially to filter just the appropriate 2019 and subsequent data

#Need additional flow data

#Angostura reach
start.date <- "2021-01-01"
end.date <- "2024-12-31"
siteAng <- "08330000"
siteSan <- "08354900"
pCode <- "00060"

library(dataRetrieval)

new_ang_data <- readNWISdv(siteNumbers = siteAng,
                           parameterCd = pCode,
                           startDate = start.date,
                           endDate = end.date) %>%
  mutate(year = as.numeric(format(Date, "%Y")),
         month = as.numeric(format(Date, "%m")),
         day = as.numeric(format(Date, "%d")),
         Date = paste(month, day, year, sep = "/")) %>%
  rename(cfs = X_00060_00003) %>%
  select(names(angQ))

new_SanA_data <- readNWISdv(siteNumbers = siteSan,
                            parameterCd = pCode,
                            startDate = start.date,
                            endDate = end.date)%>%
  mutate(year = as.numeric(format(Date, "%Y")),
         month = as.numeric(format(Date, "%m")),
         day = as.numeric(format(Date, "%d")),
         Date = paste(month, day, year, sep = "/")) %>%
  rename(cfs = X_00060_00003) %>%
  select(names(sanaQ))

angQ_all <- bind_rows(angQ, new_ang_data)
sanaQ_all <- bind_rows(sanaQ, new_SanA_data)

# Re-start original oos code; change 
oosc<-subset(oosc,oosc$Habitat!=""&oosc$Gear!="larval")
ghab<-c("","BW","MCPLPO","MCPO","MCSHPLPO","MCSHPO","PO","SCPLPO","SCPO","SCSHPLPO","SCSHPO","SHPO","MCED","MCSHED","SCED","SCSHED")
oosc$cHab<-ifelse(is.na(match(oosc$Habitat,ghab))==FALSE,1,3)
oosc$strata<-findInterval(oosc$RM_Start,vgrid)
oosc$period<-8*(oosc$year-2002)+oosc$month-3
extractdom<-function(x){
  lx<-nchar(x)
  tx<-numeric()
  for (i in 1:lx){
    tx[i]<-substr(x,i,i)}
  wx<-which(tx=="/")
  as.numeric(substr(x,wx[1]+1,wx[2]-1))}
oosc$dom<-NA
oosc$jul<-NA
oosc$dry<-NA
#for (j in 1:length(oosc[,1])){
for (j in 1:nrow(oosc)){
  oosc$dom[j]<-extractdom(paste(oosc$Date.Collected[j]))
  oosc$jul[j]<-julian(as.Date(paste(oosc$year[j],oosc$month[j],oosc$dom[j],sep="-")),as.Date(paste(oosc$year[j],"03","31",sep="-")))[[1]]
  tf1<-fdryday_oos[(oosc$year[j]-2018),oosc$strata[j]]
  tl1<-ldryday_oos[(oosc$year[j]-2018),oosc$strata[j]]
  oosc$dry[j]<-ifelse(is.na(tf1)==TRUE|(tf1>oosc$jul[j])|(tl1<oosc$jul[j]),0,1)}
oosc$uni_id<-paste(oosc$year,oosc$jul,oosc$RM_Start,oosc$Haul)
#
toosc<-data.frame(uni=sort(unique(oosc$uni_id)))
toosc$RMStart<-oosc$RM_Start[match(toosc[,1],oosc$uni_id)]
toosc$period<-oosc$period[match(toosc[,1],oosc$uni_id)]
toosc$jul<-oosc$jul[match(toosc[,1],oosc$uni_id)]
toosc$year<-oosc$year[match(toosc[,1],oosc$uni_id)]
toosc$month<-oosc$month[match(toosc[,1],oosc$uni_id)]
toosc$date<-oosc$Date.Collected[match(toosc[,1],oosc$uni_id)]
toosc$angQ<-angQ_all$cfs[match(toosc$date,angQ_all$Date)]
toosc$sanaQ<-sanaQ_all$cfs[match(toosc$date,sanaQ_all$Date)]
toosc$HaulNo<-oosc$Haul[match(toosc[,1],oosc$uni_id)]
toosc$type<-ifelse(oosc$cHab[match(toosc[,1],oosc$uni_id)]==1,1,2)
toosc$effort<-as.numeric(paste(oosc$Effort_m.2[match(toosc[,1],oosc$uni_id)]))
toosc$strata<-oosc$strata[match(toosc[,1],oosc$uni_id)]
toosc$cQ<-ifelse(toosc$strata>573,toosc$angQ,toosc$sanaQ)
toosc$hybama<-0
thybama<-subset(oosc,oosc$Species_Codes=="HYBAMA")
for (i in 1:length(toosc[,1])){
  temp<-subset(thybama,thybama$uni_id==toosc$uni[i])
  toosc$hybama[i]<-sum(temp$SumOfSPEC)}
toosc$Cstrata<-crossgrid[toosc$strata]
toosc$c_id<-paste(toosc$year,toosc$jul,toosc$Cstrata,toosc$type)
##
oos_C<-data.frame(cid=sort(unique(toosc$c_id)))
oos_C$jul<-toosc$jul[match(oos_C[,1],toosc$c_id)]
oos_C$year<-toosc$year[match(oos_C[,1],toosc$c_id)]
oos_C$month<-toosc$month[match(oos_C[,1],toosc$c_id)]
oos_C$Cstrata<-toosc$Cstrata[match(oos_C[,1],toosc$c_id)]
oos_C$type<-toosc$type[match(oos_C[,1],toosc$c_id)]
oos_C$effort<-0
oos_C$cQ<-0
oos_C$hybama<-0
oos_C$period<-toosc$period[match(oos_C[,1],toosc$c_id)]
mean_names<-c("cQ")
sum_names<-c("effort","hybama")
tsm<-match(sum_names,names(toosc))
sm<-match(sum_names,names(oos_C))
tmn<-match(mean_names,names(toosc))
smn<-match(mean_names,names(oos_C))
for (i in 1:length(oos_C[,1])){
  temp<-subset(toosc[,tsm],toosc$c_id==oos_C$cid[i]) #subset data to catch and effort in a single survey
  oos_C[i,sm]<-colSums(temp) #sum these values
  temp<-subset(toosc[,tmn],toosc$c_id==oos_C$cid[i]) #subset to discharge
  oos_C[i,smn]<-mean(temp) # mean discharge
}

Nobs_oosC<-dim(oos_C)[1]

save(cum_nd_oos, cum_phiR_oos, Nobs_oosC, Nstrata, oos_C, StrataLen_oos, w_m_oos,
     file = "output/oos_data_new.RData")
