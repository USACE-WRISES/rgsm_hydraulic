library(tidyverse)
library(rstan)
library(zoo)
library(hrbrthemes) # For plotting
library(dataRetrieval)

# Read in data
load(file = "output/input_data.RData")

# Create functions to calculate larval carrying capacity
## Helper functions ----

# Calculate standard error for expert informed covariates from quantile, upper, lower, mean info
calcsig <- function(up, down, mean, q, INT = c(0, 1000)) {
  sig <- function(x) {
    abs(pnorm(up, mean, x) - pnorm(down, mean, x) - q)
  }
  optimize(sig, interval = INT)$minimum
}

# Linear interpolation (used for Q-habitat relationship)
lin_int <- function(x, xlow, xhi, ylow, yhi) {
  (yhi - ylow) * (x - xlow) / (xhi - xlow) + ylow
} 

## Functions to calculate expert informed covariate ----

#Amount of habitat

#Original function
q2hab<-function(q,pars){
  p<-c(0,pars,pars[length(pars)])
  qs<-c(0,5,50,100,150,200,250,1000,1500,2000,2500,3000,4000,5000,6000,7000,2*10^4)
  t1<-findInterval(q,qs)
  lin_int(q,qs[t1],qs[(t1+1)],p[t1],p[(t1+1)])
}

# Function modified for eco-value curves
q2hab_2d <- function(q, pars) {
  p <- c(pars, pars[length(pars)])
  qs <- c(0, 700, 2000, 3000, 4000, 5000, 6000,7000, 2 * 10 ^ 4)
  t1 <- findInterval(q, qs)
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}

# Function modified for eco-value curves
q2hab_combined <- function(q, pars) {
  p <- c(pars, pars[length(pars)])
  qs <- c(0, 5, 50,  100,  150,  200, 1000, 1500,  2000, 3000, 4000, 5000, 6000,7000, 2 * 10 ^ 4)
  t1 <- findInterval(q, qs)
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}


#Calculate proportion of eggs ready on specific days based on the expert elicited values 
#Calculates it for all days in spawning window (March 1 = Day 1 through July 11)
prop <- function(pars) {
  t <- seq(2, 132, 10)
  tout <- pars[1]
  for (i in 2:131) {
    t1 <- findInterval(i, t)
    tout[i] <- lin_int(i, t[t1], t[(t1 + 1)], pars[t1], pars[(t1 + 1)])
  }
  tout[132] <- pars[14]
  tout / sum(tout)
}

#Calculate the larval carrying capacity from Yackulic et al.

# prPARS = Proportion of eggs laid
# kappa  = spawning cue
# q2fpars = larval habitat
# D = Duration

#Original function
calc_cov<-function(Q,q2fPARS,kappa,D,prPARS){
  thab<-numeric()
  for (i in 1:214){thab[i]<-q2hab(Q[i],pars=q2fPARS)} #finds rel. amount of habitat on given day, based on flow-hab "curve"
  tprop<-prop(prPARS) #proportional egg laying on each day based on expert elicitation; days are interpolated
  t2hab<-numeric()
  for (i in 1:132){t2hab[i]<-min(thab[i:(i+D)])} #Find minimum habitat within required Duration (D) on each day in spawning window
  tstart<-min(c(which(Q>kappa), #date when spawning is cued; first date where:q>flow cue 
                which(Q[-1]-Q[-length(Q)]>100), #discharge changes by 100 cfs in a day
                132)) #or last day of spawning window
  out<-sum(tprop[1:tstart])*t2hab[tstart]+ #Proportion of eggs available pre-spawning *amount of habitat on first day
    sum(tprop[tstart:132]*t2hab[tstart:132]) #Proportion of eggs on each subsequent day * amount of hab. on each subsequent day through day 132  
  return(out)}

#Altered function to use the q2hab_2d function
calc_cov_2d <- function(Q, q2fPARS, kappa, D, prPARS) {
  thab <- numeric()
  for (i in 1:214) {
    thab[i] <- log(q2hab_2d(Q[i], pars = q2fPARS) + 1)
  } #finds rel. amount of habitat on given day, based on flow-hab "curve"
  tprop <- prop(prPARS) #proportional egg laying on each day based on expert elicitation; days are interpolated
  t2hab <- numeric()
  for (i in 1:132) {
    t2hab[i] <- min(thab[i:(i + D)])
  } #Find minimum habitat within required Duration (D) on each day in spawning window
  tstart <- min(c(which(Q > kappa), #date when spawning is cued; first date where:q>flow cue
                  which(Q[-1] - Q[-length(Q)] > 100), #discharge changes by 100 cfs in a day
                  132)) #or last day of spawning window
  out <- sum(tprop[1:tstart]) * t2hab[tstart] + #Proportion of eggs available pre-spawning *amount of habitat on first day
    sum(tprop[tstart:132] * t2hab[tstart:132]) #Proportion of eggs on each subsequent day * amount of hab. on each subsequent day through day 132
  return(out)
}

#Altered function to calculate composite covariate function
calc_cov_combined <- function(Q, q2fPARS, kappa, D, prPARS) {
  thab <- numeric()
  for (i in 1:214) {
    thab[i] <- log(q2hab_combined(Q[i], pars = q2fPARS) + 1)
  } #finds rel. amount of habitat on given day, based on flow-hab "curve"
  tprop <- prop(prPARS) #proportional egg laying on each day based on expert elicitation; days are interpolated
  t2hab <- numeric()
  for (i in 1:132) {
    t2hab[i] <- min(thab[i:(i + D)])
  } #Find minimum habitat within required Duration (D) on each day in spawning window
  tstart <- min(c(which(Q > kappa), #date when spawning is cued; first date where:q>flow cue
                  which(Q[-1] - Q[-length(Q)] > 100), #discharge changes by 100 cfs in a day
                  132)) #or last day of spawning window
  out <- sum(tprop[1:tstart]) * t2hab[tstart] + #Proportion of eggs available pre-spawning *amount of habitat on first day
    sum(tprop[tstart:132] * t2hab[tstart:132]) #Proportion of eggs on each subsequent day * amount of hab. on each subsequent day through day 132
  return(out)
}

######################	
#### STEP X: Calculate covariates for different models
######################	

#t1 = prPARS = Proportion of eggs laid
#t2 = kappa  = spawning cue
#t3 = q2fpars = larval habitat
#t4 = D = Duration

# Original larval carry capacity indices (based on expert 3)
preds2<-array(NA,dim=c(17,3))
riversegment<-c("3e","3c","3a")

temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
t1<-numeric()
for (t in 1:6){
  se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
  t1[t]<-temp[t,6]-se
  se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
  t1[(t+7)]<-temp[(t+6),6]+se
}
t1[7]<-1
t1[14]<-0
temp<-subset(ee[[3]],ee[[3]][,1]=="2")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t2<-temp[1,6]
temp<-subset(ee[[3]],ee[[3]][,1]=="4")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t4<-temp[1,6]-se
#
for (r in 1:3){
  temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
  t3<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t3[t]<-temp[t,6]+se
  }	
  for (t in 8:15){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t3[t]<-temp[t,6]-se
  }
  t3[7]<-temp[7,6]
  for (j in 1:17){
    if (r==3) {q<-aQ[,j]} else {q<-sQ[,j]}
    preds2[j,r]<-calc_cov(q,t3,t2,t4,t1)
  }}

## Covariates using average flow in May-June
mjflow_sana<-subset(sanaQ,(sanaQ$month==5|sanaQ$month==6)&sanaQ$year>2001&sanaQ$year<2019)
mjflow_ang<-subset(angQ,(angQ$month==5|angQ$month==6)&angQ$year>2001&angQ$year<2019)
simpflow<-cbind(tapply(mjflow_sana$cfs,mjflow_sana$year,mean),tapply(mjflow_sana$cfs,mjflow_sana$year,mean),tapply(mjflow_ang$cfs,mjflow_ang$year,mean))/1000

## Covariate using inundated acres

q_inund_curve <- read_csv("data/ecoval2d.csv") %>% # 2D
  select(q, inund_lower = tot_outside_chan_lower, 
         inund_upper = tot_outside_chan_upper, wua_lower, wua_upper,reach, reach_num)

reaches <-c("San Acacia", "Isleta", "Angostura")

q_hab_lookup <- data.frame()

for (i in 1:length(reaches)){
  q_inund_i <- filter(q_inund_curve, reach == reaches[i])
  
  hab_lookup_i <- tibble(cfs = 1:7e3,
                         reach = rep(reaches[i], 7e3)) %>%
    mutate(inund_lower = q2hab_2d(cfs, q_inund_i$inund_lower),
           inund_upper = q2hab_2d(cfs, q_inund_i$inund_upper),
           wua_lower = q2hab_2d(cfs, q_inund_i$wua_lower),
           wua_upper = q2hab_2d(cfs, q_inund_i$wua_upper))
  
  q_hab_lookup <- bind_rows(q_hab_lookup, hab_lookup_i)
  
}

Qdur = 20 # Elicited value of kappa, rounded to nearest whole number

inund_cov_list <- list()

for (i in 1:4) {
  maxmins <- array(NA, dim = c(17,3))
  
  for (j in 1:3) {
    if (j == 1) {
      q <- filter(sanaQ, year %in% 2002:2018)
    } else {
      q <- filter(angQ, year %in% 2002:2018)
    }
    
    hab_lookup_ij <- filter(q_hab_lookup, 
                            reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup)[i+2])
      
    maxmins[,j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
    
  }
  inund_cov_list[[i]] <- maxmins
  #names(inund_cov_list[i]) <- names(q_hab_lookup)[i+2]
}

## Covariate using combined inundated acres/elicited values

elicited_vals <- data.frame()

for (r in 1:3){
  temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
  t3<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t3[t]<-temp[t,6]+se
  }	
  for (t in 8:15){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t3[t]<-temp[t,6]-se
  }
  t3[7]<-temp[7,6]
  
  temp$best <- t3
  temp$reach <- reaches[r]
 
  elicited_vals <- bind_rows(elicited_vals, temp)
   
}

ee_vals_sub <- select(elicited_vals,
                      q=X,
                      best,
                      reach) %>%
  filter(!q %in% c(250, 2500)) %>%
  mutate(q = as.numeric(q))

q_inund_curve_mod <- q_inund_curve %>%
  select(-reach_num) %>%
  filter(!q %in% c(700)) %>%
  full_join(ee_vals_sub) %>%
  arrange(reach, q) %>%
  pivot_longer(inund_lower:wua_upper, names_to  = "inund_type", values_to  = "hab_orig")

inund_ratio <- filter(q_inund_curve_mod, q == 2000) %>%
  select(-c(q)) %>%
  mutate(ratio = hab_orig/best) %>%
  select(reach, inund_type, ratio)

q_inund_curve_combined <- q_inund_curve_mod  %>%
  right_join(inund_ratio) %>%
  mutate(log_hab_orig = case_when(is.na(hab_orig) ~ NA,
                                  TRUE ~ log(hab_orig)),
         hab = case_when(is.na(hab_orig) ~ ratio*best,
                         TRUE ~ hab_orig)) %>%
  arrange(reach, inund_type, q) %>%
  select(q, inund_type, hab, reach) %>% 
  pivot_wider(names_from = inund_type, values_from = hab)
  
q_hab_lookup_combined <- data.frame()

for (i in 1:length(reaches)){
  q_inund_i <- filter(q_inund_curve_combined, reach == reaches[i])
  
  hab_lookup_i <- tibble(cfs = 1:10e3,
                         reach = rep(reaches[i], 10e3)) %>%
    mutate(inund_lower = q2hab_combined(cfs, q_inund_i$inund_lower),
           inund_upper = q2hab_combined(cfs, q_inund_i$inund_upper),
           wua_lower = q2hab_combined(cfs, q_inund_i$wua_lower),
           wua_upper = q2hab_combined(cfs, q_inund_i$wua_upper))
  
  q_hab_lookup_combined <- bind_rows(q_hab_lookup_combined, hab_lookup_i)
  
}

#Save hab_lookup_combined for output scenarios
#write_csv(q_hab_lookup_combined, "data/q_hab_lookup_combined.csv")

Qdur = 20 # Elicited value of kappa, rounded to nearest whole number

inund_cov_combined_list <- list()

for (i in 1:4) {
  maxmins <- array(NA, dim = c(17,3))
  
  for (j in 1:3) {
    if (j == 1) {
      q <- filter(sanaQ, year %in% 2002:2018)
    } else {
      q <- filter(angQ, year %in% 2002:2018)
    }
    
    hab_lookup_ij <- filter(q_hab_lookup_combined, 
                            reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup_combined)[i+2])
    
    maxmins[,j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
    
  }
  inund_cov_combined_list[[i]] <- maxmins
  #names(inund_cov_list[i]) <- names(q_hab_lookup)[i+2]
}

# Model that incorporates inundation into larval carrying capacity

# calculate new set of larval carry capacity indices based on inundation

inund_lcc_list <- list()

for (i in 1:4){

preds<-array(NA,dim=c(17,3))

temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
t1<-numeric()
for (t in 1:6){
  se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
  t1[t]<-temp[t,6]-se
  se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
  t1[(t+7)]<-temp[(t+6),6]+se
}
t1[7]<-1
t1[14]<-0
temp<-subset(ee[[3]],ee[[3]][,1]=="2")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t2<-temp[1,6]
temp<-subset(ee[[3]],ee[[3]][,1]=="4")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t4<-temp[1,6]-se
#
for (r in 1:3){
  t3 <- filter(q_inund_curve, reach_num == r) %>%
    #mutate(inund_acres = log(inund_acres+1)) %>%
    pull(i+1)
  for (t in 1:17){
    if (r==3) {q<-aQ[,t]} else {q<-sQ[,t]}
    preds[t,r]<-calc_cov_2d(q,t3,t2,t4,t1)
  }}

inund_lcc_list[[i]] <- preds

}

# Model that incorporates combined inundation into larval carrying capacity

# calculate new set of larval carry capacity indices based on inundation

#rearrange q_inund_combined to match q_inund_curve

q_inund_combined_format <- q_inund_curve_combined %>%
  relocate(reach, .after = wua_upper) %>%
  mutate(reach_num = rep(c(3,2,1), each = 14))

#write_csv(q_inund_combined_format, "data/q_inund_combined.csv")

inund_lcc_list_combined <- list()

for (i in 1:4){
  
  preds<-array(NA,dim=c(17,3))
  
  temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
  t1<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t1[t]<-temp[t,6]-se
    se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
    t1[(t+7)]<-temp[(t+6),6]+se
  }
  t1[7]<-1
  t1[14]<-0
  temp<-subset(ee[[3]],ee[[3]][,1]=="2")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t2<-temp[1,6]
  temp<-subset(ee[[3]],ee[[3]][,1]=="4")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t4<-temp[1,6]-se
  #
  for (r in 1:3){
    t3 <- filter(q_inund_combined_format, reach_num == r) %>%
      #mutate(inund_acres = log(inund_acres+1)) %>%
      pull(i+1)
    for (t in 1:17){
      if (r==3) {q<-aQ[,t]} else {q<-sQ[,t]}
      preds[t,r]<-calc_cov_combined(q,t3,t2,t4,t1)
    }}
  
  inund_lcc_list_combined[[i]] <- preds
  
}

# Fit candidate models using different covariates ####

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
ct<-function(x){x-mean(x)}

# Original LCC model
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
                  ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,
                  X=ct(preds2),
                  mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
                  w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)

M2_1re <- stan("code/int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

# Original springflows model
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
                  ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,
                  X=ct(simpflow),
                  mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
                  w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)

M2_2re <- stan("code/int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

# Raw Inundation models
# Inund lower
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
                  ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(log(inund_cov_list[[1]]+1)),
                  mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
                  w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)

m_inund_i_l <- stan("code/int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

#Inundation - combined cov. models
# Inund lower
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
                  ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(log(inund_cov_combined_list[[1]]+1)),
                  mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
                  w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)

m_inund_i_l_combined <- stan("code/int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

# LCC Inundation models
# LCC inund lower
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
                  ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(inund_lcc_list[[1]]),
                  mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
                  w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)

m_lcc_i_l <- stan("code/int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

# LCC Inundation Combined models
# LCC
rgsm.data <- list(Nyears=Nyears,Nstrata=Nstrata,NAVsamps=NAVsamps,Nobs_mon1=Nobs_mon1,Nobs_mon0=Nobs_mon0,Nobs_mon01=Nobs_mon01,Nobs_monV=Nobs_monV,
                  ewidths=ewidths, Nexperts=Nexperts,emove=emove,refQ=refQ,X=ct(inund_lcc_list_combined[[1]]),
                  mon1=mon1,mon0=mon0,mon01=mon01,monV=monV,StrataLen=StrataLen,lNz=lNz,lCz=lCz,R0=R0,R1=R1,cum_nd=cum_nd,Ntotjul=Ntotjul,cum_phiR=cum_phiR,NR0=NR0,NR1=NR1,
                  w_um=w_um,mAV=mAV,nmons=nmons,w_m=w_m,mon1_effQ=mon1_effQ,mon01_effQ=mon01_effQ,mon0_effQ=mon0_effQ,monV_effQ=monV_effQ,prop_nd=prop_nd)

m_lcc_i_l_combined <- stan("code/int3re.stan",data = rgsm.data,chains = 3, iter = 10000) 

# Create and save output data ####

# Load select models, if needed
M2_1re <- readRDS("output/full/orig_lcc.rds")
M2_2re <- readRDS("output/full/springflow.rds")
m_inund_i_l <- readRDS("output/full/inund_low.rds")
m_inund_i_l_combined <- readRDS("output/full/inund_low_combined.rds")
m_lcc_i_l <- readRDS("output/full/lcc_low.rds")
m_lcc_i_l_combined <- readRDS("output/full/lcc_low_combined.rds")

# Make model dataframe to store model names, metadata, etc.
model_df <- data.frame(
         model = c("orig_lcc", "springflow",
                  "inund_low", 
                  "inund_low_combined",
                  "lcc_low", 
                  "lcc_low_combined"),
         name = c("LCC IPM", "May-June IPM",
                   "Inundation IPM",
                   "Composite Inundation IPM",
                   "Hydraulic LCC IPM",
                   "Composite Hydraulic LCC IPM"),
         r2_sa = NA,
         r2_isl = NA,
         r2_ang = NA)

# Function to create and save output
output_fun <- function(mod, mod_name){
  print(mod_name)
  #Save model
  saveRDS(mod, paste0("output/full/", mod_name, ".rds"))
  
  #Samples (all iterations)
  #mcmc_samps <-  as.data.frame(mod)
  
  #Samples (just utilized iterations)
  mcmc_sub <-  as.data.frame(mod, pars = c("a", "mu_lbeta", "sd_lbeta", "B_lbeta", "effS"))
  
  #Summary stats
  mod_summ <- as.data.frame(summary(mod)$summary)
  
  #Save to csv
  #write_csv(mcmc_samps, "output/full/", mod_name, "_mcmc.csv")
  write_csv(mcmc_sub, paste0("output/", mod_name, "_mcmc_sub.csv"))
  write_csv(mod_summ, paste0("output/full/", mod_name, "_summ.csv"))
  
}

# Make model list to convert and save model output
model_list <- list(M2_1re, 
                   M2_2re, 
                   m_inund_i_l, 
                   m_inund_i_l_combined,
                   m_lcc_i_l, 
                   m_lcc_i_l_combined)

# Loop to run output function for all models in list
for (i in 1:length(model_list)) {
  mod_i <- model_list[[i]]
  mod_name_i <- model_df$model[i]
  
  output_fun(mod_i, mod_name_i)
  
}

#Save covariates so forecasting in other scripts can "center" predictor variables
saveRDS(inund_cov_list, "output/inund_cov_list.RData")
saveRDS(inund_cov_combined_list, "output/inund_cov_combined_list.RData")
saveRDS(inund_lcc_list, "output/inund_lcc_list.RData")
saveRDS(inund_lcc_list_combined, "output/inund_lcc_list_combined.RData")

# Model diagnostics #####

### R2 - All models ####

### r2 ####
calc_R2<-function(mod){
  eps1<-apply(rstan::extract(mod,"lbeta_eps")[[1]][,1,],1,var)
  lpred1<-apply(log(rstan::extract(mod,"Rf")[[1]][,,1]),1,var)
  eps2<-apply(rstan::extract(mod,"lbeta_eps")[[1]][,2,],1,var)
  lpred2<-apply(log(rstan::extract(mod,"Rf")[[1]][,,2]),1,var)
  eps3<-apply(rstan::extract(mod,"lbeta_eps")[[1]][,3,],1,var)
  lpred3<-apply(log(rstan::extract(mod,"Rf")[[1]][,,3]),1,var)
  r2_1<-1-eps1/lpred1
  r2_2<-1-eps2/lpred2
  r2_3<-1-eps3/lpred3
  return(c(mean(r2_1),mean(r2_2),mean(r2_3)))}

# Loop to run function for all models
for (i in 1:length(model_list)) {
  model_df[i,3:5] <- round(calc_R2(model_list[[i]]),2)
}

r2_table <- model_df %>%
  mutate(mean_r2 = (r2_sa+r2_isl+r2_ang)/3) %>%
  select("Model" = name,
         "San Acacia" = r2_sa, 
         "Isleta" = r2_isl,
         "Angostura" = r2_ang, 
         "Mean" = mean_r2)
  
knitr::kable(r2_table)

r2s_long <- r2_table %>%
  mutate(Model = factor(Model,
                        levels = rev(c("LCC IPM", 
                                       "May-June IPM", 
                                       "Inundation IPM", 
                                       "Composite Inundation IPM", 
                                       "Hydraulic LCC IPM", "Composite Hydraulic LCC IPM")))) %>%
  pivot_longer(cols = 2:5, names_to = "fit_type", values_to = "r2") %>%
  mutate(fit_type = factor(fit_type, levels = c("San Acacia","Isleta","Angostura", "Mean")),
         r2_label = sprintf("%.2f", r2))

# Code for Figure 3

sprintf("%02d", label_val)

ggplot(r2s_long, aes(fit_type, Model))+
  geom_tile(aes(fill = r2))+
  geom_hline(aes(yintercept = 4.5))+
  geom_text(aes(label = r2_label))+
  scale_fill_viridis_c()+
  labs(x = "Reach", fill = expression(R^2))+
  theme_ipsum()+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16, face = "plain"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

#ggsave("plots/r2_mods.jpeg", height = 2, width = 6, units = "in")

### Traceplots ####
# Check replace model name and parameters as needed
traceplot(m_inund_i_l_combined, pars = "B_lbeta", inc_warmup = TRUE)

### Model checks from Yackulic et al. 2022####
mod_check_fun <- function(model) {

lp0<-rstan::extract(model,"lpmon0")[[1]]
lp1<-rstan::extract(model,"lpmon1")[[1]]
lp01<-rstan::extract(model,"lpmon01")[[1]]
lpV<-rstan::extract(model,"lpmonV")[[1]]
lR0<-rstan::extract(model,"lpR0")[[1]]
lR1<-rstan::extract(model,"lpR1")[[1]]
sz<-rstan::extract(model,"sz")[[1]]
rsz<-rstan::extract(model,"rsz")[[1]]
#
pc0<-lp0
pc1<-lp1
pc01<-lp01
pcV<-lpV
pr0<-lR0
pr1<-lR1
#
q_pc0<-numeric()
q_pc1<-numeric()
q_pc01<-numeric()
q_pcV<-numeric()
q_pr0<-numeric()
q_pr1<-numeric()
qf<-function(obs,pred){
  temp<-which(obs==sort(pred))
  if (length(temp)>0) {sample(temp,1)/length(pred)} else {
    temp2<-findInterval(obs,sort(pred))
    if (temp2==0) {1/length(pred)} else { temp2/length(pred)}}}
#
for (i in 1:(dim(lp0)[2])){pc0[,i]<-rnbinom(15000,mu=exp(lp0[,i]),size=sz)
q_pc0[i]<-qf(mon0[i,5],pc0[,i])}
for (i in 1:(dim(lp1)[2])){pc1[,i]<-rnbinom(15000,mu=exp(lp1[,i]),size=sz)
q_pc1[i]<-qf(mon1[i,5],pc1[,i])}
for (i in 1:(dim(lp01)[2])){pc01[,i]<-rnbinom(15000,mu=exp(lp01[,i]),size=sz)
q_pc01[i]<-qf(mon01[i,5],pc01[,i])}
for (i in 1:(dim(lpV)[2])){pcV[,i]<-rnbinom(15000,mu=exp(lpV[,i]),size=sz)
q_pcV[i]<-qf(monV[i,5],pcV[,i])}
for (i in 1:(dim(lR0)[2])){pr0[,i]<-rnbinom(15000,mu=exp(lR0[,i]),size=rsz)
q_pr0[i]<-qf(R0[i,4],pr0[,i])}
for (i in 1:(dim(lR1)[2])){pr1[,i]<-rnbinom(15000,mu=exp(lR1[,i]),size=rsz)
q_pr1[i]<-qf(R1[i,4],pr1[,i])}
#combine monitoring data for plotting
q_mon<-c(q_pc0,q_pc1,q_pc01,q_pcV)
pc_mon<-cbind(pc0,pc1,pc01,pcV)
obsc_mon<-c(mon0[,5],mon1[,5],mon01[,5],monV[,5])
Q_mon<-c(mon0_effQ[,2],mon1_effQ[,2],mon01_effQ[,2],monV_effQ[,2])
#ditto for rescue
q_res<-c(q_pr0,q_pr1)
pc_res<-cbind(pr0,pr1)
obsc_res<-c(R0[,4],R1[,4])
jul_res<-c(R0[,2],R1[,2])
# code to make actual Fig S1 plot
par(mfrow=c(4,2))
par(mar=c(4,4,1,1))
xymax<-max(c(obsc_mon,apply(pc_mon,2,mean)))
plot(1+obsc_mon,1+apply(pc_mon,2,mean),ylim=c(1,1+xymax),xlim=c(1,1+xymax),axes=FALSE,xlab="",ylab="",pch=19,col=rgb(0,0,0,.1),main="Monitoring data",log="xy")
axis(1,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1)
axis(2,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1,las=T)
mtext("Observed catch",1,2,cex=0.7)
mtext("Predicted catch",2,2,cex=0.7)
text(1000,3,"R2 = 0.76")
text(.3,3000,"A)",xpd=T)
#
xymax<-max(c(obsc_res,apply(pc_res,2,mean)))
plot(1+obsc_res,1+apply(pc_res,2,mean),ylim=c(1,1+xymax),xlim=c(1,1+xymax),axes=FALSE,xlab="",ylab="",pch=19,col=rgb(0,0,0,.1),main="Rescue data",log="xy")
axis(1,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1)
axis(2,at=c(1,2,11,101,1001,10001),labels=c(0,1,10,100,1000,10000),pos=1,las=T)
mtext("Observed catch",1,2,cex=0.7)
mtext("Predicted catch",2,2,cex=0.7)
text(8000,3,"R2 = 0.77")
text(.15,30000,"B)",xpd=T)
# % zeros
iszero<-function(x){ifelse(x==0,1,0)}
hist(apply(iszero(pc_mon),1,mean),axes=FALSE,xlab="",ylab="",pch=19,main="")
abline(v=mean(iszero(obsc_mon)),col="red",lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=min(apply(iszero(pc_mon),1,mean)),at=c(0,750,1500,2250),labels=seq(0,0.15,.05),las=T)
mtext("Proportion of catch equal to zero",1,2,cex=0.7)
mtext("% of simulations",2,2,cex=0.7)
text(mean(iszero(obsc_mon))+.01,1500,"Observed",col="red")
text(0.52,2600,"C)",xpd=T)
# % zeros
hist(apply(iszero(pc_res),1,mean),axes=FALSE,xlab="",ylab="",pch=19,main="")
abline(v=mean(iszero(obsc_res)),col="red",lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=min(apply(iszero(pc_res),1,mean)),at=c(0,1500,3000,4500),labels=seq(0,0.3,.1),las=T)
mtext("Proportion of catch equal to zero",1,2,cex=0.7)
mtext("% of simulations",2,2,cex=0.7)
text(mean(iszero(obsc_res))-0.03,4000,"Observed",col="red")
text(0.14,4500,"D)",xpd=T)
#qq
plot(sort(q_mon),c(1:length(q_mon))/length(q_mon),pch=19,col=rgb(0,0,0,.01),main="",xlab="",ylab="",axes=FALSE)
curve(1*x,add=T,lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=0,las=T)
mtext("Observed quantile",1,2,cex=0.7)
mtext("Expected quantile",2,2,cex=0.7)
text(-0.16,1.1,"E)",xpd=T)
#qq
plot(sort(q_res),c(1:length(q_res))/length(q_res),pch=19,col=rgb(0,0,0,.1),main="",xlab="",ylab="",axes=FALSE)
curve(1*x,add=T,lty=2,lwd=2)
axis(1,pos=0)
axis(2,pos=0,las=T)
mtext("Observed quantile",1,2,cex=0.7)
mtext("Expected quantile",2,2,cex=0.7)
text(-0.16,1.1,"F)",xpd=T)
#
plot(Q_mon,q_mon,pch=19,col=rgb(0,0,0,.1),main="",xlab="",ylab="",axes=FALSE)
axis(1,at=c(0,0.5,1),labels=c(0,500,1000),pos=0)
axis(2,pos=0,las=T)
mtext("Discharge (cfs)",1,2,cex=0.7)
mtext("Observed quantile",2,2,cex=0.7)
text(-0.16,1.1,"G)",xpd=T)
#
plot(jul_res,q_res,pch=19,col=rgb(0,0,0,.1),main="",xlab="",ylab="",axes=FALSE)
axis(1,pos=0)
axis(2,pos=0,las=T)
mtext("Days after April 1st",1,2,cex=0.7)
mtext("Observed quantile",2,2,cex=0.7)
text(-27,1.1,"H)",xpd=T)

}

mod_check_fun(M2_2re)

# checking priors vs. posteriors
# code for figs S2 and S3 - can do for any fitted model with time variation in mortality rates - set mod = M2_1re for actual plots in appendix
checkpp_rep<-function(mod){
  par(mfrow=c(1,3))
  plot(density(rstan::extract(mod,"a")[[1]]),xlim=c(300,1500),axes=FALSE,main="",xlab="a",col="red")
  curve(dnorm(x,675,135),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=300,las=T)
  plot(density(rstan::extract(mod,"beta_stk")[[1]]),xlim=c(0.5,3),axes=FALSE,main="",xlab="beta_stk",col="red")
  curve(dnorm(x,1,0.1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0.5,las=T)
  plot(density(rstan::extract(mod,"beta_2")[[1]]),xlim=c(0.5,3),axes=FALSE,main="",xlab="beta_2",col="red")
  curve(dnorm(x,2,0.1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0.5,las=T)
} 

checkpp_rep(M2_2re)
checkpp_rep(m_inund_i_l)

checkpp_sd<-function(mod){
  par(mfrow=c(4,4))
  for (k in 1:3){
    plot(density(rstan::extract(mod,"mu_lM0")[[1]][,k]),xlim=c(-8,-3),axes=FALSE,main="",xlab=paste("mu_lM0",k),col="red")
    curve(dunif(x,-8,-3),col="blue",add=T,lty=2)
    axis(1,pos=0)
    axis(2,pos=-8,las=T)
  }
  for (k in 1:3){
    plot(density(rstan::extract(mod,"mu_lM1")[[1]][,k]),xlim=c(-8,-3),axes=FALSE,main="",xlab=paste("mu_lM1",k),col="red")
    curve(dunif(x,-8,-3),col="blue",add=T,lty=2)
    axis(1,pos=0)
    axis(2,pos=-8,las=T)
  }
  for (k in 1:3){
    plot(density(rstan::extract(mod,"mu_lMw")[[1]][,k]),xlim=c(-8,-3),axes=FALSE,main="",xlab=paste("mu_lMw",k),col="red")
    curve(dunif(x,-8,-3),col="blue",add=T,lty=2)
    axis(1,pos=0)
    axis(2,pos=-8,las=T)
  }
  plot(density(rstan::extract(mod,"irphi")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="irphi",col="red")
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0,las=T)
  plot(density(rstan::extract(mod,"p0")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="p0",col="red")
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0,las=T)
  plot(density(rstan::extract(mod,"p1")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="p1",col="red")
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0,las=T)
  plot(density(rstan::extract(mod,"rp0")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="rp0",col="red")
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0,las=T)
  plot(density(rstan::extract(mod,"rp1")[[1]]),xlim=c(0,1),axes=FALSE,main="",xlab="rp1",col="red")
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0,las=T)
  plot(density(rstan::extract(mod,"sz")[[1]]),main="",xlab="sz",col="red",axes=FALSE,xlim=c(0.25,0.6))
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0.25,las=T)
  plot(density(rstan::extract(mod,"rsz")[[1]]),main="",xlab="rsz",col="red",axes=FALSE,xlim=c(0.25,0.6))
  curve(dunif(x,0,1),col="blue",add=T,lty=2)
  axis(1,pos=0)
  axis(2,pos=0.25,las=T)
}

#Check different models here
checkpp_sd(M2_1re)

# Out-of-sample data ####

# Read in data
load(file = "output/oos_data_new.RData")

#Angostura reach
start.date <- "2021-01-01"
end.date <- "2024-12-31"
siteAng <- "08330000"
siteSan <- "08354900"
pCode <- "00060"

# Import new flow data
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

# Function to make out of sample forecasts
forecast_oos_re<-function(mod,Xout){
  Sp_N0<-rstan::extract(mod,"Sp_N")[[1]][,18,,]#last dimension is age / aug
  irphi<-rstan::extract(mod,"irphi")[[1]]
  beta_2<-rstan::extract(mod,"beta_2")[[1]]
  beta_stk<-rstan::extract(mod,"beta_stk")[[1]]
  a<-rstan::extract(mod,"a")[[1]]
  mu_lM0<-rstan::extract(mod,"mu_lM0")[[1]]
  mu_lM1<-rstan::extract(mod,"mu_lM1")[[1]]
  mu_lMw<-rstan::extract(mod,"mu_lMw")[[1]]
  sd_lM<-rstan::extract(mod,"sd_lM")[[1]]
  move<-rstan::extract(mod,"move")[[1]]
  rp0<-rstan::extract(mod,"rp0")[[1]]
  rp1<-rstan::extract(mod,"rp1")[[1]]
  mu_lbeta<-rstan::extract(mod,"mu_lbeta")[[1]]
  sd_lbeta<-rstan::extract(mod,"sd_lbeta")[[1]]
  B_lbeta<-rstan::extract(mod,"B_lbeta")[[1]]
  A0_perpool<-rstan::extract(mod,"A0_perpool")[[1]]
  AtQ_perpool<-rstan::extract(mod,"AtQ_perpool")[[1]]
  sl_width<-rstan::extract(mod,"sl_width")[[1]]
  bankfull<-rstan::extract(mod,"bankfull")[[1]]
  alpha0_max<-rstan::extract(mod,"alpha0_max")[[1]]
  alpha1_max<-rstan::extract(mod,"alpha1_max")[[1]]
  alpha0_int<-rstan::extract(mod,"alpha0_int")[[1]]
  alpha1_int<-rstan::extract(mod,"alpha1_int")[[1]]
  p0<-rstan::extract(mod,"p0")[[1]]
  p1<-rstan::extract(mod,"p1")[[1]]
  sz<-rstan::extract(mod,"sz")[[1]]
  sd_lbeta<-rstan::extract(mod,"sd_lbeta")[[1]]
  iter<-length(p0)
  #
  Mw<-array(NA,dim=c(iter,3,6))
  M0<-array(NA,dim=c(iter,3,6))
  M1<-array(NA,dim=c(iter,3,6))
  for (k in 1:Nstrata){
    Mw[,k,1]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
    Mw[,k,2]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
    Mw[,k,3]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
    Mw[,k,4]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
    Mw[,k,5]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
    Mw[,k,6]<-exp(rnorm(iter,mu_lMw[,k],sd_lM))
    
    M0[,k,1]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
    M0[,k,2]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
    M0[,k,3]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
    M0[,k,4]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
    M0[,k,5]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
    M0[,k,6]<-exp(rnorm(iter,mu_lM0[,k],sd_lM))
    
    M1[,k,1]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
    M1[,k,2]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
    M1[,k,3]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
    M1[,k,4]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
    M1[,k,5]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
    M1[,k,6]<-exp(rnorm(iter,mu_lM1[,k],sd_lM))
    
  }
  #
  FN<-array(NA,dim=c(iter,6,2,3))# iter, years,size classes, river segments		
  pC<-matrix(NA,ncol=Nobs_oosC,nrow=iter)
  C<-matrix(NA,ncol=Nobs_oosC,nrow=iter)
  SpN<-array(NA,dim=c(iter,6,3,3))# iter, years,size classes, river segments
  effS<-array(NA,dim=c(iter,6,3))# iter, years, river segments
  
  #add holder for number of age 0 recruits
  jul1_age0 <- array(NA, dim = c(iter, 6, 3))
  tRf_arr <- array(NA, dim = c(iter, 6, 3))

    for (k in 1:Nstrata){
    # Year 1
    SpN[,1,1,k]<-Sp_N0[,k,1] # Unique to year 1
    SpN[,1,2,k]<-Sp_N0[,k,2] # Unique to year 1
    SpN[,1,3,k]<-(w_m_oos[1,k,1]*exp(-60*Mw[,k,1])+w_m_oos[1,k,2]*exp(-151*Mw[,k,1]))*irphi
    effS[,1,k]=SpN[,1,1,k]+SpN[,1,2,k]*beta_2+SpN[,1,3,k]*beta_stk
    tN1<-rowSums(SpN[,1,1:2,k])
    tRf<-exp(rnorm(iter,mu_lbeta[,k],sd_lbeta)+B_lbeta*Xout[1,k])
    tRf_arr[,1,k] <- tRf # New holder to keep track of carrying capacity through time
    tN0=a*(effS[,1,k])/(1+a*(effS[,1,k])/tRf)		#this is the number of fish recruited based on covariates
    
    # new holder for # of recruits on Jul 1
    jul1_age0[,1,k] <- tN0
    
    FN[,1,1,k]= tN0*exp(-M0[,k,1]*124)*((1-cum_nd_oos[1,91,k])+(move-1)*(cum_nd_oos[1,215,k]-cum_nd_oos[1,91,k])+(1-move)*rp0*(cum_phiR_oos[1,215,k]-cum_phiR_oos[1,91,k]))
    FN[,1,2,k]=tN1*exp(-M1[,k,1]*215)*(1+(move-1)*cum_nd_oos[1,215,k]+(1-move)*rp1*cum_phiR_oos[1,215,k])			
    
    # Loop through years 2:6
    for (t in 2:6){
      SpN[,t,1:3,k]=cbind(FN[,t-1,1,k]*exp(-Mw[,k,t]*150),FN[,t-1,2,k]*exp(-Mw[,k,t]*150),(w_m_oos[t,k,1]*exp(-60*Mw[,k,t])+w_m_oos[t,k,2]*exp(-151*Mw[,k,t]))*irphi)
      effS[,t,k]=SpN[,t,1,k]+SpN[,t,2,k]*beta_2+SpN[,t,3,k]*beta_stk
      tN1<-rowSums(SpN[,t,1:2,k])
      tRf<-exp(rnorm(iter,mu_lbeta[,k],sd_lbeta)+B_lbeta*Xout[t,k])
      tRf_arr[,t,k] <- tRf
      tN0=a*(effS[,t,k])/(1+a*(effS[,t,k])/tRf)
      jul1_age0[,t,k] <- tN0
      FN[,t,1,k]= tN0*exp(-M0[,k,t]*124)*((1-cum_nd_oos[t,91,k])+(move-1)*(cum_nd_oos[t,215,k]-cum_nd_oos[t,91,k])+(1-move)*rp0*(cum_phiR_oos[t,215,k]-cum_phiR_oos[t,91,k]))
      FN[,t,2,k]=tN1*exp(-M1[,k,t]*215)*(1+(move-1)*cum_nd_oos[t,215,k]+(1-move)*rp1*cum_phiR_oos[t,215,k])
      
    }
    
  }			
  for (i in 1:Nobs_oosC){
    q<-oos_C$cQ[i]/1000
    totpool=exp(A0_perpool+AtQ_perpool*q)*200*sl_width[oos_C[i,5]]*q/(1+sl_width[oos_C[i,5]]*q/bankfull[oos_C[i,5]])
    totrun=(1-exp(A0_perpool+AtQ_perpool*q))*200*sl_width[oos_C[i,5]]*q/(1+sl_width[oos_C[i,5]]*q/bankfull[oos_C[i,5]])
    talpha0=c(alpha0_int*q/(1+alpha0_int*q/alpha0_max),1)
    talpha1=c(alpha1_int*q/(1+alpha1_int*q/alpha1_max),1)
    #y<-ifelse(oos_C[i,3]==2019,1,2)
    y<-oos_C[i,3] - 2018 # (i.e., 2019 will be y = 1)
    pC[,i]=(p0*oos_C$effort[i]*talpha0[oos_C$type[i]]*.2*FN[,y,1,oos_C[i,5]]*exp(-M0[,oos_C[i,5],(oos_C[i,3]-2018)]*(oos_C[i,2]-91))*
              ((1-cum_nd_oos[y,91,oos_C[i,5]])+(move-1)*(cum_nd_oos[y,oos_C[i,2],oos_C[i,5]]-cum_nd_oos[y,91,oos_C[i,5]])+
                 (1-move)*rp0*(cum_phiR_oos[y,oos_C[i,2],oos_C[i,5]]-cum_phiR_oos[y,91,oos_C[i,5]])))/
      (StrataLen_oos[y,oos_C[i,2],oos_C[i,5]]*(talpha0[1]*totpool+totrun))+
      (p1*oos_C$effort[i]*talpha1[oos_C$type[i]]*.2*FN[,y,2,oos_C[i,5]]*
         exp(-M1[,oos_C[i,5],(oos_C[i,3]-2018)]*oos_C[i,2])*(1+(move-1)*cum_nd_oos[y,oos_C[i,2],oos_C[i,5]]+(1-move)*rp1*cum_phiR_oos[y,oos_C[i,2],oos_C[i,5]]))/
      (StrataLen_oos[y,oos_C[i,2],oos_C[i,5]]*(talpha1[1]*totpool+totrun))
    C[,i]<-rnbinom(iter,mu=pC[,i],size=sz)
  }
  w19<-which(oos_C$year==2019)
  w20<-which(oos_C$year==2020)
  w21<-which(oos_C$year==2021)
  w22<-which(oos_C$year==2022)
  w23<-which(oos_C$year==2023)
  w24<-which(oos_C$year==2024)
  
  wSA19 <- which(oos_C$year==2019 & oos_C$Cstrata == 1)
  wIsl19 <- which(oos_C$year==2019 & oos_C$Cstrata == 2)
  wAng19 <- which(oos_C$year==2019 & oos_C$Cstrata == 3)
  wSA20 <- which(oos_C$year==2020 & oos_C$Cstrata == 1)
  wIsl20 <- which(oos_C$year==2020 & oos_C$Cstrata == 2)
  wAng20 <- which(oos_C$year==2020 & oos_C$Cstrata == 3)
  wSA21 <- which(oos_C$year==2021 & oos_C$Cstrata == 1)
  wIsl21 <- which(oos_C$year==2021 & oos_C$Cstrata == 2)
  wAng21 <- which(oos_C$year==2021 & oos_C$Cstrata == 3)
  wSA22 <- which(oos_C$year==2022 & oos_C$Cstrata == 1)
  wIsl22 <- which(oos_C$year==2022 & oos_C$Cstrata == 2)
  wAng22 <- which(oos_C$year==2022 & oos_C$Cstrata == 3)
  wSA23 <- which(oos_C$year==2023 & oos_C$Cstrata == 1)
  wIsl23 <- which(oos_C$year==2023 & oos_C$Cstrata == 2)
  wAng23 <- which(oos_C$year==2023 & oos_C$Cstrata == 3)
  wSA24 <- which(oos_C$year==2024 & oos_C$Cstrata == 1)
  wIsl24 <- which(oos_C$year==2024 & oos_C$Cstrata == 2)
  wAng24 <- which(oos_C$year==2024 & oos_C$Cstrata == 3)
  
  pc_6yr <- data.frame(pc19 = rowSums(pC[,w19]), pc20 = rowSums(pC[,w20]), pc21 = rowSums(pC[,w21]), 
                       pc22 = rowSums(pC[,w22]),pc23 = rowSums(pC[,w23]), pc24 = rowSums(pC[,w24]))
  
  ef_6yr <- c(sum(oos_C$effort[w19]), sum(oos_C$effort[w20]), sum(oos_C$effort[w21]), 
              sum(oos_C$effort[w22]), sum(oos_C$effort[w23]), sum(oos_C$effort[w24]))
  
  c_6yr <- data.frame(c19 = rowSums(C[,w19]), c20 = rowSums(C[,w20]), c21 = rowSums(C[,w21]), 
                      c22 = rowSums(C[,w22]),c23 = rowSums(C[,w23]), c24 = rowSums(C[,w24]))
  
  
  c_SA <- data.frame(c19 = rowSums(C[,wSA19]), c20 = rowSums(C[,wSA20]), c21 = rowSums(C[,wSA21]), 
                     c22 = rowSums(C[,wSA22]),c23 = rowSums(C[,wSA23]), c24 = rowSums(C[,wSA24]))
  
  c_Ang <- data.frame(c19 = rowSums(C[,wAng19]), c20 = rowSums(C[,wAng20]), c21 = rowSums(C[,wAng21]), 
                     c22 = rowSums(C[,wAng22]),c23 = rowSums(C[,wAng23]), c24 = rowSums(C[,wAng24]))
  
  c_Isl <- data.frame(c19 = rowSums(C[,wIsl19]), c20 = rowSums(C[,wIsl20]), c21 = rowSums(C[,wIsl21]), 
                     c22 = rowSums(C[,wIsl22]),c23 = rowSums(C[,wIsl23]), c24 = rowSums(C[,wIsl24]))
  
  return(list(ef_6yr=ef_6yr, pc_6yr=pc_6yr, c_6yr=c_6yr,SpN=SpN,FN=FN,effS=effS, jul1_age0 = jul1_age0, tRf_arr = tRf_arr,
              c_SA=c_SA, c_Isl=c_Isl, c_Ang =c_Ang))}

# Calculate flow data for out-of-sample forecasts
aQ_out <- angQ_all %>%
  filter(year %in% 2019:2024,
         month %in% 4:9) %>%
  select(year, cfs, month, day) %>%
  pivot_wider(names_from = year, values_from = cfs) %>%
  select(-month, -day) %>%
  as.matrix

sQ_out <- sanaQ_all %>%
  filter(year %in% 2019:2024,
         month %in% 4:9) %>%
  select(year, cfs, month, day) %>%
  pivot_wider(names_from = year, values_from = cfs) %>%
  select(-month, -day) %>%
  as.matrix

# Calcualte out-of-sample predictions for the original lcc model
preds_out<-array(NA,dim=c(6,3))
temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
t1<-numeric()
for (t in 1:6){
  se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
  t1[t]<-temp[t,6]-se
  se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
  t1[(t+7)]<-temp[(t+6),6]+se
}
t1[7]<-1
t1[14]<-0
temp<-subset(ee[[3]],ee[[3]][,1]=="2")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t2<-temp[1,6]
temp<-subset(ee[[3]],ee[[3]][,1]=="4")
se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
t4<-temp[1,6]-se
for (r in 1:3){
  temp<-subset(ee[[3]],ee[[3]][,1]==riversegment[r]&is.na(ee[[3]][,5])==FALSE)
  t3<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t3[t]<-temp[t,6]+se
  }	
  for (t in 8:15){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t3[t]<-temp[t,6]-se
  }
  t3[7]<-temp[7,6]
  for (j in 1:6){
    if (r==3) {q<-aQ_out[,j]} else {q<-sQ_out[,j]}
    preds_out[j,r]<-calc_cov(q,t3,t2,t4,t1) # only uses present habitat
  }}

#### out of sample forecast from integrated model with larval carrying capacity
f1<-forecast_oos_re(M2_1re,preds_out-mean(preds2))

# out of sample forecasts from models based on may-june average flow
mjflow_sana_out<-subset(sanaQ_all,(sanaQ_all$month==5|sanaQ_all$month==6)&sanaQ_all$year>2018)
mjflow_ang_out<-subset(angQ_all,(angQ_all$month==5|angQ_all$month==6)&angQ_all$year>2018)
simpflow_out<-cbind(tapply(mjflow_sana_out$cfs,mjflow_sana_out$year,mean),tapply(mjflow_sana_out$cfs,mjflow_sana_out$year,mean),tapply(mjflow_ang_out$cfs,mjflow_ang_out$year,mean))/1000
f2<-forecast_oos_re(M2_2re,simpflow_out-mean(simpflow))

#### out of sample forecast model using raw inundation
# Calculates 4 sets of covariates but we only use the first set (raw inundation)
inund_list_out <- list()

for (i in 1:4) {
  maxmins <- array(NA, dim = c(6,3))
  
  for (j in 1:3) {
    if (j == 1) {
      q <- filter(sanaQ_all, year %in% 2019:2024)
    } else {
      q <- filter(angQ_all, year %in% 2019:2024)
    }
    
    hab_lookup_ij <- filter(q_hab_lookup, 
                            reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup)[i+2])
    
    maxmins[,j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
    
  }
  inund_list_out[[i]] <- maxmins
}

f3<-forecast_oos_re(m_inund_i_l, log(inund_list_out[[1]]+1)-mean(log(inund_cov_list[[1]]+1)))

#### out of sample forecast from model using combined inundation and elicited hab. values
# Calculates 4 sets of covariates but we only use the first set (raw inundation hydraulic output)

inund_combined_list_out <- list()

for (i in 1:4) {
  maxmins <- array(NA, dim = c(6,3))
  
  for (j in 1:3) {
    if (j == 1) {
      q <- filter(sanaQ_all, year %in% 2019:2024)
    } else {
      q <- filter(angQ_all, year %in% 2019:2024)
    }
    
    hab_lookup_ij <- filter(q_hab_lookup_combined, 
                            reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup_combined)[i+2])
    
    maxmins[,j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
    
  }
  inund_combined_list_out[[i]] <- maxmins
}

f4<-forecast_oos_re(m_inund_i_l_combined, log(inund_combined_list_out[[1]]+1)-mean(log(inund_cov_combined_list[[1]]+1)))

# Model that incorporates inundation into larval carrying capacity

# calculate new set of larval carry capacity indices based on inundation

lcc_list_out <- list()

for (i in 1:4){
  
  preds<-array(NA,dim=c(6,3))
  
  temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
  t1<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t1[t]<-temp[t,6]-se
    se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
    t1[(t+7)]<-temp[(t+6),6]+se
  }
  t1[7]<-1
  t1[14]<-0
  temp<-subset(ee[[3]],ee[[3]][,1]=="2")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t2<-temp[1,6]
  temp<-subset(ee[[3]],ee[[3]][,1]=="4")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t4<-temp[1,6]-se
  #
  for (r in 1:3){
    t3 <- filter(q_inund_curve, reach_num == r) %>%
      #mutate(inund_acres = log(inund_acres+1)) %>%
      pull(i+1)
    for (t in 1:6){
      if (r==3) {q<-aQ_out[,t]} else {q<-sQ_out[,t]}
      preds[t,r]<-calc_cov_2d(q,t3,t2,t4,t1)
    }}
  
  lcc_list_out[[i]] <- preds
}

f5<-forecast_oos_re(m_lcc_i_l,lcc_list_out[[1]]-mean(inund_lcc_list[[1]]))

# calculate new set of larval carry capacity indices based on inund/expert elicitation

lcc_list_out_combined <- list()

for (i in 1:4){
  
  preds<-array(NA,dim=c(6,3))
  
  temp<-subset(ee[[3]],ee[[3]][,1]=="1"&is.na(ee[[3]][,5])==FALSE)
  t1<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t1[t]<-temp[t,6]-se
    se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
    t1[(t+7)]<-temp[(t+6),6]+se
  }
  t1[7]<-1
  t1[14]<-0
  temp<-subset(ee[[3]],ee[[3]][,1]=="2")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t2<-temp[1,6]
  temp<-subset(ee[[3]],ee[[3]][,1]=="4")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t4<-temp[1,6]-se
  #
  for (r in 1:3){
    t3 <- filter(q_inund_combined_format, reach_num == r) %>%
      #mutate(inund_acres = log(inund_acres+1)) %>%
      pull(i+1)
    for (t in 1:6){
      if (r==3) {q<-aQ_out[,t]} else {q<-sQ_out[,t]}
      preds[t,r]<-calc_cov_combined(q,t3,t2,t4,t1)
    }}
  
  lcc_list_out_combined[[i]] <- preds
}

f6<-forecast_oos_re(m_lcc_i_l_combined,lcc_list_out_combined[[1]]-mean(inund_lcc_list_combined[[1]]))

# calculate observation
obs_catch<-tapply(oos_C$hybama,oos_C$year,sum)

# calculate prediction using biop
biop_covs <- angQ_all %>%
  filter(month %in% 5:6, year %in% 2019:2024) %>%
  group_by(year) %>%
  summarize(tbQ = mean(cfs)) %>%
  mutate(effort = f1$ef_6yr)

biop_preds <- matrix(NA, nrow = 1000000, ncol = 6)

for (i in 1:nrow(biop_covs)){
  biop_preds[,i] <- (10^(rnorm(1000000,-0.1477 + (0.0004*biop_covs$tbQ[i])-(biop_covs$tbQ[i]^2) * 0.000000014284,.1)) - 1)*biop_covs$effort[i]/100
  
}

## Full subset of models ####

# out of sample figure
q025<-function(x){quantile(x,.025)}
q975<-function(x){quantile(x,.975)}
q10<-function(x){quantile(x,.1)}
q90<-function(x){quantile(x,.9)}

## Final subset of models #### 

mod_list_oos <- list("LCC IPM"= f1$c_6yr,
                     "May-June IPM" = f2$c_6yr,
                     "Inund. IPM" = f3$c_6yr,
                     "Compos. Inund. IPM" = f4$c_6yr,
                     "Hydraul. LCC" = f5$c_6yr,
                     "Compos. Hydraul. LCC" = f6$c_6yr,
                                     "BioOp Lin. Mod." = as.data.frame(biop_preds))

mod_out_oos <- data.frame()

for (i in 1:length(mod_list_oos)){
  
  output <- mod_list_oos[[i]]
  
  df <- data.frame(model = rep(names(mod_list_oos[i]), 6),
                   year = 2019:2024,
                   mean = colMeans(output),
                   median = sapply(output, median),
                   min = sapply(output, min),
                   l025 = sapply(output, q025),
                   l10 = sapply(output, q10),
                   u90 = sapply(output, q90),
                   u975 = sapply(output, q975))

  mod_out_oos <- bind_rows(mod_out_oos, df)
}

obs_catch_df <- data.frame(year = 2019:2024, 
                           catch = obs_catch)

mod_out_catch <- left_join(mod_out_oos, obs_catch_df) %>%
  mutate(model = factor(model,
                           levels = c("LCC IPM", "May-June IPM", "BioOp Lin. Mod.",
                                      "Hydraul. LCC", "Compos. Hydraul. LCC", 
                                      "Inund. IPM", "Compos. Inund. IPM")))

mod_out_catch %>%
  ggplot(aes(x = model))+
  geom_point(aes(y = mean), color = "black")+
  geom_errorbar(aes(ymin = l025,ymax = u975), 
               width = 0, color="red")+
  geom_errorbar(aes(ymin = l10,ymax = u90), 
                width = 0, color="black")+
  geom_hline(data = obs_catch_df, 
             aes(yintercept = catch), 
             lty = 2)+
  geom_rect(aes(xmin=3.5, xmax=7.5, ymin=0, ymax=Inf), 
            fill = "lightblue", alpha = .1)+
  facet_wrap(~year, scales = "free_y", nrow = 2)+
  labs(x = "Model", y = "RGSM Catch")+
  theme_ipsum(axis_title_just = 1,
              base_size = 14,
              strip_text_size = 18,
              axis_title_size = 16) +
  theme(axis.text.x = element_text(angle = 90),
        panel.spacing.y = unit(0.1, "lines"))

ggsave(filename = "plots/fig3_oos_results_new_mod_order.jpeg", bg = "white", width = 8, height = 7, units = "in")
