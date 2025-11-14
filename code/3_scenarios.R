library(tidyverse) 
library(zoo) # Rolling mins.
library(hrbrthemes) # Plotting
library(sf) # Spatial analysis
library(patchwork) # Composite figures
library(ggthemes) # For scale_color_colorblind()

# ## Helper functions ----
# Calculate significance level for expert informed covariate
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

# Calculate quantiles
q025 <- function(x) {quantile(x, .025)}
q975 <- function(x) {quantile(x, .975)}
q10 <- function(x) {quantile(x, .1)}
q90 <- function(x) {quantile(x, .9)}

## Selecting hydrographs for simulation ----
q_sim <- function(flow_input) {
  q_all <- array(NA, dim = c(214, length(flow_input), 2))
  
  if (length(flow_input[flow_input < 1]) > 0) {
    for (i in 1:length(flow_input)) {
      #Get year that most closely matches that flow
      aQ_quant <- quantile(aQ_med_spr$tot_vol, flow_input[i])
      a_year_sel <- aQ_med_spr$year[which.min(abs(aQ_quant - aQ_med_spr$tot_vol))]
      
      sQ_quant <- quantile(sQ_med_spr$tot_vol, flow_input[i])
      s_year_sel <- sQ_med_spr$year[which.min(abs(sQ_quant - sQ_med_spr$tot_vol))]
      
      #Plot discharge from select year vs. median daily flow
      aQ_temp <- filter(angQ, year == a_year_sel, month > 2, month < 10)
      
      sQ_temp <- filter(sanaQ, year == s_year_sel, month > 2, month < 10)
      
      q_all[, i, 1] <- aQ_temp$cfs
      q_all[, i, 2] <- sQ_temp$cfs
    }
  } else {
    for (i in 1:length(flow_input)) {
      #Plot discharge from select year vs. median daily flow
      aQ_temp <- filter(angQ, year == flow_input[i], month > 2, month < 10)
      
      sQ_temp <- filter(sanaQ, year == flow_input[i], month > 2, month < 10)
      
      q_all[, i, 1] <- aQ_temp$cfs
      q_all[, i, 2] <- sQ_temp$cfs
      
    }
  }
  return(q_all)
}

## Functions for calculating the larval carrying capacity index covariate ----

# Calculates larval habitat based on Q value
#This was altered to match hydraulic breakpoints
q2hab_2d <- function(q, pars, input) {
  p <- c(pars, pars[length(pars)])
  qs <- c(input, 2 * 10 ^ 4)
  t1 <- findInterval(q, qs)
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}

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

calc_cov_combined <- function(Q, hab_lookup_reach, kappa, D, prPARS) {
  thab <- numeric()
  for (i in 1:214) {
    thab[i] <- log(hab_lookup_reach$hab[hab_lookup_reach$cfs == round(Q[i])] + 1)
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

preds_calc <- function(exp, flows, hab_lookup) {
  preds <- array(NA, dim = c(dim(flows)[2], 3))
  
  temp<-subset(ee[[exp]],ee[[exp]][,1]=="1"&is.na(ee[[exp]][,5])==FALSE)
  t1<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t1[t]<-temp[t,6]-se
    se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
    t1[(t+7)]<-temp[(t+6),6]+se
  }
  t1[7]<-1
  t1[14]<-0
  temp<-subset(ee[[exp]],ee[[exp]][,1]=="2")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t2<-temp[1,6]
  temp<-subset(ee[[exp]],ee[[exp]][,1]=="4")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t4<-temp[1,6]-se
  #
  for (r in 1:3){
    t3 <- filter(hab_lookup, reach_num == r) #%>%
    for (j in 1:dim(flows)[2]) {
      if (r==3) {q<-flows[, j, 1]} else {q<-flows[,j,2]}
      preds[j,r]<-calc_cov_combined(q,t3,t2,t4,t1)
    }}
  
  return(preds)
  
}

## Functions for calculating inundation covariate ----

# Q-habitat lookup function
hab_lookup_fun <- function(inund_curve, reaches) {
  hab_lookup <- tibble()
  
  for (i in 1:length(reaches)) {
    q_hab_i <- filter(inund_curve, reach == reaches[i])
    
    sites <- unique(q_hab_i$site)
    
    hab_lookup_i <- tibble()
    
    for (j in 1:length(sites)) {
      q_hab_ij <- if ("site" %in% colnames(q_hab_i)) {
        q_hab_i %>% filter(site %in% sites[j])
      } else {
        q_hab_i  # Return the original dataframe unfiltered if the column doesn't exist
      }
      
      hab_lookup_ij <- tibble(cfs = 1:10e3) %>%
        mutate(hab = q2hab_2d(cfs, q_hab_ij$hab, q_hab_ij$cfs))
      
      hab_lookup_i <- bind_rows(hab_lookup_i, hab_lookup_ij)
      
    }
    
    hab_lookup_i <- group_by(hab_lookup_i, cfs) %>% summarize(hab = sum(hab)) %>%
      mutate(reach = reaches[i])
    
    hab_lookup <- bind_rows(hab_lookup, hab_lookup_i)
  }
  return(hab_lookup)
  
}

# Function for calculating inundation covariate

maxmins_fun <- function(lookup, sim_q, reaches) {
  
  maxmins_hab <- array(NA, dim = c(dim(sim_q)[2], 3))
  maxmins_cfs <- array(NA, dim = c(dim(sim_q)[2], 3))
  
  for (i in 1:length(reaches)){
    hab_lookup_i <- filter(lookup, reach == reaches[i])
    
    if (i == 1) {
      q <- round(sim_q[, , 2], 0) %>% as_tibble() %>% set_names("cfs")
    } else {
      q <- round(sim_q[, , 1], 0) %>% as_tibble() %>% set_names("cfs")
    }
    
    for (t in 1:dim(sim_q)[2]) {
      x <-  q[, t]
      
      maxmins_hab[t, i] <- max(rollapplyr(hab_lookup_i$hab[match(x$cfs, hab_lookup_i$cfs)],
                                          Qdur, mean, fill = NA)[32:122], na.rm = T)
      maxmins_cfs[t,i] <- max(rollapplyr(x$cfs,
                                         Qdur, mean, fill = NA)[32:122], na.rm = T)
      
    }
  }
  
  return(list(maxmins_hab = maxmins_hab, maxmins_cfs = maxmins_cfs))
  
}

## Flow augmentation function ----

aug_q_fun <- function(st_date, en_date, vol_AcFt, q_array) {
  simq_dates = data.frame(date = seq(as.Date("2000-03-01"), as.Date("2000-09-30"), 1),
                          day_ind = 1:214) %>%
    mutate(date = as.character(format(date, "%m-%d")))
  
  start_ind <- simq_dates$day_ind[simq_dates$date == st_date]
  end_ind <- simq_dates$day_ind[simq_dates$date == en_date]
  
  aug_vec <- rep(0, 214)
  
  ndays_aug <- length(start_ind:end_ind)
  
  augCFS = round(vol_AcFt * 43559.935 / ndays_aug / (24 * 60 * 60), 0)
  
  aug_vec[start_ind:end_ind] <- rep(augCFS, ndays_aug)
  
  aug_q_array <- q_array + aug_vec
  
  return(aug_q_array)
}


## Forecast recruit function ----
#Input: Model of choice, Covariates on Recruitment
forecast_recruit <- function(mod, 
                             Xout, 
                             years_sim,  
                             maxminsBase, 
                             maxminsAug, 
                             maxminsRest, 
                             reaches) {
  
  #Extract model parameters
  a <- pull(mod, a)
  mu_lbeta <- select(mod, starts_with("mu_lbeta")) %>% as.matrix()
  sd_lbeta <- pull(mod, sd_lbeta)
  B_lbeta <- pull(mod, B_lbeta)
  effS_all <- select(mod, starts_with("effS"))
  iter <- length(a)
  
  #Extract effective spawner numbers for all years and select desired level of effS rank

  #Create long version
  effS_long <- effS_all %>%
    mutate(sample = 1:15000) %>%
    gather(-sample, key = "param", value = "effS") %>%
    mutate(year = rep(rep(c(1:17), each = iter), 3),
           reach = rep(c(1:3), each = iter * 17))
  
  effS_medians <- effS_long %>%
    group_by(year, reach) %>%
    summarize(median_effS = median(effS)) %>%
    spread(key = reach, value = median_effS) %>%
    ungroup() %>%
    select(-year)
  
  effS_rank <- apply(effS_medians, 2, rank)
  effS_ind <- which(effS_rank == spawn_rank, arr.ind = TRUE)[, 1]
  
  #Holder for recruitment predictions
  age0_base <- array(NA, dim = c(iter, nrow(Xout), 3))
  age0_rest <- array(NA, dim = c(iter, nrow(Xout), 3))
  age0_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  diff_rest <- array(NA, dim = c(iter, nrow(Xout), 3))
  prop_rest <- array(NA, dim = c(iter, nrow(Xout), 3))
  diff_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  prop_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  ratio_rest <- array(NA, dim = c(iter, nrow(Xout), 3))
  ratio_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  effS_sim <- array(NA, dim = c(iter, 3))
  
  for (k in 1:3) {
    # Loop through 3 reaches
    
    #Selects year with desired spawner characteristics
    effS <- filter(effS_long, year == effS_ind[k], reach == k) %>%
      pull(effS)
    
    for (i in 1:nrow(Xout)) {
      tRf0 <- exp(rnorm(iter, mu_lbeta[, k], sd_lbeta)) # mu_lbeta = river specific intercept
      
      tRf_base <- exp(B_lbeta * Xout[i, k, 1])*tRf0 #tRf = Larval carrying capacity; 
      tRf_rest <- exp(B_lbeta * Xout[i, k, 2])*tRf0
      tRf_hydro <- exp(B_lbeta * Xout[i, k, 3])*tRf0
      
      age0_base[, i, k] = a * (effS) / (1 + a * (effS) / tRf_base)		#this is the baseline number of fish recruited plus those as function of inundation
      age0_rest[, i, k] = a * (effS) / (1 + a * (effS) / tRf_rest)
      age0_hydro[, i, k] = a * (effS) / (1 + a * (effS) / tRf_hydro)
      diff_rest[,i,k] = age0_rest[, i, k] - age0_base[, i, k]
      diff_hydro[,i,k] = age0_hydro[, i, k] - age0_base[, i, k]
      prop_hydro[, i, k] = diff_hydro[, i, k]/age0_base[, i, k]
      prop_rest[, i, k] = diff_rest[, i, k]/age0_base[, i, k]
      ratio_hydro[, i, k] = age0_hydro[, i, k]/age0_base[, i, k]
      ratio_rest[, i, k] = age0_rest[, i, k]/age0_base[, i, k]
    }
    
    effS_sim[, k] <- effS
  }
  
  # Total recruits w/ and w/o scenarios
  age0_baseline <- apply(age0_base, c(2, 3), median) %>% as.vector()
  lower_age0 <- apply(age0_base, c(2, 3), q10) %>% as.vector()
  upper_age0 <- apply(age0_base, c(2, 3), q90) %>% as.vector()
  
  age0_rest_med <- apply(age0_rest, c(2, 3), median) %>% as.vector()
  lower_age0_rest <- apply(age0_rest, c(2, 3), q10) %>% as.vector()
  upper_age0_rest <- apply(age0_rest, c(2, 3), q90) %>% as.vector()
  
  age0_hydro_med <- apply(age0_hydro, c(2, 3), median) %>% as.vector()
  lower_age0_hydro <- apply(age0_hydro, c(2, 3), q10) %>% as.vector()
  upper_age0_hydro <- apply(age0_hydro, c(2, 3), q90) %>% as.vector()
  
  # Raw difference in recruits 
  diff_rest_med <- apply(diff_rest, c(2, 3), median) %>% as.vector()
  lower_diff_rest <- apply(diff_rest, c(2, 3), q10) %>% as.vector()
  upper_diff_rest <- apply(diff_rest, c(2, 3), q90) %>% as.vector()
  
  diff_hydro_med <- apply(diff_hydro, c(2, 3), median) %>% as.vector()
  lower_diff_hydro <- apply(diff_hydro, c(2, 3), q10) %>% as.vector()
  upper_diff_hydro <- apply(diff_hydro, c(2, 3), q90) %>% as.vector()
  
  # Diff. in recruits over baseline recruits
  prop_rest_med <- apply(prop_rest, c(2, 3), median) %>% as.vector()
  lower_prop_rest <- apply(prop_rest, c(2, 3), q10) %>% as.vector()
  upper_prop_rest <- apply(prop_rest, c(2, 3), q90) %>% as.vector()
  
  prop_hydro_med <- apply(prop_hydro, c(2, 3), median) %>% as.vector()
  lower_prop_hydro <- apply(prop_hydro, c(2, 3), q10) %>% as.vector()
  upper_prop_hydro <- apply(prop_hydro, c(2, 3), q90) %>% as.vector()
  
  # Ratio of recruits w/ and w/o scenarios
  ratio_rest_med <- apply(ratio_rest, c(2, 3), median) %>% as.vector()
  lower_ratio_rest <- apply(ratio_rest, c(2, 3), q10) %>% as.vector()
  upper_ratio_rest <- apply(ratio_rest, c(2, 3), q90) %>% as.vector()
  
  ratio_hydro_med <- apply(ratio_hydro, c(2, 3), median) %>% as.vector()
  lower_ratio_hydro <- apply(ratio_hydro, c(2, 3), q10) %>% as.vector()
  upper_ratio_hydro <- apply(ratio_hydro, c(2, 3), q90) %>% as.vector()
  
  # Compile results dataframe
  age0_dat <- tibble(
    reach = rep(reaches, each = nrow(Xout)),
    flow_years = rep(years_sim, 3), #toggle this line off for Kyle plots
    cfs_baseline = as.vector(maxminsBase$maxmins_cfs),
    cfs_augmentation = as.vector(maxminsAug$maxmins_cfs),
    acres_inund_baseline = as.vector(maxminsBase$maxmins_hab),
    acres_inund_rest =  as.vector(maxminsRest$maxmins_hab),
    acres_inund_aug =  as.vector(maxminsAug$maxmins_hab),
    recruits_baseline = age0_baseline,
    lower_age0 = lower_age0,
    upper_age0 = upper_age0,
    
    recruits_rest = age0_rest_med,
    lower_age0_rest = lower_age0_rest,
    upper_age0_rest = upper_age0_rest,
    
    recruits_hydro = age0_hydro_med,
    lower_age0_hydro = lower_age0_hydro,
    upper_age0_hydro = upper_age0_hydro,
    
    diff_recruits_rest = diff_rest_med,
    lower_diff_rest = lower_diff_rest,
    upper_diff_rest = upper_diff_rest,
    
    diff_recruits_hydro = diff_hydro_med,
    lower_diff_hydro = lower_diff_hydro,
    upper_diff_hydro = upper_diff_hydro,
    
    prop_rest_med = prop_rest_med,
    lower_prop_rest = lower_prop_rest,
    upper_prop_rest = upper_prop_rest,
    
    prop_hydro_med = prop_hydro_med,
    lower_prop_hydro = lower_prop_hydro,
    upper_prop_hydro = upper_prop_hydro,
    
    ratio_rest_med = ratio_rest_med,
    lower_ratio_rest = lower_ratio_rest,
    upper_ratio_rest = upper_ratio_rest,
    
    ratio_hydro_med = ratio_hydro_med,
    lower_ratio_hydro = lower_ratio_hydro,
    upper_ratio_hydro = upper_ratio_hydro
  )
  
  return(age0_dat)
  
}

# Import input data
reaches <-c("San Acacia", "Isleta", "Angostura")

# Read in USGS gage data
angQ <- read.csv("data/yackulic2022_data/abq_gage_08330000.csv") %>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))
sanaQ <- read.csv("data/yackulic2022_data/SanAcacia_gage_08354900.csv")%>% 
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

## Prepare flow data ----

#Get median daily flows from 2 gages (March - September) for plotting
aQ_med_dly <- angQ %>%
  filter(month > 2, month < 10) %>%
  group_by(month, day) %>%
  summarise(med_cfs = median(cfs))

sQ_med_dly <- sanaQ %>%
  filter(month > 2, month < 10) %>%
  group_by(month, day) %>%
  summarise(med_cfs = median(cfs))

#Get median Apr-Jun flows by year to rank wet vs. dry years
aQ_med_spr <- angQ %>%
  mutate(vol = cfs*86400) %>%
  filter(month > 3, month < 7) %>%
  group_by(year) %>%
  summarize(med_cfs = median(cfs),
            tot_vol = sum(vol))

sQ_med_spr <- sanaQ %>%
  mutate(vol = cfs*86400) %>%
  filter(month > 3, month < 7) %>%
  group_by(year) %>%
  summarize(med_cfs = median(cfs),
            tot_vol = sum(vol))

# Import model covariates (so prediction variables can be centered)
inund_lcc_list_combined <- readRDS("output/inund_lcc_list_combined.RData")

# Import and re-arrange baseline flow-inundation curves
q_inund_curve <- read_csv("data/ecoval2d.csv") %>% # 2D
  select(cfs = q, 
         inund_lower = tot_outside_chan_lower,
         inund_upper = tot_outside_chan_upper,
         wua_lower, wua_upper,
         reach, reach_num) %>%
  pivot_longer(cols = inund_lower:wua_upper,
               names_to = "hab_type",
               values_to = "hab")

### Code for Figure 2 ####

#https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
pal_df <- data.frame(reach = reaches,
                     color = c("#7570b3","#1b9e77","#d95f02"),
                     ltype = c(4,1,3))

cols <- setNames(pal_df$color, pal_df$reach)
ltype <- setNames(pal_df$ltype, pal_df$reach)

# Plot q-hab curves
# Converting to metric for pub.
q_inund_curve %>%
  filter(hab_type == "inund_lower") %>%
  ggplot()+
  geom_line(aes(cfs*0.0283168466, 
                hab/247.1, color = reach, linetype = reach), size = 1.5)+
  scale_color_manual(values = cols)+
  scale_linetype_manual(values = ltype)+
  labs(x=expression("Discharge (m"^{3}*"/s)"), y = expression("Floodplain inundation (km"^{2}*")"))+
  theme_ipsum()+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "lines"))


#ggsave(filename = "plots/ecovalue_curve.jpeg", width = 8, height = 6, units = "in")

# Import habitat lookup from combined hydraulic/expert elicitation
q_hab_lookup_combined <- read_csv("data/q_hab_lookup_combined.csv")

#Select desired spawner rank (1 to 17, with 1 = minimum, 9 = median, 17 = max)
spawn_rank <- 9

# Select desired mod
# This selects what type of covariate is used for recruitment (cc, inundation, tihm)
mod_name <- "lcc_low_combined"

# Select inundation covariate type; must correspond with model that was chosen
hab_cov <- "inund_lower"

# Select flow duration if using inundation covariate; this is the value from EE 3
Qdur <- 20

# Load select fitted model -------------------------------------------------------------
sim_mod <- read_csv(paste0("output/",
                           mod_name, 
                           "_mcmc_sub.csv"))

# Select cov. list that corresponds to model for standardizing; index value for picking l, h, wua-l, wua-h
cov_list <- inund_lcc_list_combined[[1]]

# Reach num dataframe; reach numbers needed for functins
reach_df <- data.frame(reach = reaches, 
                       reach_num = c(1,2,3))

# No flow habitat df; need to add habitat at 0 cfs
no_flow_hab <- reach_df %>%
  mutate(hab = 0, cfs = 0)

# Baseline flow-habitat lookup table, no with reach numbers and 0 cfs habitat
hab_lookup_baseline <- q_hab_lookup_combined %>%
  select(cfs, reach, hab = any_of(hab_cov)) %>%
  left_join(reach_df) %>%
  bind_rows(no_flow_hab) %>%
  arrange(reach, cfs)

# Restoration scenario #####

# Import raw restoration hydraulic data
q_hab_rest_raw <- read_csv("data/san_acacia_rest_ecovals.csv") %>%
  filter(terrain == "AB") %>%
  select(-terr_inund_notes, raw_hab = hab)

# Code for Figure S1 
ggplot(q_hab_rest_raw) +
  geom_line(aes(cfs*0.0283168466, raw_hab*4046.86/1000, color = as.factor(site)))+
  scale_color_colorblind()+
  labs(x=expression("Discharge (m"^{3}*"/s)"), y = expression("Floodplain inundation (in thousands of m"^{2}*")"))+
  theme_ipsum()+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "lines"))

#ggsave(filename = "plots/rest_q_hab.jpeg", width = 8, height = 6, units = "in")

# How much inundated habitat is in San Acacia at 7000 cfs
q_inund_curve %>%
  filter(hab_type == "inund_lower", 
         reach == "San Acacia",
         cfs == 7000)

# Normalize by area of reach at 7,000 cfs: create ratio of habitat at restoration site @7000 to habitat at reach-scale at   
sana_site_ratio <- q_hab_rest_raw %>%
  filter(cfs==7000) %>%
  select(site, size, raw_hab) %>%
  mutate(ratio = raw_hab/6284)

# Q-hab curve for the non restoration reaches; required for hab_lookup function
no_rest_reaches = reaches[c(2,3)]

# Need to add q-hab for Angostura and Isleta; recruit functions have all three reaches
q_hab_curves_non_rest <- tibble(cfs = rep(c(0,7000), length(no_rest_reaches)),
                                reach = rep(no_rest_reaches, each = 2),
                                hab = 0)

# Create dataframe of all lookup values
hab_lookup_rest_all <- tibble()
rest_sites <- unique(q_hab_rest_raw$site)

# Initial inundation at each reach; needed if merging effect
init_inund <- q_hab_rest_raw %>% 
  select(site, cfs, raw_hab) %>%
  filter(cfs != 0) %>%
  group_by(site) %>%
  slice_min(raw_hab) %>%
  slice_max(cfs) %>%
  select(site, cfs)

restReach <- "San Acacia"

for (i in 1:length(rest_sites)){
  #Extract reach prop
  reach_prop_i <- filter(sana_site_ratio, site == rest_sites[i]) %>% pull(ratio)
  
  init_inund_i <- filter(init_inund, site == rest_sites[i]) %>% pull()
  
  #Baseline habitat scaled to size of reach
  baseline_hab_prop_i <- hab_lookup_baseline %>%
    mutate(hab_site_base = hab*reach_prop_i)
  
  # Hab @ start of inund.; add this above initial inund
  hab_channel_max <- baseline_hab_prop_i %>% filter(reach == restReach,
                                 cfs == init_inund_i) %>%
    pull(hab_site_base)
  
  q_hab_i <- filter(q_hab_rest_raw, site == rest_sites[i]) %>%
    select(cfs, reach, hab = raw_hab, site)
  
  curves_i <- bind_rows(q_hab_curves_non_rest, q_hab_i)
  
  # Adding channel habitat
  lookup_i <- hab_lookup_fun(curves_i, reaches) %>% mutate(site = rest_sites[i]) %>%
    rename(raw_hab = hab) %>%
    left_join(baseline_hab_prop_i) %>%
    mutate(rest_hab = case_when(cfs > init_inund_i ~ raw_hab + hab_channel_max,
                                TRUE ~ hab_site_base)) %>%
    mutate(rest_hab_norm = rest_hab - hab_site_base)

  # Not adding channel habitat
  lookup_i_2 <- hab_lookup_fun(curves_i, reaches) %>% mutate(site = rest_sites[i]) %>%
    rename(raw_hab = hab) %>%
    left_join(baseline_hab_prop_i) %>%
    mutate(rest_hab_norm = raw_hab - hab_site_base) %>%
    mutate(rest_hab_norm = case_when(rest_hab_norm < 0 ~ 0,
                                  TRUE ~ rest_hab_norm))
  
  # Currently using code that does not add channel habitat
  hab_lookup_rest_all <- bind_rows(hab_lookup_rest_all, lookup_i_2)    
  
}

hab_lookup_rest <- hab_lookup_rest_all %>%
  group_by(cfs, reach) %>%
  summarize(hab = sum(rest_hab_norm)) %>%
  left_join(reach_df) %>%
  bind_rows(no_flow_hab) %>%
  arrange(reach, cfs) %>%
  mutate(hab = case_when(!reach %in% restReach  ~ 0,
                         TRUE ~ hab))

hab_lookup_rest_total <- bind_rows(hab_lookup_baseline, hab_lookup_rest) %>%
  group_by(reach, cfs, reach_num) %>%
  summarize(hab = sum(hab))

# Data for restoration simulation
years <- c(1993:2020)
simQ <- q_sim(years)

start_date <- "05-20"
end_date <- "06-14"
augInput <- 0

# Need augmentation data for functions to work, but set to 0
augQ <- aug_q_fun(start_date, end_date, augInput, simQ)

# New functions to calculate LCC covariate for scenarios

ee <- readRDS("output/ee_list.RData")
expert_num <- 3

# Calculate lcc predictions
lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline) 
lcc_rest <- preds_calc(expert_num, simQ, hab_lookup_rest_total) 
lcc_hydro <- preds_calc(expert_num, augQ, hab_lookup_baseline) 

#Calculate maximum of minimum floodplain inundation to get discharge and habitat projections (not to calculate recruitment) 
maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
maxmins_rest <- maxmins_fun(hab_lookup_rest, simQ, reaches)
maxmins_hydro <- maxmins_fun(hab_lookup_baseline, augQ, reaches)

# Re-scale data for making forecast (done by subtracting the mean of the original covariates)
#Make sure the correct covariates are select for centering
preds_array <- array(NA, dim = c(length(years), length(reaches), 3)) 
preds_array[,,1] <- lcc_baseline - mean(cov_list)
preds_array[,,2] <- lcc_rest - mean(cov_list)
preds_array[,,3] <- lcc_hydro - mean(cov_list)

age0_dat <- forecast_recruit(sim_mod, preds_array, years,
                             maxminsBase = maxmins_baseline, 
                             maxminsAug = maxmins_hydro, 
                             maxminsRest = maxmins_rest, 
                             reaches)

# Filter to just San Acacia
age0_dat_rest <- age0_dat %>% filter(reach == "San Acacia")  %>% mutate(scenario = "rest_sites")

# Berm scenario ####

# Notes on berm changes to hydrology
#Degradation shift the WUA@ 1500 cfs to 4000 cfs).
#Aggredation, shift the WUA @ 1500 cfs to 750 cfs.

# Calculation Aggredation and Degradation q-hab curves

q_inund_curve_combined <- read_csv("data/q_inund_combined.csv") %>%
  select(cfs = q, inund_lower,
         inund_upper,
         wua_lower, wua_upper,
         reach) %>%
  pivot_longer(cols = inund_lower:wua_upper,
               names_to = "hab_type",
               values_to = "hab") %>%
  mutate(reach_num = case_when(reach == "San Acacia" ~ 1,
                               reach == "Isleta" ~ 2,
                               reach == "Angostura" ~ 3))
  
# Q-inund curve using composite curve; shifted for aggradation
q_inund_aggdeg <- q_inund_curve_combined %>% 
  mutate(#deg_cfs = cfs + 2500,
         agg_cfs = case_when (cfs -750 < 0 ~ 0,
                              TRUE ~ cfs - 750))

# What proportion of the reach gets inundated by berm?
#~0.0007-0.001 is an appropriate slope estimate, so 1000-1300 ft of effect per 1 ft berm.
riv_slope <- .0007
berm_height <- 1
n_berms <- 5
length_affected <- 5*berm_height/riv_slope/3.28084

# Read in reach data
sana_reach <- sf::st_read("data/ReachSegments/ReachSegments.shp") %>%
  filter(ReachName == "San Acacia Reach") %>%
  sf::st_transform(crs = 4326)

# Proportion of reach inundated
prop_inund <- length_affected/st_length(sana_reach) %>% as.numeric()

# Habitat in inundated reach
q_agg_initial <- q_inund_aggdeg %>% 
  filter(hab_type == hab_cov) %>%
  select(reach, hab, cfs = agg_cfs) %>%
  mutate(hab = prop_inund*hab)

# Add flat habitat between 7k and 10k cfs (occurs extremely rarely in recent period of record)
q_agg <- q_agg_initial %>% 
  group_by(reach) %>%
  summarize(hab = max(hab)) %>%
  mutate(cfs = 10000) %>%
  bind_rows(q_agg_initial) %>%
  arrange(reach, cfs)

# Get baseline habitat for the proportion of the reach inundated
hab_lookup_baseline_inund_prop <- mutate(hab_lookup_baseline, base_hab = prop_inund*hab) %>% select(-hab)

#Get the total berm-related habitat
hab_lookup_berm_tot <- hab_lookup_fun(q_agg, reaches) %>%
  rename(tot_hab = hab)

# Get berm habitat above baseline level
hab_lookup_berm = left_join(hab_lookup_baseline_inund_prop, hab_lookup_berm_tot) %>%
  mutate(berm_hab = case_when(base_hab > tot_hab ~ 0,
                              is.na(tot_hab) ~ 0,
                              TRUE ~ tot_hab - base_hab))

# Add normalized berm hab. with baseline hab. to calculate overall reach hab. with berms 
hab_lookup_berm_with_baseline <- hab_lookup_berm %>%
  select(cfs, reach, hab = berm_hab, reach_num) %>%
  bind_rows(hab_lookup_baseline) %>%
  group_by(reach, cfs, reach_num) %>%
  summarize(berm_hab = sum(hab)) %>%
  left_join(hab_lookup_baseline)

# Generate flow and augmented flow data
simQ <- q_sim(years)
augQ <- aug_q_fun(start_date, end_date, augInput, simQ)

# Special Berm function
maxmins_hab_temp <- array(NA, dim = c(dim(simQ)[2], 3))
maxmins_cfs_temp <- array(NA, dim = c(dim(simQ)[2], 3))

raw_hab_array <- array(NA, dim = c(dim(simQ)[1], dim(simQ)[2], 3))

# Calculates raw habitat; washes out > 5000 cfs
for (i in 1:length(reaches)){
  hab_lookup_berm_i <- filter(hab_lookup_berm_with_baseline, reach == reaches[i])
  
  if (i == 1) {
    q_temp <- round(simQ[, , 2], 0) %>% as_tibble() %>% set_names("cfs")
  } else {
    q_temp <- round(simQ[, , 1], 0) %>% as_tibble() %>% set_names("cfs")
  }
  
  for (t in 1:dim(simQ)[2]) {
    x_temp <-  q_temp[, t]
    
    hab_vec <- as.data.frame(x_temp) %>%
      left_join(hab_lookup_berm_i) %>%
      mutate(
        below_5k = cumsum(cfs > 5000),
        hab_i = if_else(below_5k >= 1, hab, berm_hab)) %>%
      pull(hab_i)
    
    raw_hab_array[,t,i] <- hab_vec
    
    maxmins_hab_temp[t, i] <- max(rollapplyr(hab_vec, Qdur, mean, fill = NA)[32:122], na.rm = T)
    
    maxmins_cfs_temp[t,i] <- max(rollapplyr(x_temp, Qdur, mean, fill = NA)[32:122], na.rm = T)
    
  }
}

# Alter cov. calc. and pred. calc. function to use raw habitat array
calc_cov_combined_berm <- function(Q, raw_hab_reach_yr, kappa, D, prPARS) {
  thab <- log(raw_hab_reach_yr + 1)
  #for (i in 1:214) {
    #thab[i] <- log(hab_lookup_reach$hab[hab_lookup_reach$cfs == round(Q[i])] + 1)
  #} #finds rel. amount of habitat on given day, based on flow-hab "curve"
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

preds_calc_berm <- function(exp, flows, hab_array) {
  preds <- array(NA, dim = c(dim(flows)[2], 3))
  
  temp<-subset(ee[[exp]],ee[[exp]][,1]=="1"&is.na(ee[[exp]][,5])==FALSE)
  t1<-numeric()
  for (t in 1:6){
    se<-calcsig(temp[t,4],temp[t,3],temp[t,6],temp[t,5]/100)
    t1[t]<-temp[t,6]-se
    se<-calcsig(temp[(t+6),4],temp[(t+6),3],temp[(t+6),6],temp[(t+6),5]/100)
    t1[(t+7)]<-temp[(t+6),6]+se
  }
  t1[7]<-1
  t1[14]<-0
  temp<-subset(ee[[exp]],ee[[exp]][,1]=="2")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t2<-temp[1,6]
  temp<-subset(ee[[exp]],ee[[exp]][,1]=="4")
  se<-calcsig(temp[1,4],temp[1,3],temp[1,6],temp[1,5]/100)
  t4<-temp[1,6]-se
  #
  for (r in 1:3){
    raw_hab_reach <- hab_array[,,r]
    for (j in 1:dim(flows)[2]) {
      t3 <- raw_hab_reach[,j]
      if (r==3) {q<-flows[, j, 1]} else {q<-flows[,j,2]}
      preds[j,r]<-calc_cov_combined_berm(q,t3,t2,t4,t1)
    }}
  
  return(preds)
  
}

lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline) 
lcc_rest <- preds_calc_berm(expert_num, simQ, raw_hab_array) 
lcc_hydro <- preds_calc(expert_num, augQ, hab_lookup_baseline) 

#Calculate maximum of minimum floodplain inundation 
maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
maxmins_rest <- list(maxmins_hab = maxmins_hab_temp, maxmins_cfs = maxmins_cfs_temp) # this is different than normal
maxmins_hydro <- maxmins_fun(hab_lookup_baseline, augQ, reaches)

#Make sure the correct covariates are select for centering
preds_array <- array(NA, dim = c(length(years), length(reaches), 3)) 
preds_array[,,1] <- lcc_baseline - mean(cov_list)
preds_array[,,2] <- lcc_rest - mean(cov_list)
preds_array[,,3] <- lcc_hydro - mean(cov_list)

age0_dat_berm <- forecast_recruit(sim_mod, preds_array, years,
                               maxminsBase = maxmins_baseline, 
                               maxminsAug = maxmins_hydro, 
                               maxminsRest = maxmins_rest, 
                               reaches)

age0_dat_berm_scen <- age0_dat_berm %>% filter(reach == "San Acacia") %>% mutate(scenario = "berm")

# Augmentation scenario ####
# Data for restoration simulation
years <- c(1993:2020)
simQ <- q_sim(years)

# Flow augmentation scenario, volume in acre-ft
augInput <- 10000

# Date windows for adding water
dates <- data.frame(start = seq(as.Date('2011-05-01'),as.Date('2011-05-31'),by = 1)) %>%
  mutate(end = start + 26,
         start_dm = format(start, "%m-%d"),
         end_dm = format(end, "%m-%d")) 

aug_recruit_all <- data.frame()

# Check effects across all flow windows
for (i in 1:nrow(dates)){

augQ <- aug_q_fun(dates$start_dm[i], dates$end_dm[i], augInput, simQ)

lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline) 
lcc_rest <- preds_calc(expert_num, simQ, hab_lookup_rest_total) 
lcc_hydro <- preds_calc(expert_num, augQ, hab_lookup_baseline) 

#Calculate maximum of minimum floodplain inundation 
maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
maxmins_rest <- maxmins_fun(hab_lookup_rest, simQ, reaches)
maxmins_hydro <- maxmins_fun(hab_lookup_baseline, augQ, reaches)

#Make sure the correct covariates are select for centering
preds_array <- array(NA, dim = c(length(years), length(reaches), 3)) 
preds_array[,,1] <- lcc_baseline - mean(cov_list)
preds_array[,,2] <- lcc_rest - mean(cov_list)
preds_array[,,3] <- lcc_hydro - mean(cov_list)

age0_dat <- forecast_recruit(sim_mod, preds_array, years,
                             maxminsBase = maxmins_baseline, 
                             maxminsAug = maxmins_hydro, 
                             maxminsRest = maxmins_rest, 
                             reaches) %>%
  mutate(start_date = dates$start_dm[i])

aug_recruit_all <- bind_rows(aug_recruit_all, age0_dat)

}

# Boxplots of different flow windows
filter(aug_recruit_all, 
       reach == "San Acacia",
       ) %>%
ggplot(aes(as.factor(cfs_baseline), y = prop_hydro_med)) +
  geom_boxplot()

# Figure S2: Plot median flow props for augmentation scenarios
aug_recruit_all %>%
  mutate(year_set = case_when(flow_years %in% c(2013,2020) ~ "2013 and 2020",
                               TRUE ~ "All other years"),
         start_date_num = as.numeric(str_sub(start_date,start = -2, end = -1))) %>%
filter(reach == "San Acacia") %>%
  ggplot(aes(cfs_baseline*0.0283168466, y = ratio_hydro_med, color = start_date_num)) +
  geom_point()+
  facet_wrap(~year_set, scales = "free_y", nrow = 2) +
labs(y = "Proportion of\nbaseline recruitment", 
     x = expression("Spring flow index (m"^{3}*"/s)"),
     color = "May\nstarting\ndate")+
  theme_ipsum()+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        strip.text = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        plot.title = element_text(size = 16, face = "plain"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

#ggsave(filename = "plots/flow_aug_dates.jpeg", width = 6, height = 6, units = "in")

## Assess variability of flow windows ####
# Boxplot of ranks of each year by start date  
aug_recruit_all %>%
  filter(reach == "San Acacia") %>%
  group_by(cfs_baseline) %>%
  mutate(rank = rank(prop_hydro_med)) %>%
ggplot(aes(start_date, rank)) +
  geom_boxplot()

# Find the date with the median effect on recruitment for plotting
med_aug_recruit <- aug_recruit_all %>%
  filter(reach == "San Acacia") %>%
  group_by(cfs_baseline) %>%
  mutate(rank = rank(prop_hydro_med)) %>%
  group_by(start_date) %>%
  summarize(med_rank = median(rank)) %>% 
  arrange(med_rank) %>%
  slice(16) # the 16th value of 31 entries will be the median

med_aug_recruit$start_date

# This is the example flow augmentation data frame for plotting 
# Here we remove the restoration results in the data frame and replacing with flow augmentation scenarios
# This makes later plotting easier when we combine this with berm and restoration data
aug_examp <- aug_recruit_all %>%
  filter(start_date == med_aug_recruit$start_date,
    reach == "San Acacia") %>% 
  mutate(scenario = "Aug") %>%
  select(-c(recruits_rest, diff_recruits_rest, prop_rest_med,
            lower_prop_rest, upper_prop_rest, ratio_rest_med, 
            lower_ratio_rest, upper_ratio_rest)) %>%
  rename(recruits_rest = recruits_hydro, 
         diff_recruits_rest = diff_recruits_hydro, 
         prop_rest_med=prop_hydro_med, 
         lower_prop_rest=lower_prop_hydro, 
         upper_prop_rest=upper_prop_hydro,
         ratio_rest_med=ratio_hydro_med, 
         lower_ratio_rest=lower_ratio_hydro, 
         upper_ratio_rest=upper_ratio_hydro)

# Combine all three scenario dataframes
berm_rest_hydro <- bind_rows(age0_dat_rest, age0_dat_berm_scen, aug_examp) %>%
  mutate(Scenario = case_when(scenario == "rest_sites" ~ "Floodplain\nRestoration",
                                scenario == "Aug" ~ "Flow\nAugmentation",
                                TRUE ~ "Temporary\nBerms"))

# Figure 5 code
pal_df2 <- data.frame(Scenario = unique(berm_rest_hydro$Scenario),
                     color = c("#a6cee3", "#1f78b4","#b2df8a"))

cols2 <- setNames(pal_df2$color, pal_df2$Scenario)

raw_plot <- ggplot(berm_rest_hydro,
       aes(cfs_baseline*0.0283168466, recruits_rest/1000000, color = Scenario)) +
  geom_errorbar(aes(ymin = recruits_baseline/1000000, ymax = recruits_rest/1000000), color = "gray50")+
  geom_point(aes(y = recruits_baseline/1000000), color = "gray50")+
  geom_point()+
  #scale_x_continuous(breaks = c(1000,3000,5000))+
  scale_color_manual(values = cols2)+
  labs(y = "Recruitment\n(millions)",
       x = "Flox index")+
  theme_ipsum()+
  facet_wrap(~Scenario)+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        plot.title = element_text(size = 16, face = "plain"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

raw_plot

diff_plot <- ggplot(berm_rest_hydro,
       aes(cfs_baseline*0.0283168466, y = diff_recruits_rest/1000000, color = Scenario)) +
  geom_point()+
  #scale_x_continuous(breaks = c(1000,3000, 5000))+
  scale_color_manual(values = cols2)+
  labs(y = "Recruitment above\nbaseline (millions)",
       x = "Flox index")+
  theme_ipsum()+
  facet_wrap(~Scenario)+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        plot.title = element_text(size = 16, face = "plain"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

diff_plot

prop_plot <- ggplot(berm_rest_hydro,
       aes(cfs_baseline*0.0283168466, ratio_rest_med, color = Scenario)) +
  geom_point()+
  geom_errorbar(aes(ymin = lower_ratio_rest, ymax = upper_ratio_rest))+
  #scale_x_continuous(breaks = c(1000,3000, 5000))+
  scale_y_continuous(limits = c(1,2.72))+
  scale_color_manual(values = cols2)+
  labs(y = "Proportion of\nbaseline recruitment",
       x = expression("Spring flow index (m"^{3}*"/s)"))+
  theme_ipsum()+
  facet_wrap(~Scenario)+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        legend.position = "none",
        strip.text = element_blank(),
        axis.title.x = element_text(size = 12, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        plot.title = element_text(size = 16, face = "plain"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

prop_plot

raw_plot/plot_spacer()/diff_plot/plot_spacer()/prop_plot + plot_layout(heights =c(4.5, -2.1 ,4.5, -2.1, 4.5))

#ggsave(filename = "plots/scenario_plot_v4_no_channel.jpeg", width = 6, height = 7, units = "in")

# Create Table S1
output <- berm_rest_hydro %>%
  mutate(Discharge = cfs_baseline*0.0283168466) %>%
  select(Scenario, Year = flow_years, Discharge, 
         `Recruits\n(baseline)` = recruits_baseline, `Recruits\n(scenario)` = recruits_rest, 
         `Recruits\nover\nbaseline` = diff_recruits_rest, `Recruit\nratio` = ratio_rest_med) %>%
  mutate(across(`Recruits\n(baseline)`:`Recruits\nover\nbaseline`, \(x) round(x, digits = 0)),
         `Recruit\nratio` = round(`Recruit\nratio`, 2),
         Discharge = round(Discharge, 1))

#write_csv(output, "output/supp_tab_m3.csv")

# Recruitment at individual restoration sites ####
# Same scenario parameters
years <- c(1993:2020)
simQ <- q_sim(years)

start_date <- "05-20"
end_date <- "06-14"
augInput <- 0

augQ <- aug_q_fun(start_date, end_date, augInput, simQ)

restReach <- "San Acacia"

age0_dat_all <- data.frame()

# Loop through restoration sites to calculate recruitment
for (i in 1:length(rest_sites)){
  
  hab_lookup_rest_i <- hab_lookup_rest_all %>%
    filter(site == rest_sites[i]) %>%
    select(cfs, reach, hab = rest_hab_norm) %>%
    left_join(reach_df) %>%
    bind_rows(no_flow_hab) %>%
    arrange(reach, cfs) %>%
    mutate(hab = case_when(!reach %in% restReach  ~ 0,
                           TRUE ~ hab)) %>%
  bind_rows(hab_lookup_baseline) %>%
    group_by(reach, cfs, reach_num) %>%
    summarize(hab = sum(hab))
  
  lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline) 
  lcc_rest <- preds_calc(expert_num, simQ, hab_lookup_rest_i) 
  lcc_hydro <- preds_calc(expert_num, augQ, hab_lookup_baseline) 
  
  #Calculate maximum of minimum floodplain inundation 
  maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
  maxmins_rest <- maxmins_fun(hab_lookup_rest_i, simQ, reaches)
  maxmins_hydro <- maxmins_fun(hab_lookup_baseline, augQ, reaches)
  
  #Make sure the correct covariates are select for centering
  preds_array <- array(NA, dim = c(length(years), length(reaches), 3)) 
  preds_array[,,1] <- lcc_baseline - mean(cov_list)
  preds_array[,,2] <- lcc_rest - mean(cov_list)
  preds_array[,,3] <- lcc_hydro - mean(cov_list)
  
  age0_dat <- forecast_recruit(sim_mod, preds_array, years,
                               maxminsBase = maxmins_baseline, 
                               maxminsAug = maxmins_hydro, 
                               maxminsRest = maxmins_rest, 
                               reaches) %>%
    mutate(site = rest_sites[i])
  
  age0_dat_all <- bind_rows(age0_dat_all, age0_dat)
  
}

# Get info about restoration sites
rest_site_info <- select(q_hab_rest_raw, site, rest_desc, rest_type, n_rest_features, size) %>% distinct()

# Dataframe of San Acacia, with site info
recruits <- age0_dat_all %>%
  filter(reach %in% c("San Acacia")) %>%
  left_join(rest_site_info) %>%
  mutate(site = as.factor(site)) #paste0("RM ", site))

# Plot effects through time
ggplot(recruits)+
  geom_path(aes(flow_years, 
                diff_recruits_rest, 
                color = site))

recruits_min <- slice_min(recruits, recruits_baseline) %>% pull(recruits_baseline)

# Figure S3 Recruits at sites compared to minimum recruitment
ggplot(recruits) +
  geom_violin(aes(site, diff_recruits_rest/1000), outliers = FALSE)+
  geom_jitter(aes(site, diff_recruits_rest/1000), width = .1)+
  geom_hline(aes(yintercept = recruits_min/1000), lty = 2)+
  labs(x = "Restoration site",
       y = "Recruits over baseline (thousands)")+
  theme_ipsum()+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, hjust = 0.5))

#ggsave(filename = "plots/recruit_rest_sites_no_channel.jpeg", width = 6, height = 5, units = "in")

# Figure 6
ggplot(recruits, 
       aes(cfs_baseline*0.0283168466, ratio_rest_med, color = site))+#, shape = rest_type))+
  geom_point(aes(size = size))+
  geom_errorbar(aes(ymin = lower_ratio_rest, ymax = upper_ratio_rest))+
  theme_ipsum()+
  guides(color = "none")+
  facet_wrap(~rest_type)+
  scale_color_colorblind()+
  labs(x=expression("Spring flow index (m"^{3}*"/s)"), y = "Proportion of baseline recruitment", 
       size = "Number of\nrestoration\nfeatures", color = "Site")+
  theme(plot.margin = unit(c(1,0,1,1), "cm"),
        strip.text = element_text(size = 14),
        axis.title.x = element_text(size = 12, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggsave(filename = "plots/recruit_sites_ratio_no_channel.jpeg", width = 6, height = 5, units = "in")
