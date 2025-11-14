# rgsm_hydraulic
Analysis of Rio Grande Silvery minnow recruitment using an integrated population model 
and hydraulic model output.

This analysis modifies a Rio Grande Silvery minnow population model developed by 
Yackulic et al. 2022 (https://doi.org/10.1002/ecs2.4240; https://github.com/cyack/rgsm_integrated)
by integrating it with a novel hydaulic model to estimate recruitment at the site-scale.

Below are brief descriptions of the code and data in this project.

## R and Stan code required to reproduce results and figures in manuscript
1_data_setup.R – R code for cleaning and organizing silvery minnow and environmental data
2_run_models.R - R code for fitting population models and evaluating using out-of-sample data
3_scenarios.R - R code for assessing recruitment in different scenarios using fitted models
int3re.stan – Stan code for fitting integrated model with time varying mortality rates

## Data files containing information required to reproduce results and figures in manuscript

### Data from original analysis for model fitting (Yackulic et al. 2022)
abq_gage_08330000.csv – daily flow data for USGS gage 08330000 in long format
aQ.csv – daily flow data for USGS gage 08330000 from March 1st to September 30th organized by year
ASIRQuery4Combined.csv – RGSM monitoring dataset - subset for analysis, updated version available online
Fish Rescue 2009-2020.csv – Fish rescue data
fws2br_rmconverter.csv - file for converting between river mile system used by USFWS and USBOR
mesohab_sum.csv – summary of habitat measurements from Braun et al., 2015
oos_catch.csv – catch data for October 2019 and October 2020.
oos_releases.csv – information on augmentation relevant to out of sample analysis
oos_rivereyes.csv – information on river drying during 2019 and 2020
releases 8_26_2019.csv – information on augmentation from 2002 to 2018
rescue_surv.csv – information on survival of rescued fish from Archdeacon et al., 2020
rivereyes_v2.csv - information on river drying from 2002 to 2018
SanAcacia_gage_0835490.csv 0 – daily flow data for USGS gage 08354900 in long format
sQ.csv – daily flow data for USGS gage 08354900from March 1st to September 30th organized by year
sum_ee1.csv - Expert elicitation results from expert 1
sum_ee2.csv - Expert elicitation results from expert 2
sum_ee3.csv - Expert elicitation results from expert 3
sum_ee4.csv - Expert elicitation results from expert 4
sum_ee5.csv - Expert elicitation results from expert 5

### Hydraulic-model related data

ecoval2d.csv - Hydraulic model output relating discharge to inundated acres of floodplain
q_inund_combined.csv - Relationship between discharge and floodplain habitat with low flow data informed by expert elicitation     
q_hab_lookup_combined.csv - Interpolated lookup table of q_inund_combined.csv   
ReachSegments.shp - Shapefile of Middle Rio Grande reaches              
san_acacia_rest_ecovals.csv - Hydraulic model output for 8 San Acacia restoration sites, from Harris et al. 2021 
(https://webapps.usgs.gov/mrgescp/documents/Harris_2021_Hydraulic-Habitat-Suitability-for-RGSM-at-San-Acacia-Restoration-Sites.pdf)

### New out-of-sample data
dry_eyes_download_reyes.csv - Updated river drying data
oosR_raw.csv - Updated rescue data
PopMon_MonthlyHaulUSBR.csv - Updated RGSM monitoring data

## Output from data preparation and modeling

### Output derived from input data preparation
ee_list.RData                                           
input_data.RData                
oos_data_new.RData   

### Covariates derived for different model versions for scenarios
inund_cov_combined_list.RData  
inund_cov_list.RData            
inund_lcc_list_combined.RData  
inund_lcc_list.RData            
 
### Model results to be used in scenarios
inund_low_combined_mcmc_sub.csv
inund_low_mcmc_sub.csv          
lcc_low_combined_mcmc_sub.csv  
lcc_low_mcmc_sub.csv                      
orig_lcc_mcmc_sub.csv           
springflow_mcmc_sub.csv        

### Table of scenario results
supp_tab_m3.csv                

