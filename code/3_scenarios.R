# =============================================================================
# Script: 3_scenarios.R
# Purpose: Generate RGSM age-0 recruitment projections under three management
#          scenarios — (1) floodplain restoration, (2) temporary berms, and
#          (3) environmental flow augmentation — and produce manuscript figures
#          comparing scenario outcomes across all historical flow years.
#
# Overview of workflow:
#   1. Define helper functions used throughout the script
#   2. Load and prepare flow, hydraulic, and model data
#   3. Restoration scenario: build site-specific habitat lookups, calculate
#      the larval carrying capacity (LCC) covariate, run recruitment model,
#      save results
#   4. Berm scenario: compute berm-modified habitat with flow-dependent
#      washout rule, run modified covariate and recruitment functions, save
#      results
#   5. Flow augmentation scenario: loop over 31 possible start dates in May,
#      run full model for each, identify the median-effect date for plotting
#   6. Combine all three scenarios and produce manuscript figures (Fig. 4,
#      Fig. 5, Fig. S1, Fig. S2, Fig. S3, Table S1)
#   7. Loop over individual restoration sites to assess site-level recruitment
#
# Key dependencies (all must be installed before running):
#   tidyverse, zoo, hrbrthemes, sf, patchwork, ggthemes
#
# Key input files (paths relative to project root):
#   data/yackulic2022_data/abq_gage_08330000.csv       -- Albuquerque USGS gage
#   data/yackulic2022_data/SanAcacia_gage_08354900.csv -- San Acacia USGS gage
#   data/ecoval2d.csv                                  -- 2D hydraulic model output
#   data/q_hab_lookup_combined.csv                     -- Interpolated hab lookup
#   data/san_acacia_rest_ecovals.csv                   -- Restoration site hydraulics
#   data/q_inund_combined.csv                          -- Combined inundation curves
#   data/ReachSegments/ReachSegments.shp               -- Reach geometry shapefile
#   output/inund_lcc_list_combined.RData               -- Historical LCC covariates
#   output/lcc_low_combined_mcmc_sub.csv               -- MCMC posterior samples
#   output/ee_list.RData                               -- Expert elicitation values
#
# Key output files:
#   plots/fig2_ecovalue_curve.jpeg
#   plots/figS1_rest_q_hab.jpeg
#   plots/figS2_flow_aug_dates.jpeg
#   plots/figS3_recruit_rest_sites_no_channel.jpeg
#   plots/fig4_scenario_plot_v4_no_channel.jpeg
#   plots/fig5_recruit_sites_ratio_no_channel.jpeg
#   output/supp_tab_m3_new_v2.csv
# =============================================================================

library(tidyverse)
library(zoo)        # rollapplyr() for rolling window calculations
library(hrbrthemes) # theme_ipsum() for publication-style plot themes
library(sf)         # Spatial operations on reach geometry shapefile
library(patchwork)  # Compositing multiple ggplots into a single figure
library(ggthemes)   # scale_color_colorblind() for accessible color palettes


# =============================================================================
# SECTION 1: HELPER FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# calcsig()
# Purpose: Recover the implied standard deviation from an expert-elicited
#   probability statement. Given that an expert states their central estimate
#   (mean) and that the true value falls within (down, up) with probability q,
#   this function finds the standard deviation of the normal distribution
#   consistent with those three inputs using numerical optimization.
#
#   This is used in preds_calc() to reconstruct the full distributional shape
#   of each elicited parameter from the elicitation data in ee_list.RData.
#
# Args:
#   up   -- upper bound of elicited confidence interval
#   down -- lower bound of elicited confidence interval
#   mean -- elicited central estimate
#   q    -- stated probability that true value falls within (down, up), e.g. 0.9
#   INT  -- search interval for the standard deviation; default c(0, 1000)
#
# Returns: numeric scalar -- the recovered standard deviation
# -----------------------------------------------------------------------------
calcsig <- function(up, down, mean, q, INT = c(0, 1000)) {
  sig <- function(x) {
    # Find sigma such that P(down < X < up) = q under N(mean, sigma)
    abs(pnorm(up, mean, x) - pnorm(down, mean, x) - q)
  }
  optimize(sig, interval = INT)$minimum
}

# -----------------------------------------------------------------------------
# lin_int()
# Purpose: Linear interpolation between two points. Used throughout to
#   interpolate habitat values at discharge levels between hydraulic model
#   breakpoints (e.g., between 700 and 1500 cfs).
#
# Args:
#   x           -- the x value at which to interpolate
#   xlow, xhi   -- x bounds of the known interval
#   ylow, yhi   -- y values at xlow and xhi respectively
#
# Returns: numeric -- interpolated y value at x
# -----------------------------------------------------------------------------
lin_int <- function(x, xlow, xhi, ylow, yhi) {
  (yhi - ylow) * (x - xlow) / (xhi - xlow) + ylow
}

# Shorthand quantile functions used in forecast_recruit() to summarize
# posterior distributions. q10/q90 define the 80% credible interval.
q025 <- function(x) {quantile(x, .025)}
q975 <- function(x) {quantile(x, .975)}
q10  <- function(x) {quantile(x, .1)}
q90  <- function(x) {quantile(x, .9)}


# -----------------------------------------------------------------------------
# q_sim()
# Purpose: Build a 3D array of daily discharge values for simulation.
#   Dimension 1: day index (214 days, March 1 to September 30)
#   Dimension 2: flow scenario (one column per element of flow_input)
#   Dimension 3: gage (1 = Angostura/Albuquerque gage, 2 = San Acacia gage)
#
# Accepts two types of input:
#   - Values between 0 and 1 (quantiles): selects the historical year whose
#     April-June total volume most closely matches that quantile of the full
#     historical distribution. Enables stylized "wet" / "dry" year simulations.
#   - Integer year values (>= 1): extracts the observed hydrograph for that
#     calendar year directly from the gage records.
#
# The March-September window (month > 2, month < 10) covers the full RGSM
# spawning and larval rearing period.
#
# Args:
#   flow_input -- numeric vector of years (e.g., 1993:2020) or quantiles
#                 (e.g., c(0.1, 0.5, 0.9))
#
# Returns: 3D array [214 days x n_flows x 2 gages]
# -----------------------------------------------------------------------------
q_sim <- function(flow_input) {
  q_all <- array(NA, dim = c(214, length(flow_input), 2))
  
  if (length(flow_input[flow_input < 1]) > 0) {
    # Quantile mode: match each quantile to the closest historical year
    for (i in 1:length(flow_input)) {
      # Find the year whose Apr-Jun volume is closest to the requested quantile,
      # separately for each gage (gages may map to different years)
      aQ_quant   <- quantile(aQ_med_spr$tot_vol, flow_input[i])
      a_year_sel <- aQ_med_spr$year[which.min(abs(aQ_quant - aQ_med_spr$tot_vol))]
      
      sQ_quant   <- quantile(sQ_med_spr$tot_vol, flow_input[i])
      s_year_sel <- sQ_med_spr$year[which.min(abs(sQ_quant - sQ_med_spr$tot_vol))]
      
      aQ_temp <- filter(angQ,  year == a_year_sel, month > 2, month < 10)
      sQ_temp <- filter(sanaQ, year == s_year_sel, month > 2, month < 10)
      
      q_all[, i, 1] <- aQ_temp$cfs
      q_all[, i, 2] <- sQ_temp$cfs
    }
    
  } else {
    # Year mode: extract observed hydrograph directly for each year
    for (i in 1:length(flow_input)) {
      aQ_temp <- filter(angQ,  year == flow_input[i], month > 2, month < 10)
      sQ_temp <- filter(sanaQ, year == flow_input[i], month > 2, month < 10)
      
      q_all[, i, 1] <- aQ_temp$cfs
      q_all[, i, 2] <- sQ_temp$cfs
    }
  }
  return(q_all)
}


# =============================================================================
# SECTION 1B: LARVAL CARRYING CAPACITY (LCC) COVARIATE FUNCTIONS
#
# These three functions work together to compute the LCC covariate for a given
# year and reach. The LCC covariate is the key link between flow conditions and
# RGSM recruitment in the Beverton-Holt model.
#
# Conceptual flow:
#   daily Q --> q2hab_2d() --> daily habitat (thab, log-scale)
#                                    |
#          expert elicitation --> prop() --> daily egg-laying weights (tprop)
#                                    |
#                         calc_cov_combined() --> LCC covariate (scalar)
# =============================================================================

# -----------------------------------------------------------------------------
# q2hab_2d()
# Purpose: Linearly interpolate inundated habitat (acres) at an arbitrary
#   discharge value, given the hydraulic model breakpoints from ecoval2d.csv.
#   "2d" refers to the two-dimensional HEC-RAS model that produced the values.
#
#   Appends the last pars value and 20,000 cfs as a ceiling to hold habitat
#   flat above the highest modeled discharge, preventing extrapolation.
#
# Args:
#   q     -- discharge value(s) in cfs to interpolate at
#   pars  -- habitat values (acres) at each breakpoint (e.g., inund_lower)
#   input -- breakpoint discharge values in cfs (e.g., c(700,1000,1500,...))
#
# Returns: numeric -- interpolated habitat in acres at discharge q
# -----------------------------------------------------------------------------
q2hab_2d <- function(q, pars, input) {
  p  <- c(pars, pars[length(pars)])  # hold flat above highest breakpoint
  qs <- c(input, 2 * 10 ^ 4)         # 20,000 cfs ceiling
  
  t1 <- findInterval(q, qs)          # find which interval q falls in
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}

# -----------------------------------------------------------------------------
# prop()
# Purpose: Compute the daily proportional egg-laying distribution across
#   the 132-day spawning window (days 1-132, approximately March 1 to July 10)
#   based on expert elicitation of RGSM spawning phenology.
#
#   Experts elicited relative spawning amounts at 14 time points spaced 10
#   days apart (days 2, 12, 22, ..., 132). This function interpolates between
#   those values to produce a smooth daily profile, normalized to sum to 1.
#
# Args:
#   pars -- 14-element vector of elicited relative spawning values at 10-day
#           breakpoints (from ee_list.RData, question type "1")
#
# Returns: numeric vector of length 132 -- proportional egg-laying by day
# -----------------------------------------------------------------------------
prop <- function(pars) {
  t    <- seq(2, 132, 10)  # 14 breakpoint days
  tout <- pars[1]           # day 1 initialized to first elicited value
  
  for (i in 2:131) {
    t1      <- findInterval(i, t)
    tout[i] <- lin_int(i, t[t1], t[(t1 + 1)], pars[t1], pars[(t1 + 1)])
  }
  tout[132] <- pars[14]    # day 132 uses last elicited value directly
  
  tout / sum(tout)          # normalize to proportions summing to 1
}

# -----------------------------------------------------------------------------
# calc_cov_combined()
# Purpose: Compute the LCC covariate for a single year and reach. This is the
#   core calculation linking daily discharge to the recruitment model.
#
# Steps:
#   1. Convert daily Q to log(habitat + 1) via the reach-specific lookup table
#   2. For each day in the spawning window, find the minimum habitat over the
#      next D days (the "duration bottleneck" from expert elicitation)
#   3. Identify tstart: the first day spawning is triggered, defined as the
#      earliest of (a) Q > kappa, (b) Q rises >100 cfs in one day, (c) day 132
#   4. Compute a weighted sum of minimum habitat, where weights are the daily
#      egg-laying proportions from prop()
#
# The resulting scalar (LCC covariate) is centered by subtracting mean(cov_list)
# before entering the Beverton-Holt model as the predictor of carrying capacity.
#
# Args:
#   Q               -- numeric vector [214] of daily discharge in cfs
#   hab_lookup_reach -- tibble (cfs integer 0-10000, hab acres); flow-habitat
#                       lookup for this reach at 1-cfs resolution
#   kappa           -- flow cue threshold (cfs); from EE question type "2"
#   D               -- duration requirement (days); Qdur = 20 in deployed model
#   prPARS          -- 14-element elicited spawning proportions; input to prop()
#
# Returns: numeric scalar -- the LCC covariate value
# -----------------------------------------------------------------------------
calc_cov_combined <- function(Q, hab_lookup_reach, kappa, D, prPARS) {
  
  # Step 1: Look up log(habitat+1) for each of the 214 simulation days.
  # Integer rounding of Q matches the 1-cfs resolution of the lookup table.
  thab <- numeric()
  for (i in 1:214) {
    thab[i] <- log(hab_lookup_reach$hab[hab_lookup_reach$cfs == round(Q[i])] + 1)
  }
  
  tprop <- prop(prPARS)  # daily egg-laying proportions
  
  # Step 2: For each spawning day, find the minimum habitat over the next D days.
  # This captures the "bottleneck" constraint -- larvae need D consecutive days
  # of inundation for sufficient survival.
  t2hab <- numeric()
  for (i in 1:132) {
    t2hab[i] <- min(thab[i:(i + D)])
  }
  
  # Step 3: Spawning trigger -- first day any of three conditions is met:
  #   (a) discharge exceeds kappa (flow cue threshold)
  #   (b) discharge rises by more than 100 cfs in a single day (rising limb)
  #   (c) day 132 (end of window; spawning occurs regardless)
  tstart <- min(c(which(Q > kappa),
                  which(Q[-1] - Q[-length(Q)] > 100),
                  132))
  
  # Step 4: Weighted sum of minimum habitat across the spawning window.
  # Eggs laid before tstart all "experience" the habitat at tstart (they hatch
  # when conditions are first met). Eggs laid from tstart onward experience
  # habitat on their own day.
  out <- sum(tprop[1:tstart]) * t2hab[tstart] +
    sum(tprop[tstart:132] * t2hab[tstart:132])
  return(out)
}

# -----------------------------------------------------------------------------
# preds_calc()
# Purpose: Compute the LCC covariate for all flow years and all reaches,
#   returning a matrix [n_years x 3 reaches].
#
#   Reconstructs expert elicitation parameters (t1, t2, t4) from ee_list.RData
#   using calcsig(). Each elicited value is stored with its central estimate,
#   confidence bounds, and stated containment probability; calcsig() recovers
#   the implied SD, which is added or subtracted to get the point estimate.
#
#   Reach-gage assignment:
#     r == 1 (San Acacia): San Acacia gage [array dim 2]
#     r == 2 (Isleta):     San Acacia gage [array dim 2]
#     r == 3 (Angostura):  Albuquerque gage [array dim 1]
#
# Args:
#   exp        -- integer; expert index (expert_num = 3 in deployed model)
#   flows      -- 3D array [214 x n_years x 2] from q_sim()
#   hab_lookup -- tibble (cfs, reach, hab, reach_num); all-reach lookup table
#
# Returns: matrix [n_years x 3] of LCC covariate values
# -----------------------------------------------------------------------------
preds_calc <- function(exp, flows, hab_lookup) {
  preds <- array(NA, dim = c(dim(flows)[2], 3))
  
  # Reconstruct t1: 14-element spawning phenology vector (EE question type "1")
  # For 6 interior breakpoints, calcsig() recovers the implied SD; subtracted
  # for lower breakpoints, added for upper breakpoints.
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "1" & is.na(ee[[exp]][, 5]) == FALSE)
  t1 <- numeric()
  for (t in 1:6) {
    se         <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t1[t]      <- temp[t, 6] - se
    se         <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100)
    t1[(t + 7)] <- temp[(t + 6), 6] + se
  }
  t1[7]  <- 1  # peak relative spawning = 1 (fixed maximum)
  t1[14] <- 0  # end of window = 0 (spawning ceases)
  
  # t2: flow cue kappa (EE question type "2"); central value used directly
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "2")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t2   <- temp[1, 6]
  
  # t4: duration requirement D (EE question type "4"); lower bound (central - SD)
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "4")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t4   <- temp[1, 6] - se
  
  # Compute LCC covariate for each reach and year
  for (r in 1:3) {
    t3 <- filter(hab_lookup, reach_num == r)
    for (j in 1:dim(flows)[2]) {
      # Angostura uses ABQ gage (dim 1); San Acacia and Isleta use SA gage (dim 2)
      if (r == 3) {q <- flows[, j, 1]} else {q <- flows[, j, 2]}
      preds[j, r] <- calc_cov_combined(q, t3, t2, t4, t1)
    }
  }
  
  return(preds)
}


# =============================================================================
# SECTION 1C: HABITAT LOOKUP AND INUNDATION COVARIATE FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# hab_lookup_fun()
# Purpose: Build a flow-habitat lookup table at 1-cfs resolution (1 to 10,000
#   cfs) for each reach, by linearly interpolating between the hydraulic model
#   breakpoints using q2hab_2d().
#
#   If the input curve has a "site" column (e.g., for restoration parcels),
#   loops over each site within a reach and sums habitat across sites to get
#   total reach-level habitat at each flow.
#
# Args:
#   inund_curve -- tibble with columns: cfs, hab, reach, and optionally site
#   reaches     -- character vector of reach names to process
#
# Returns: tibble with columns cfs (1:10000), hab (acres), reach
# -----------------------------------------------------------------------------
hab_lookup_fun <- function(inund_curve, reaches) {
  hab_lookup <- tibble()
  
  for (i in 1:length(reaches)) {
    q_hab_i      <- filter(inund_curve, reach == reaches[i])
    sites        <- unique(q_hab_i$site)
    hab_lookup_i <- tibble()
    
    for (j in 1:length(sites)) {
      # If site column exists, filter to current site; otherwise use all rows
      q_hab_ij <- if ("site" %in% colnames(q_hab_i)) {
        q_hab_i %>% filter(site %in% sites[j])
      } else {
        q_hab_i
      }
      
      # Interpolate to 1-cfs resolution from breakpoints
      hab_lookup_ij <- tibble(cfs = 1:10e3) %>%
        mutate(hab = q2hab_2d(cfs, q_hab_ij$hab, q_hab_ij$cfs))
      
      hab_lookup_i <- bind_rows(hab_lookup_i, hab_lookup_ij)
    }
    
    # Sum across sites within each reach
    hab_lookup_i <- group_by(hab_lookup_i, cfs) %>%
      summarize(hab = sum(hab)) %>%
      mutate(reach = reaches[i])
    
    hab_lookup <- bind_rows(hab_lookup, hab_lookup_i)
  }
  return(hab_lookup)
}

# -----------------------------------------------------------------------------
# maxmins_fun()
# Purpose: Compute the "flow index" (x-axis in recruitment plots) and peak
#   inundation metrics for each simulated year and reach.
#
#   For each year/reach, finds the maximum Qdur-day rolling mean of habitat
#   (or discharge) within the spawning window (days 32-122, ~April 1 to
#   July 1). This metric captures the best sustained inundation event in a
#   given year and is used for plotting and output CSVs. It is NOT the
#   same as the LCC covariate used for Beverton-Holt predictions.
#
# Args:
#   lookup  -- hab_lookup tibble at 1-cfs resolution
#   sim_q   -- 3D flow array [214 x n_years x 2] from q_sim()
#   reaches -- character vector of reach names
#
# Returns: list with two matrices [n_years x 3]:
#   maxmins_hab -- peak Qdur-day mean inundated habitat (acres)
#   maxmins_cfs -- peak Qdur-day mean discharge (cfs); the flow index
# -----------------------------------------------------------------------------
maxmins_fun <- function(lookup, sim_q, reaches) {
  maxmins_hab <- array(NA, dim = c(dim(sim_q)[2], 3))
  maxmins_cfs <- array(NA, dim = c(dim(sim_q)[2], 3))
  
  for (i in 1:length(reaches)) {
    hab_lookup_i <- filter(lookup, reach == reaches[i])
    
    # San Acacia (i==1) uses San Acacia gage; Isleta and Angostura use ABQ gage
    if (i == 1) {
      q <- round(sim_q[, , 2], 0) %>% as_tibble() %>% set_names("cfs")
    } else {
      q <- round(sim_q[, , 1], 0) %>% as_tibble() %>% set_names("cfs")
    }
    
    for (t in 1:dim(sim_q)[2]) {
      x <- q[, t]
      
      # Match daily discharge to habitat, compute rolling mean, take the maximum
      # within the spawning window (days 32-122)
      maxmins_hab[t, i] <- max(
        rollapplyr(hab_lookup_i$hab[match(x$cfs, hab_lookup_i$cfs)],
                   Qdur, mean, fill = NA)[32:122],
        na.rm = TRUE)
      
      maxmins_cfs[t, i] <- max(
        rollapplyr(x$cfs, Qdur, mean, fill = NA)[32:122],
        na.rm = TRUE)
    }
  }
  
  return(list(maxmins_hab = maxmins_hab, maxmins_cfs = maxmins_cfs))
}


# -----------------------------------------------------------------------------
# aug_q_fun()
# Purpose: Add an environmental flow pulse to a simulated hydrograph array.
#   Converts total volume (acre-feet) to a constant daily CFS rate added
#   uniformly across the specified date window, then adds to all years.
#
#   Conversion: vol(ac-ft) * 43,559.935 ft3/ac-ft / ndays / 86,400 sec/day = cfs
#
# Args:
#   st_date  -- start date string "MM-DD"
#   en_date  -- end date string "MM-DD"
#   vol_AcFt -- total augmentation volume in acre-feet
#   q_array  -- 3D flow array [214 x n_years x 2] from q_sim()
#
# Returns: 3D array with augmentation added to both gage dimensions
# -----------------------------------------------------------------------------
aug_q_fun <- function(st_date, en_date, vol_AcFt, q_array) {
  
  # Reference table: MM-DD strings mapped to day indices 1-214 (Mar1=1, Sep30=214)
  simq_dates <- data.frame(
    date    = seq(as.Date("2000-03-01"), as.Date("2000-09-30"), 1),
    day_ind = 1:214
  ) %>%
    mutate(date = as.character(format(date, "%m-%d")))
  
  start_ind <- simq_dates$day_ind[simq_dates$date == st_date]
  end_ind   <- simq_dates$day_ind[simq_dates$date == en_date]
  
  aug_vec   <- rep(0, 214)
  ndays_aug <- length(start_ind:end_ind)
  
  augCFS <- round(vol_AcFt * 43559.935 / ndays_aug / (24 * 60 * 60), 0)
  aug_vec[start_ind:end_ind] <- rep(augCFS, ndays_aug)
  
  aug_q_array <- q_array + aug_vec
  return(aug_q_array)
}


# =============================================================================
# SECTION 1D: RECRUITMENT FORECAST FUNCTION
# =============================================================================

# -----------------------------------------------------------------------------
# forecast_recruit()
# Purpose: Apply the Beverton-Holt recruitment model across all flow years,
#   MCMC posterior samples, and three reaches, under three scenarios.
#
# Beverton-Holt equation for reach k, year i, MCMC sample m:
#
#   age0[m,i,k] = a[m] * effS[m,k] / (1 + a[m] * effS[m,k] / tRf[m,i,k])
#
#   where tRf[m,i,k] = exp(B_lbeta * Xout[i,k,s]) * tRf0[m,k]
#         tRf0[m,k]  = exp(rnorm(iter, mu_lbeta[k], sd_lbeta))
#
# As tRf -> inf (unlimited habitat), age0 -> a * effS (density-independent)
# As tRf -> 0   (no habitat),       age0 -> tRf (fully habitat-limited)
#
# Spawner selection: spawn_rank (global, = 9 = median) selects a single
# historical year of effective spawner data to hold spawners constant across
# scenarios, isolating the effect of habitat differences on recruitment.
#
# Args:
#   mod         -- tibble of 15,000 MCMC posterior samples (sim_mod)
#   Xout        -- array [n_years x 3 x 3] of centered LCC covariates
#                  (dim3: 1=baseline, 2=restoration, 3=flow augmentation)
#   years_sim   -- numeric vector of simulated years
#   maxminsBase/Aug/Rest -- output of maxmins_fun() for each scenario
#   reaches     -- character vector of reach names
#
# Returns: tibble (n_years x n_reaches rows) with median, 10th, 90th
#   percentile recruitment under all three scenarios plus inundation metrics
# -----------------------------------------------------------------------------
forecast_recruit <- function(mod, Xout, years_sim,
                             maxminsBase, maxminsAug, maxminsRest, reaches) {
  
  # Extract MCMC parameter vectors
  a        <- pull(mod, a)
  mu_lbeta <- select(mod, starts_with("mu_lbeta")) %>% as.matrix()  # [iter x 3]
  sd_lbeta <- pull(mod, sd_lbeta)
  B_lbeta  <- pull(mod, B_lbeta)
  effS_all <- select(mod, starts_with("effS"))
  iter     <- length(a)
  
  # Reshape effS to long format and identify which historical year has the
  # target spawner rank (spawn_rank = 9 = median across 17 fitting years)
  effS_long <- effS_all %>%
    mutate(sample = 1:15000) %>%
    gather(-sample, key = "param", value = "effS") %>%
    mutate(year  = rep(rep(c(1:17), each = iter), 3),
           reach = rep(c(1:3), each = iter * 17))
  
  effS_medians <- effS_long %>%
    group_by(year, reach) %>%
    summarize(median_effS = median(effS)) %>%
    spread(key = reach, value = median_effS) %>%
    ungroup() %>%
    select(-year)
  
  effS_rank <- apply(effS_medians, 2, rank)
  effS_ind  <- which(effS_rank == spawn_rank, arr.ind = TRUE)[, 1]
  
  # Pre-allocate output arrays [iter x n_years x 3 reaches]
  age0_base   <- array(NA, dim = c(iter, nrow(Xout), 3))
  age0_rest   <- array(NA, dim = c(iter, nrow(Xout), 3))
  age0_hydro  <- array(NA, dim = c(iter, nrow(Xout), 3))
  diff_rest   <- array(NA, dim = c(iter, nrow(Xout), 3))
  prop_rest   <- array(NA, dim = c(iter, nrow(Xout), 3))
  diff_hydro  <- array(NA, dim = c(iter, nrow(Xout), 3))
  prop_hydro  <- array(NA, dim = c(iter, nrow(Xout), 3))
  ratio_rest  <- array(NA, dim = c(iter, nrow(Xout), 3))
  ratio_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  effS_sim    <- array(NA, dim = c(iter, 3))
  
  for (k in 1:3) {
    # Pull all 15,000 MCMC samples of effS for the selected spawner year/reach
    effS <- filter(effS_long, year == effS_ind[k], reach == k) %>% pull(effS)
    
    for (i in 1:nrow(Xout)) {
      # tRf0: baseline carrying capacity without inundation effect.
      # Drawn from N(mu_lbeta[k], sd_lbeta) on the log scale.
      tRf0 <- exp(rnorm(iter, mu_lbeta[, k], sd_lbeta))
      
      # tRf: carrying capacity amplified by LCC covariate (additive on log scale)
      tRf_base  <- exp(B_lbeta * Xout[i, k, 1]) * tRf0
      tRf_rest  <- exp(B_lbeta * Xout[i, k, 2]) * tRf0
      tRf_hydro <- exp(B_lbeta * Xout[i, k, 3]) * tRf0
      
      # Beverton-Holt: R = a*S / (1 + a*S/tRf)
      age0_base[, i, k]  <- a * effS / (1 + a * effS / tRf_base)
      age0_rest[, i, k]  <- a * effS / (1 + a * effS / tRf_rest)
      age0_hydro[, i, k] <- a * effS / (1 + a * effS / tRf_hydro)
      
      # Differences and ratios relative to baseline (across all MCMC iterations)
      diff_rest[, i, k]   <- age0_rest[, i, k]  - age0_base[, i, k]
      diff_hydro[, i, k]  <- age0_hydro[, i, k] - age0_base[, i, k]
      prop_hydro[, i, k]  <- diff_hydro[, i, k] / age0_base[, i, k]
      prop_rest[, i, k]   <- diff_rest[, i, k]  / age0_base[, i, k]
      ratio_hydro[, i, k] <- age0_hydro[, i, k] / age0_base[, i, k]
      ratio_rest[, i, k]  <- age0_rest[, i, k]  / age0_base[, i, k]
    }
    effS_sim[, k] <- effS
  }
  
  # Summarize posterior: median and 80% interval (q10/q90) for each year/reach.
  # apply() over c(2,3) marginalizes over MCMC iterations (dimension 1).
  age0_baseline    <- apply(age0_base,  c(2, 3), median) %>% as.vector()
  lower_age0       <- apply(age0_base,  c(2, 3), q10)    %>% as.vector()
  upper_age0       <- apply(age0_base,  c(2, 3), q90)    %>% as.vector()
  age0_rest_med    <- apply(age0_rest,  c(2, 3), median) %>% as.vector()
  lower_age0_rest  <- apply(age0_rest,  c(2, 3), q10)    %>% as.vector()
  upper_age0_rest  <- apply(age0_rest,  c(2, 3), q90)    %>% as.vector()
  age0_hydro_med   <- apply(age0_hydro, c(2, 3), median) %>% as.vector()
  lower_age0_hydro <- apply(age0_hydro, c(2, 3), q10)    %>% as.vector()
  upper_age0_hydro <- apply(age0_hydro, c(2, 3), q90)    %>% as.vector()
  diff_rest_med    <- apply(diff_rest,  c(2, 3), median) %>% as.vector()
  lower_diff_rest  <- apply(diff_rest,  c(2, 3), q10)    %>% as.vector()
  upper_diff_rest  <- apply(diff_rest,  c(2, 3), q90)    %>% as.vector()
  diff_hydro_med   <- apply(diff_hydro, c(2, 3), median) %>% as.vector()
  lower_diff_hydro <- apply(diff_hydro, c(2, 3), q10)    %>% as.vector()
  upper_diff_hydro <- apply(diff_hydro, c(2, 3), q90)    %>% as.vector()
  prop_rest_med    <- apply(prop_rest,  c(2, 3), median) %>% as.vector()
  lower_prop_rest  <- apply(prop_rest,  c(2, 3), q10)    %>% as.vector()
  upper_prop_rest  <- apply(prop_rest,  c(2, 3), q90)    %>% as.vector()
  prop_hydro_med   <- apply(prop_hydro, c(2, 3), median) %>% as.vector()
  lower_prop_hydro <- apply(prop_hydro, c(2, 3), q10)    %>% as.vector()
  upper_prop_hydro <- apply(prop_hydro, c(2, 3), q90)    %>% as.vector()
  ratio_rest_med   <- apply(ratio_rest,  c(2, 3), median) %>% as.vector()
  lower_ratio_rest <- apply(ratio_rest,  c(2, 3), q10)    %>% as.vector()
  upper_ratio_rest <- apply(ratio_rest,  c(2, 3), q90)    %>% as.vector()
  ratio_hydro_med   <- apply(ratio_hydro, c(2, 3), median) %>% as.vector()
  lower_ratio_hydro <- apply(ratio_hydro, c(2, 3), q10)    %>% as.vector()
  upper_ratio_hydro <- apply(ratio_hydro, c(2, 3), q90)    %>% as.vector()
  
  # Assemble output tibble: one row per (flow year x reach).
  # rep(reaches, each = n_years) repeats each reach name n_years times so
  # that filtering by reach gives a contiguous block for that reach.
  age0_dat <- tibble(
    reach                = rep(reaches, each = nrow(Xout)),
    flow_years           = rep(years_sim, 3),
    cfs_baseline         = as.vector(maxminsBase$maxmins_cfs),  # x-axis flow index
    cfs_augmentation     = as.vector(maxminsAug$maxmins_cfs),
    acres_inund_baseline = as.vector(maxminsBase$maxmins_hab),
    acres_inund_rest     = as.vector(maxminsRest$maxmins_hab),
    acres_inund_aug      = as.vector(maxminsAug$maxmins_hab),
    recruits_baseline    = age0_baseline,
    lower_age0           = lower_age0,
    upper_age0           = upper_age0,
    recruits_rest        = age0_rest_med,
    lower_age0_rest      = lower_age0_rest,
    upper_age0_rest      = upper_age0_rest,
    recruits_hydro       = age0_hydro_med,
    lower_age0_hydro     = lower_age0_hydro,
    upper_age0_hydro     = upper_age0_hydro,
    diff_recruits_rest   = diff_rest_med,
    lower_diff_rest      = lower_diff_rest,
    upper_diff_rest      = upper_diff_rest,
    diff_recruits_hydro  = diff_hydro_med,
    lower_diff_hydro     = lower_diff_hydro,
    upper_diff_hydro     = upper_diff_hydro,
    prop_rest_med        = prop_rest_med,
    lower_prop_rest      = lower_prop_rest,
    upper_prop_rest      = upper_prop_rest,
    prop_hydro_med       = prop_hydro_med,
    lower_prop_hydro     = lower_prop_hydro,
    upper_prop_hydro     = upper_prop_hydro,
    ratio_rest_med       = ratio_rest_med,   # KEY: scenario/baseline ratio
    lower_ratio_rest     = lower_ratio_rest,
    upper_ratio_rest     = upper_ratio_rest,
    ratio_hydro_med      = ratio_hydro_med,
    lower_ratio_hydro    = lower_ratio_hydro,
    upper_ratio_hydro    = upper_ratio_hydro
  )
  
  return(age0_dat)
}


# =============================================================================
# SECTION 2: LOAD AND PREPARE INPUT DATA
# =============================================================================

# Ordered reach name vector. Order determines reach_num indexing (1, 2, 3)
# used throughout preds_calc() and maxmins_fun(). Do not reorder without
# updating all downstream index references.
reaches <- c("San Acacia", "Isleta", "Angostura")

# USGS daily mean discharge records, 1993-2020.
# Albuquerque gage (08330000): used for Angostura and as Isleta proxy.
# San Acacia gage (08354900): used for San Acacia reach.
angQ  <- read.csv("data/yackulic2022_data/abq_gage_08330000.csv") %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))
sanaQ <- read.csv("data/yackulic2022_data/SanAcacia_gage_08354900.csv") %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

# Median daily flow by calendar day across all years (for hydrograph plots)
aQ_med_dly <- angQ %>%
  filter(month > 2, month < 10) %>%
  group_by(month, day) %>%
  summarise(med_cfs = median(cfs))

sQ_med_dly <- sanaQ %>%
  filter(month > 2, month < 10) %>%
  group_by(month, day) %>%
  summarise(med_cfs = median(cfs))

# Annual Apr-Jun total volume by year. tot_vol is used by q_sim() to match
# quantile inputs to historical years. cfs * 86400 sec/day = cubic feet/day.
aQ_med_spr <- angQ %>%
  mutate(vol = cfs * 86400) %>%
  filter(month > 3, month < 7) %>%
  group_by(year) %>%
  summarize(med_cfs = median(cfs), tot_vol = sum(vol))

sQ_med_spr <- sanaQ %>%
  mutate(vol = cfs * 86400) %>%
  filter(month > 3, month < 7) %>%
  group_by(year) %>%
  summarize(med_cfs = median(cfs), tot_vol = sum(vol))

# Historical LCC covariates from the model fitting period.
# cov_list[[1]] selects values for inund_lower (matching hab_cov below).
# CRITICAL: mean(cov_list) is subtracted from preds_array to center predictions
# on the same scale as the fitted model. Must be updated if model is re-fitted.
inund_lcc_list_combined <- readRDS("output/inund_lcc_list_combined.RData")

# Raw 2D hydraulic model output at 7 discharge breakpoints per reach.
# Pivoted to long format so a single hab_type can be selected by name.
q_inund_curve <- read_csv("data/ecoval2d.csv") %>%
  select(cfs = q,
         inund_lower = tot_outside_chan_lower,  # conservative channel boundary
         inund_upper = tot_outside_chan_upper,
         wua_lower, wua_upper,
         reach, reach_num) %>%
  pivot_longer(cols = inund_lower:wua_upper,
               names_to  = "hab_type",
               values_to = "hab")


# =============================================================================
# FIGURE 2: Baseline flow-inundation curves
# =============================================================================
# ColorBrewer Dark2 palette: https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
pal_df <- data.frame(reach = reaches,
                     color = c("#7570b3", "#1b9e77", "#d95f02"),
                     ltype = c(4, 1, 3))
cols  <- setNames(pal_df$color, pal_df$reach)
ltype <- setNames(pal_df$ltype, pal_df$reach)

# cfs * 0.0283168466 converts to m3/s; acres / 247.1 converts to km2
q_inund_curve %>%
  filter(hab_type == "inund_lower") %>%
  ggplot() +
  geom_line(aes(cfs * 0.0283168466, hab / 247.1,
                color = reach, linetype = reach), size = 1.5) +
  scale_color_manual(values = cols) +
  scale_linetype_manual(values = ltype) +
  labs(x = expression("Discharge (m"^{3}*"/s)"),
       y = expression("Floodplain inundation (km"^{2}*")")) +
  theme_ipsum() +
  theme(plot.margin = unit(c(1, 0, 1, 1), "cm"),
        axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "lines"))

#ggsave(filename = "plots/fig2_ecovalue_curve.jpeg", bg = "white", width = 8, height = 6, units = "in")


# --- Model configuration ------------------------------------------------------
# spawn_rank: 1 = minimum effS year, 9 = median, 17 = maximum across 17 years.
# Must be updated (e.g., to 14 for median of 28 years) if dataset is extended.
spawn_rank <- 9

# mod_name: selects which fitted model to load. "lcc_low_combined" uses the
# inund_lower habitat metric fit to all three reaches simultaneously.
mod_name <- "lcc_low_combined"

# hab_cov: habitat metric column name; must match the fitted model.
# "inund_lower" = floodplain area outside the active channel (conservative boundary)
hab_cov <- "inund_lower"

# Qdur: rolling window width (days) in maxmins_fun() and calc_cov_combined().
# Value of 20 comes from expert elicitation question type "4" (expert 3).
Qdur <- 20

# MCMC posterior samples: 15,000 rows subsampled from the full chain.
# Columns: a, mu_lbeta.1-.3, sd_lbeta, B_lbeta, effS.* (17 yrs x 3 reaches)
sim_mod <- read_csv(paste0("output/", mod_name, "_mcmc_sub.csv"))

# Index [[1]] selects inund_lower covariate values for centering predictions
cov_list <- inund_lcc_list_combined[[1]]

# Reach number lookup: required for joins and array indexing in preds_calc()
reach_df <- data.frame(reach = reaches, reach_num = c(1, 2, 3))

# Zero-habitat anchor: ensures the lookup is defined at 0 cfs (no inundation).
# Required to prevent lin_int() from extrapolating below the lowest breakpoint.
no_flow_hab <- reach_df %>% mutate(hab = 0, cfs = 0)

# Load interpolated lookup table (generated from ecoval2d.csv breakpoints)
q_hab_lookup_combined <- read_csv("data/q_hab_lookup_combined.csv")

# Baseline flow-habitat lookup at 1-cfs resolution for all three reaches
hab_lookup_baseline <- q_hab_lookup_combined %>%
  select(cfs, reach, hab = any_of(hab_cov)) %>%
  left_join(reach_df) %>%
  bind_rows(no_flow_hab) %>%
  arrange(reach, cfs)


# =============================================================================
# SECTION 3: RESTORATION SCENARIO
#
# Estimates recruitment if all identified San Acacia restoration parcels were
# implemented simultaneously. Key steps:
#   1. Load per-parcel hydraulic data (san_acacia_rest_ecovals.csv)
#   2. Compute each parcel's incremental habitat above existing baseline
#   3. Sum across parcels; combine with baseline for full reach lookup
#   4. Run preds_calc() and forecast_recruit() as for any other scenario
#
# NOTE: the deployed version uses lookup_i_2 (no channel habitat addition).
# lookup_i (with channel habitat) is retained for comparison but not used.
# See earlier discussion and certification document Section 3 for details.
# =============================================================================

# Raw per-parcel hydraulic data. Filter to above-bank terrain only ("AB"):
# channel inundation is excluded to avoid double-counting with the baseline.
q_hab_rest_raw <- read_csv("data/san_acacia_rest_ecovals.csv") %>%
  filter(terrain == "AB") %>%
  select(-terr_inund_notes, raw_hab = hab)

# Figure S1: site-level q-hab curves in metric units
ggplot(q_hab_rest_raw) +
  geom_line(aes(cfs * 0.0283168466, raw_hab * 4046.86 / 1000,
                color = as.factor(site))) +
  scale_color_colorblind() +
  labs(x = expression("Discharge (m"^{3}*"/s)"),
       y = expression("Floodplain inundation (in thousands of m"^{2}*")")) +
  theme_ipsum() +
  theme(plot.margin = unit(c(1, 0, 1, 1), "cm"),
        axis.title.x = element_text(size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 16, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.width = unit(2, "lines"))

#ggsave(filename = "plots/figS1_rest_q_hab.jpeg", bg = "white", width = 8, height = 6, units = "in")

# Total San Acacia inundation at 7,000 cfs = 6,284 acres (from 2D hydraulic model).
# Used as denominator to normalize each parcel's contribution to reach-level habitat.
q_inund_curve %>%
  filter(hab_type == "inund_lower", reach == "San Acacia", cfs == 7000)

# Compute each parcel's proportional size relative to full reach at 7,000 cfs
sana_site_ratio <- q_hab_rest_raw %>%
  filter(cfs == 7000) %>%
  select(site, size, raw_hab) %>%
  mutate(ratio = raw_hab / 6284)  # 6284 acres = total San Acacia at 7000 cfs

# Zero-habitat placeholders for non-restoration reaches (required by hab_lookup_fun)
no_rest_reaches       <- reaches[c(2, 3)]
q_hab_curves_non_rest <- tibble(
  cfs   = rep(c(0, 7000), length(no_rest_reaches)),
  reach = rep(no_rest_reaches, each = 2),
  hab   = 0
)

hab_lookup_rest_all <- tibble()
rest_sites          <- unique(q_hab_rest_raw$site)

# Discharge at which each parcel first inundates (used in lookup_i only)
init_inund <- q_hab_rest_raw %>%
  select(site, cfs, raw_hab) %>%
  filter(cfs != 0) %>%
  group_by(site) %>%
  slice_min(raw_hab) %>%
  slice_max(cfs) %>%
  select(site, cfs)

restReach <- "San Acacia"

for (i in 1:length(rest_sites)) {
  
  # Proportional size of this parcel relative to the full reach at 7,000 cfs
  reach_prop_i <- filter(sana_site_ratio, site == rest_sites[i]) %>% pull(ratio)
  init_inund_i <- filter(init_inund, site == rest_sites[i]) %>% pull()
  
  # Scale baseline reach habitat to this parcel's proportional footprint
  baseline_hab_prop_i <- hab_lookup_baseline %>%
    mutate(hab_site_base = hab * reach_prop_i)
  
  # Baseline channel habitat at the parcel's initial inundation flow (lookup_i only)
  hab_channel_max <- baseline_hab_prop_i %>%
    filter(reach == restReach, cfs == init_inund_i) %>%
    pull(hab_site_base)
  
  q_hab_i  <- filter(q_hab_rest_raw, site == rest_sites[i]) %>%
    select(cfs, reach, hab = raw_hab, site)
  curves_i <- bind_rows(q_hab_curves_non_rest, q_hab_i)
  
  # lookup_i: ALTERNATIVE -- adds channel habitat contribution above init_inund_i.
  # Retained for comparison; NOT used in the deployed model.
  lookup_i <- hab_lookup_fun(curves_i, reaches) %>%
    mutate(site = rest_sites[i]) %>%
    rename(raw_hab = hab) %>%
    left_join(baseline_hab_prop_i) %>%
    mutate(rest_hab = case_when(
      cfs > init_inund_i ~ raw_hab + hab_channel_max,
      TRUE               ~ hab_site_base)) %>%
    mutate(rest_hab_norm = rest_hab - hab_site_base)
  
  # lookup_i_2: DEPLOYED -- incremental restoration habitat only (no channel add).
  # Negative values clamped to 0 (restoration never reduces habitat).
  lookup_i_2 <- hab_lookup_fun(curves_i, reaches) %>%
    mutate(site = rest_sites[i]) %>%
    rename(raw_hab = hab) %>%
    left_join(baseline_hab_prop_i) %>%
    mutate(rest_hab_norm = raw_hab - hab_site_base) %>%
    mutate(rest_hab_norm = case_when(
      rest_hab_norm < 0 ~ 0,
      TRUE              ~ rest_hab_norm))
  
  hab_lookup_rest_all <- bind_rows(hab_lookup_rest_all, lookup_i_2)
}

# Sum incremental restoration habitat across all parcels within each reach/cfs.
# Non-restoration reaches forced to 0. hab_lookup_rest = incremental only.
hab_lookup_rest <- hab_lookup_rest_all %>%
  group_by(cfs, reach) %>%
  summarize(hab = sum(rest_hab_norm)) %>%
  left_join(reach_df) %>%
  bind_rows(no_flow_hab) %>%
  arrange(reach, cfs) %>%
  mutate(hab = case_when(!reach %in% restReach ~ 0, TRUE ~ hab))

# Combined lookup: baseline + restoration incremental. Used in preds_calc()
# as the "with restoration" habitat input for lcc_rest.
hab_lookup_rest_total <- bind_rows(hab_lookup_baseline, hab_lookup_rest) %>%
  group_by(reach, cfs, reach_num) %>%
  summarize(hab = sum(hab))

# Run restoration scenario across all historical years
years      <- c(1993:2020)
simQ       <- q_sim(years)
start_date <- "05-20"
end_date   <- "06-14"
augInput   <- 0  # set to 0: null augmentation required by function signatures
augQ       <- aug_q_fun(start_date, end_date, augInput, simQ)

ee         <- readRDS("output/ee_list.RData")
expert_num <- 3  # selects third expert's elicited parameter values

# Compute LCC covariates for each scenario
lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline)
lcc_rest     <- preds_calc(expert_num, simQ, hab_lookup_rest_total)
lcc_hydro    <- preds_calc(expert_num, augQ, hab_lookup_baseline)  # same as baseline (augInput=0)

# Flow index and inundation metrics for output.
# maxmins_rest uses incremental lookup (not combined) so acres_inund_rest
# reports restoration-added habitat alone.
maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
maxmins_rest     <- maxmins_fun(hab_lookup_rest,     simQ, reaches)
maxmins_hydro    <- maxmins_fun(hab_lookup_baseline, augQ, reaches)

# Center LCC covariates on the historical fitting distribution.
# Subtracting mean(cov_list) ensures predictions are on the same scale as
# the fitted Beverton-Holt parameters.
preds_array        <- array(NA, dim = c(length(years), length(reaches), 3))
preds_array[, , 1] <- lcc_baseline - mean(cov_list)  # baseline
preds_array[, , 2] <- lcc_rest     - mean(cov_list)  # restoration
preds_array[, , 3] <- lcc_hydro    - mean(cov_list)  # flow aug (null)

age0_dat <- forecast_recruit(sim_mod, preds_array, years,
                             maxminsBase = maxmins_baseline,
                             maxminsAug  = maxmins_hydro,
                             maxminsRest = maxmins_rest,
                             reaches)

# Retain San Acacia only and tag with scenario label for combined figure
age0_dat_rest <- age0_dat %>%
  filter(reach == "San Acacia") %>%
  mutate(scenario = "rest_sites")


# =============================================================================
# SECTION 4: BERM SCENARIO
#
# Estimates recruitment benefit of temporary flow-deflection berms in San Acacia.
# Berms raise local water surface, shifting the flow-inundation curve left by
# approximately 750 cfs (aggradation effect).
#
# Key parameters:
#   riv_slope = 0.0007 ft/ft: literature value for longitudinal slope
#   berm_height = 1 ft: height of each berm structure
#   n_berms = 5: number of berms installed
#   Washout threshold: 5,000 cfs -- berms assumed to fail above this flow
#
# Critical difference from other scenarios: uses calc_cov_combined_berm() and
# preds_calc_berm(), which accept a pre-computed raw habitat array rather than
# a lookup table. This is necessary because the flow-dependent washout rule
# cannot be represented in a static lookup table -- it depends on the daily
# flow sequence in each simulated year.
# =============================================================================

# Load combined inundation curve for berm hydraulics
q_inund_curve_combined <- read_csv("data/q_inund_combined.csv") %>%
  select(cfs = q, inund_lower, inund_upper, wua_lower, wua_upper, reach) %>%
  pivot_longer(cols = inund_lower:wua_upper,
               names_to  = "hab_type",
               values_to = "hab") %>%
  mutate(reach_num = case_when(reach == "San Acacia" ~ 1,
                               reach == "Isleta"     ~ 2,
                               reach == "Angostura"  ~ 3))

# Apply aggradation shift: habitat requiring Q cfs in natural channel now
# occurs at Q-750 cfs with berms installed. Values below 0 clamped to 0.
q_inund_aggdeg <- q_inund_curve_combined %>%
  mutate(agg_cfs = case_when(cfs - 750 < 0 ~ 0, TRUE ~ cfs - 750))

# Berm geometry parameters
riv_slope       <- 0.0007  # ft/ft longitudinal slope
berm_height     <- 1       # ft per berm
n_berms         <- 5       # number of berms
# Length of reach affected = n_berms * (berm_height / slope), converted ft -> m
length_affected <- 5 * berm_height / riv_slope / 3.28084

# San Acacia reach geometry from shapefile (needed to compute prop_inund)
sana_reach <- sf::st_read("data/ReachSegments/ReachSegments.shp") %>%
  filter(ReachName == "San Acacia Reach") %>%
  sf::st_transform(crs = 4326)

# Fraction of total reach length affected by berms
prop_inund <- length_affected / st_length(sana_reach) %>% as.numeric()

# Aggradation-shifted habitat scaled to the affected reach proportion
q_agg_initial <- q_inund_aggdeg %>%
  filter(hab_type == hab_cov) %>%
  select(reach, hab, cfs = agg_cfs) %>%
  mutate(hab = prop_inund * hab)

# Extend to 10,000 cfs with flat habitat (max value) to cover lookup range
q_agg <- q_agg_initial %>%
  group_by(reach) %>%
  summarize(hab = max(hab)) %>%
  mutate(cfs = 10000) %>%
  bind_rows(q_agg_initial) %>%
  arrange(reach, cfs)

# Baseline habitat for the berm-affected portion of the reach
hab_lookup_baseline_inund_prop <- mutate(hab_lookup_baseline,
                                         base_hab = prop_inund * hab) %>%
  select(-hab)

# Total habitat under berm scenario (from aggradation-shifted curve)
hab_lookup_berm_tot <- hab_lookup_fun(q_agg, reaches) %>%
  rename(tot_hab = hab)

# Incremental berm habitat above baseline, clamped to 0 where berm adds nothing
hab_lookup_berm <- left_join(hab_lookup_baseline_inund_prop, hab_lookup_berm_tot) %>%
  mutate(berm_hab = case_when(
    base_hab > tot_hab ~ 0,  # berm provides no benefit if baseline already higher
    is.na(tot_hab)     ~ 0,  # 0 outside interpolation range
    TRUE               ~ tot_hab - base_hab))

# Combined lookup: baseline + incremental berm habitat for both column names
hab_lookup_berm_with_baseline <- hab_lookup_berm %>%
  select(cfs, reach, hab = berm_hab, reach_num) %>%
  bind_rows(hab_lookup_baseline) %>%
  group_by(reach, cfs, reach_num) %>%
  summarize(berm_hab = sum(hab)) %>%
  left_join(hab_lookup_baseline)

simQ <- q_sim(years)
augQ <- aug_q_fun(start_date, end_date, augInput, simQ)

# --- Build raw habitat array with flow-dependent berm washout ----------------
# For each reach and year, compute daily habitat incorporating the washout rule:
# once cumulative days with Q > 5,000 cfs >= 1, berms are failed for all
# subsequent days and habitat reverts to the baseline level.
maxmins_hab_temp <- array(NA, dim = c(dim(simQ)[2], 3))
maxmins_cfs_temp <- array(NA, dim = c(dim(simQ)[2], 3))

# raw_hab_array [214 days x n_years x 3 reaches]: full daily habitat time series
raw_hab_array <- array(NA, dim = c(dim(simQ)[1], dim(simQ)[2], 3))

for (i in 1:length(reaches)) {
  hab_lookup_berm_i <- filter(hab_lookup_berm_with_baseline, reach == reaches[i])
  
  if (i == 1) {
    q_temp <- round(simQ[, , 2], 0) %>% as_tibble() %>% set_names("cfs")
  } else {
    q_temp <- round(simQ[, , 1], 0) %>% as_tibble() %>% set_names("cfs")
  }
  
  for (t in 1:dim(simQ)[2]) {
    x_temp <- q_temp[, t]
    
    hab_vec <- as.data.frame(x_temp) %>%
      left_join(hab_lookup_berm_i) %>%
      mutate(
        # cumsum(cfs > 5000) is 0 until first day Q > 5000, then stays >= 1.
        # This implements the one-time washout: once failed, stays failed.
        below_5k = cumsum(cfs > 5000),
        hab_i    = if_else(below_5k >= 1, hab, berm_hab)  # revert to baseline after washout
      ) %>%
      pull(hab_i)
    
    raw_hab_array[, t, i] <- hab_vec
    
    maxmins_hab_temp[t, i] <- max(rollapplyr(hab_vec, Qdur, mean, fill = NA)[32:122], na.rm = TRUE)
    maxmins_cfs_temp[t, i] <- max(rollapplyr(x_temp, Qdur, mean, fill = NA)[32:122], na.rm = TRUE)
  }
}

# --- Berm-specific covariate and prediction functions ------------------------
# These are identical to calc_cov_combined() and preds_calc() except they
# accept the pre-computed daily habitat vector directly, bypassing the lookup
# table step. This allows the flow-dependent washout to propagate into the LCC
# covariate calculation.

calc_cov_combined_berm <- function(Q, raw_hab_reach_yr, kappa, D, prPARS) {
  # raw_hab_reach_yr is a 214-element vector of daily habitat already
  # incorporating the washout rule; log-transform directly.
  thab  <- log(raw_hab_reach_yr + 1)
  tprop <- prop(prPARS)
  
  t2hab <- numeric()
  for (i in 1:132) {
    t2hab[i] <- min(thab[i:(i + D)])
  }
  
  tstart <- min(c(which(Q > kappa),
                  which(Q[-1] - Q[-length(Q)] > 100),
                  132))
  
  out <- sum(tprop[1:tstart]) * t2hab[tstart] +
    sum(tprop[tstart:132] * t2hab[tstart:132])
  return(out)
}

preds_calc_berm <- function(exp, flows, hab_array) {
  preds <- array(NA, dim = c(dim(flows)[2], 3))
  
  # Expert elicitation reconstruction: identical to preds_calc()
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "1" & is.na(ee[[exp]][, 5]) == FALSE)
  t1 <- numeric()
  for (t in 1:6) {
    se         <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t1[t]      <- temp[t, 6] - se
    se         <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100)
    t1[(t + 7)] <- temp[(t + 6), 6] + se
  }
  t1[7]  <- 1
  t1[14] <- 0
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "2")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t2   <- temp[1, 6]
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "4")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t4   <- temp[1, 6] - se
  
  for (r in 1:3) {
    raw_hab_reach <- hab_array[, , r]  # 214 x n_years habitat matrix for reach r
    for (j in 1:dim(flows)[2]) {
      t3 <- raw_hab_reach[, j]         # 214-element daily habitat vector for year j
      if (r == 3) {q <- flows[, j, 1]} else {q <- flows[, j, 2]}
      preds[j, r] <- calc_cov_combined_berm(q, t3, t2, t4, t1)
    }
  }
  return(preds)
}

# Run berm scenario
lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline)
lcc_rest     <- preds_calc_berm(expert_num, simQ, raw_hab_array)  # berm-specific
lcc_hydro    <- preds_calc(expert_num, augQ, hab_lookup_baseline)

maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
# NOTE: maxmins_rest for berms uses the pre-computed arrays from the raw_hab
# loop, NOT maxmins_fun(), because the washout rule cannot be represented
# in a static lookup table.
maxmins_rest     <- list(maxmins_hab = maxmins_hab_temp, maxmins_cfs = maxmins_cfs_temp)
maxmins_hydro    <- maxmins_fun(hab_lookup_baseline, augQ, reaches)

preds_array        <- array(NA, dim = c(length(years), length(reaches), 3))
preds_array[, , 1] <- lcc_baseline - mean(cov_list)
preds_array[, , 2] <- lcc_rest     - mean(cov_list)
preds_array[, , 3] <- lcc_hydro    - mean(cov_list)

age0_dat_berm <- forecast_recruit(sim_mod, preds_array, years,
                                  maxminsBase = maxmins_baseline,
                                  maxminsAug  = maxmins_hydro,
                                  maxminsRest = maxmins_rest,
                                  reaches)

age0_dat_berm_scen <- age0_dat_berm %>%
  filter(reach == "San Acacia") %>%
  mutate(scenario = "berm")


# =============================================================================
# SECTION 5: FLOW AUGMENTATION SCENARIO
#
# Estimates recruitment benefit of adding 10,000 acre-feet as an environmental
# flow over a 27-day window, tested across all May start dates (May 1-31).
# 10,000 ac-ft over 27 days ~= 170 cfs constant addition.
#
# Purpose: identify which timing produces the most consistent benefit across
# the full range of historical flow years.
# =============================================================================

years    <- c(1993:2020)
simQ     <- q_sim(years)
augInput <- 10000  # acre-feet

# Generate all 31 possible 27-day windows starting in May
dates <- data.frame(start = seq(as.Date('2011-05-01'), as.Date('2011-05-31'), by = 1)) %>%
  mutate(end      = start + 26,
         start_dm = format(start, "%m-%d"),
         end_dm   = format(end,   "%m-%d"))

aug_recruit_all <- data.frame()

# Loop over all start dates: run full model pipeline for each
for (i in 1:nrow(dates)) {
  augQ <- aug_q_fun(dates$start_dm[i], dates$end_dm[i], augInput, simQ)
  
  lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline)
  lcc_rest     <- preds_calc(expert_num, simQ, hab_lookup_rest_total)
  lcc_hydro    <- preds_calc(expert_num, augQ, hab_lookup_baseline)
  
  maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
  maxmins_rest     <- maxmins_fun(hab_lookup_rest,     simQ, reaches)
  maxmins_hydro    <- maxmins_fun(hab_lookup_baseline, augQ, reaches)
  
  preds_array        <- array(NA, dim = c(length(years), length(reaches), 3))
  preds_array[, , 1] <- lcc_baseline - mean(cov_list)
  preds_array[, , 2] <- lcc_rest     - mean(cov_list)
  preds_array[, , 3] <- lcc_hydro    - mean(cov_list)
  
  age0_dat <- forecast_recruit(sim_mod, preds_array, years,
                               maxminsBase = maxmins_baseline,
                               maxminsAug  = maxmins_hydro,
                               maxminsRest = maxmins_rest,
                               reaches) %>%
    mutate(start_date = dates$start_dm[i])  # tag with start date
  
  aug_recruit_all <- bind_rows(aug_recruit_all, age0_dat)
}

# Exploratory boxplot: prop_hydro_med vs flow index across date windows
filter(aug_recruit_all, reach == "San Acacia") %>%
  ggplot(aes(as.factor(cfs_baseline), y = prop_hydro_med)) +
  geom_boxplot()

# Figure S2: recruitment ratio vs flow index colored by start date
aug_recruit_all %>%
  mutate(
    # Separate anomalous years (2013, 2020 had unusually high ratios)
    year_set       = case_when(flow_years %in% c(2013, 2020) ~ "2013 and 2020",
                               TRUE                          ~ "All other years"),
    start_date_num = as.numeric(str_sub(start_date, start = -2, end = -1))
  ) %>%
  filter(reach == "San Acacia") %>%
  ggplot(aes(cfs_baseline * 0.0283168466, y = ratio_hydro_med,
             color = start_date_num)) +
  geom_point() +
  facet_wrap(~ year_set, scales = "free_y", nrow = 2) +
  labs(y     = "Proportion of\nbaseline recruitment",
       x     = expression("Spring flow index (m"^{3}*"/s)"),
       color = "May\nstarting\ndate") +
  theme_ipsum() +
  theme(plot.margin  = unit(c(1, 0, 1, 1), "cm"),
        strip.text   = element_text(size = 14, hjust = 0.5),
        axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14))

#ggsave(filename = "plots/figS2_flow_aug_dates.jpeg", bg = "white", width = 6, height = 6, units = "in")

# Rank each year by prop_hydro_med within start_date to assess consistency
aug_recruit_all %>%
  filter(reach == "San Acacia") %>%
  group_by(cfs_baseline) %>%
  mutate(rank = rank(prop_hydro_med)) %>%
  ggplot(aes(start_date, rank)) +
  geom_boxplot()

# Find the start date with the median effect on recruitment.
# For each date, compute the median rank across years; slice(16) selects
# the 16th of 31 sorted entries = the median start date.
med_aug_recruit <- aug_recruit_all %>%
  filter(reach == "San Acacia") %>%
  group_by(cfs_baseline) %>%
  mutate(rank = rank(prop_hydro_med)) %>%
  group_by(start_date) %>%
  summarize(med_rank = median(rank)) %>%
  arrange(med_rank) %>%
  slice(16)

med_aug_recruit$start_date  # print selected median date

# Prepare augmentation data for combined figure: rename hydro columns to "rest"
# naming convention so all three scenarios can share the same plotting code.
aug_examp <- aug_recruit_all %>%
  filter(start_date == med_aug_recruit$start_date, reach == "San Acacia") %>%
  mutate(scenario = "Aug") %>%
  select(-c(recruits_rest, diff_recruits_rest, prop_rest_med,
            lower_prop_rest, upper_prop_rest, ratio_rest_med,
            lower_ratio_rest, upper_ratio_rest)) %>%
  rename(recruits_rest      = recruits_hydro,
         diff_recruits_rest = diff_recruits_hydro,
         prop_rest_med      = prop_hydro_med,
         lower_prop_rest    = lower_prop_hydro,
         upper_prop_rest    = upper_prop_hydro,
         ratio_rest_med     = ratio_hydro_med,
         lower_ratio_rest   = lower_ratio_hydro,
         upper_ratio_rest   = upper_ratio_hydro)


# =============================================================================
# SECTION 6: COMBINED SCENARIO FIGURES (Fig. 4 / Table S1)
# =============================================================================

# Combine all three scenario data frames for joint plotting
berm_rest_hydro <- bind_rows(age0_dat_rest, age0_dat_berm_scen, aug_examp) %>%
  mutate(Scenario = case_when(
    scenario == "rest_sites" ~ "Floodplain\nRestoration",
    scenario == "Aug"        ~ "Flow\nAugmentation",
    TRUE                     ~ "Temporary\nBerms"))

# Colorblind-safe palette (Wong 2011 Nature Methods):
# #56B4E9 sky blue, #D55E00 vermillion, #009E73 green
pal_df2 <- data.frame(Scenario = unique(berm_rest_hydro$Scenario),
                      color    = c("#56B4E9", "#D55E00", "#009E73"))
cols2   <- setNames(pal_df2$color, pal_df2$Scenario)

# Panel 1: Raw recruitment -- error bars span baseline to scenario
# Vertical dashed lines at 28.3 and 113 m3/s mark inundation onset / saturation
raw_plot <- ggplot(berm_rest_hydro,
                   aes(cfs_baseline * 0.0283168466,
                       recruits_rest / 1000000, color = Scenario)) +
  geom_errorbar(aes(ymin = recruits_baseline / 1000000,
                    ymax = recruits_rest      / 1000000), color = "gray50") +
  geom_point(aes(y = recruits_baseline / 1000000), color = "gray50") +
  geom_point() +
  geom_vline(aes(xintercept = 28.3),  linetype = 2) +
  geom_vline(aes(xintercept = 113.0), linetype = 2) +
  scale_color_manual(values = cols2) +
  labs(y = "Recruitment\n(millions)", x = "Flow index") +
  theme_ipsum() +
  facet_wrap(~ Scenario) +
  theme(plot.margin     = unit(c(1, 0, 1, 1), "cm"),
        legend.position = "none",
        axis.title.x    = element_blank(),
        axis.text.x     = element_blank(),
        strip.text      = element_text(size = 14, hjust = 0.5),
        axis.title.y    = element_text(size = 12, hjust = 0.5))

raw_plot

# Panel 2: Recruits above baseline (difference)
diff_plot <- ggplot(berm_rest_hydro,
                    aes(cfs_baseline * 0.0283168466,
                        y = diff_recruits_rest / 1000000, color = Scenario)) +
  geom_point() +
  scale_color_manual(values = cols2) +
  geom_vline(aes(xintercept = 28.3),  linetype = 2) +
  geom_vline(aes(xintercept = 113.0), linetype = 2) +
  labs(y = "Recruitment above\nbaseline (millions)", x = "Flow index") +
  theme_ipsum() +
  facet_wrap(~ Scenario) +
  theme(plot.margin     = unit(c(1, 0, 1, 1), "cm"),
        legend.position = "none",
        axis.title.x    = element_blank(),
        strip.text      = element_blank(),
        axis.text.x     = element_blank(),
        axis.title.y    = element_text(size = 12, hjust = 0.5))

diff_plot

# Panel 3: Recruitment ratio (scenario / baseline) with 80% credible intervals.
# y-axis capped at 2.72 (observed maximum across all scenarios and years).
prop_plot <- ggplot(berm_rest_hydro,
                    aes(cfs_baseline * 0.0283168466,
                        ratio_rest_med, color = Scenario)) +
  geom_point() +
  geom_vline(aes(xintercept = 28.3),  linetype = 2) +
  geom_vline(aes(xintercept = 113.0), linetype = 2) +
  geom_errorbar(aes(ymin = lower_ratio_rest, ymax = upper_ratio_rest)) +
  scale_y_continuous(limits = c(1, 2.72)) +
  scale_color_manual(values = cols2) +
  labs(y = "Proportion of\nbaseline recruitment",
       x = expression("Spring flow index (m"^{3}*"/s)")) +
  theme_ipsum() +
  facet_wrap(~ Scenario) +
  theme(plot.margin     = unit(c(1, 0, 1, 1), "cm"),
        legend.position = "none",
        strip.text      = element_blank(),
        axis.title.x    = element_text(size = 12, hjust = 0.5),
        axis.title.y    = element_text(size = 12, hjust = 0.5))

prop_plot

# Composite Figure 4: stack panels with negative patchwork heights to remove
# vertical gaps between panels sharing the same x-axis range
raw_plot / plot_spacer() / diff_plot / plot_spacer() / prop_plot +
  plot_layout(heights = c(4.5, -2.1, 4.5, -2.1, 4.5))

#ggsave(filename = "plots/fig4_scenario_plot_v4_no_channel.jpeg", bg = "white", width = 6, height = 6, units = "in")

# Table S1: formatted summary table in metric units
output <- berm_rest_hydro %>%
  mutate(Discharge = cfs_baseline * 0.0283168466) %>%
  select(Scenario,
         Year                       = flow_years,
         Discharge,
         `Recruits\n(baseline)`     = recruits_baseline,
         `Recruits\n(scenario)`     = recruits_rest,
         `Recruits\nover\nbaseline` = diff_recruits_rest,
         `Recruit\nratio`           = ratio_rest_med) %>%
  mutate(across(`Recruits\n(baseline)`:`Recruits\nover\nbaseline`,
                \(x) round(x, digits = 0)),
         `Recruit\nratio` = round(`Recruit\nratio`, 2),
         Discharge        = round(Discharge, 1))

#write_csv(output, "output/supp_tab_m3_new_v2.csv")


# =============================================================================
# SECTION 7: SITE-LEVEL RECRUITMENT ANALYSIS (Figures 5 and S3)
#
# Estimates the recruitment benefit of each individual restoration parcel
# separately, allowing comparison across sites by area, type, and flow regime.
# Uses the same model pipeline as Section 3 but looped once per site.
# =============================================================================

years      <- c(1993:2020)
simQ       <- q_sim(years)
start_date <- "05-20"
end_date   <- "06-14"
augInput   <- 0
augQ       <- aug_q_fun(start_date, end_date, augInput, simQ)
restReach  <- "San Acacia"

age0_dat_all <- data.frame()

for (i in 1:length(rest_sites)) {
  
  # Combined lookup for this site only: incremental site habitat + baseline
  hab_lookup_rest_i <- hab_lookup_rest_all %>%
    filter(site == rest_sites[i]) %>%
    select(cfs, reach, hab = rest_hab_norm) %>%
    left_join(reach_df) %>%
    bind_rows(no_flow_hab) %>%
    arrange(reach, cfs) %>%
    mutate(hab = case_when(!reach %in% restReach ~ 0, TRUE ~ hab)) %>%
    bind_rows(hab_lookup_baseline) %>%
    group_by(reach, cfs, reach_num) %>%
    summarize(hab = sum(hab))
  
  lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline)
  lcc_rest     <- preds_calc(expert_num, simQ, hab_lookup_rest_i)
  lcc_hydro    <- preds_calc(expert_num, augQ, hab_lookup_baseline)
  
  maxmins_baseline <- maxmins_fun(hab_lookup_baseline,  simQ, reaches)
  maxmins_rest     <- maxmins_fun(hab_lookup_rest_i,    simQ, reaches)
  maxmins_hydro    <- maxmins_fun(hab_lookup_baseline,  augQ, reaches)
  
  preds_array        <- array(NA, dim = c(length(years), length(reaches), 3))
  preds_array[, , 1] <- lcc_baseline - mean(cov_list)
  preds_array[, , 2] <- lcc_rest     - mean(cov_list)
  preds_array[, , 3] <- lcc_hydro    - mean(cov_list)
  
  age0_dat <- forecast_recruit(sim_mod, preds_array, years,
                               maxminsBase = maxmins_baseline,
                               maxminsAug  = maxmins_hydro,
                               maxminsRest = maxmins_rest,
                               reaches) %>%
    mutate(site = rest_sites[i])
  
  age0_dat_all <- bind_rows(age0_dat_all, age0_dat)
}

# Join site metadata and filter to San Acacia
rest_site_info <- select(q_hab_rest_raw, site, rest_desc, rest_type,
                         n_rest_features, size) %>% distinct()

recruits <- age0_dat_all %>%
  filter(reach %in% c("San Acacia")) %>%
  left_join(rest_site_info) %>%
  mutate(site = as.factor(site))

# Exploratory: time series of recruitment gains by site across years
ggplot(recruits) +
  geom_path(aes(flow_years, diff_recruits_rest, color = site))

# Minimum baseline recruitment across all years; reference line for plots
recruits_min <- slice_min(recruits, recruits_baseline) %>% pull(recruits_baseline)

# Figure S3: Violin + jitter of site-level recruitment gains.
# Horizontal dashed line = minimum baseline recruitment (ecological threshold).
ggplot(recruits) +
  geom_violin(aes(site, diff_recruits_rest / 1000), outliers = FALSE) +
  geom_jitter(aes(site, diff_recruits_rest / 1000), width = .1) +
  geom_hline(aes(yintercept = recruits_min / 1000), lty = 2) +
  labs(x = "Restoration site", y = "Recruits over baseline (thousands)") +
  theme_ipsum() +
  theme(plot.margin  = unit(c(1, 0, 1, 1), "cm"),
        axis.title.x = element_text(size = 14, hjust = 0.5),
        axis.title.y = element_text(size = 14, hjust = 0.5))

#ggsave(filename = "plots/figS3_recruit_rest_sites_no_channel.jpeg", width = 6, height = 5, units = "in")

# Figure 5: Recruitment ratio by site, flow index, and restoration type.
# Point size = restoration area (acres); color = site; facet = restoration type.
ggplot(recruits,
       aes(cfs_baseline * 0.0283168466, ratio_rest_med, color = site)) +
  geom_point(aes(size = size)) +
  geom_errorbar(aes(ymin = lower_ratio_rest, ymax = upper_ratio_rest)) +
  theme_ipsum() +
  facet_wrap(~ rest_type) +
  scale_color_colorblind() +
  labs(x     = expression("Spring flow index (m"^{3}*"/s)"),
       y     = "Proportion of baseline recruitment",
       size  = "Number of\nrestoration\nfeatures",
       color = "Site") +
  theme(plot.margin  = unit(c(1, 0, 1, 1), "cm"),
        strip.text   = element_text(size = 14),
        axis.title.x = element_text(size = 12, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12))

#ggsave(filename = "plots/fig5_recruit_sites_ratio_no_channel.jpeg", bg = "white", width = 6, height = 5, units = "in")
