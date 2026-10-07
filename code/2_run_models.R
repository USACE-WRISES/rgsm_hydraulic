# =============================================================================
# Script: 2_run_models.R
# Purpose: Build and compare six candidate Bayesian population models for
#          predicting Rio Grande Silvery Minnow (RGSM) age-0 recruitment,
#          each using a different covariate to represent larval habitat
#          conditions. Fit models to 17 years of in-sample data (2002-2018),
#          save outputs, perform diagnostics, then generate out-of-sample
#          forecasts for 2019-2024 and compare model predictions against
#          observed RGSM catch data.
#
# Overview of workflow:
#   1. Define helper and covariate functions
#   2. Calculate six types of habitat covariates:
#        (a) Original LCC from Yackulic et al. 2022 (expert elicitation only)
#        (b) May-June mean flow (simple flow index)
#        (c) Raw inundation area from 2D hydraulic model (4 variants)
#        (d) Combined inundation + expert elicitation (4 variants)
#        (e) Hydraulic-based LCC covariate (4 variants)
#        (f) Combined hydraulic + expert elicitation LCC (4 variants)
#   3. Fit six Stan models (int3re.stan) to in-sample data using each covariate
#   4. Save model outputs, subsampled MCMC posteriors, and summary statistics
#   5. Run model diagnostics: R2, traceplots, posterior predictive checks
#   6. Load 2019-2024 USGS flow data for out-of-sample forecasting
#   7. Run forecast_oos_re() for each model; compare to observed RGSM catch
#   8. Produce Figure 3 (out-of-sample comparison) and R2 table
#
# The six models and their short names used throughout:
#   M2_1re             / "orig_lcc"            -- Original LCC (Yackulic et al.)
#   M2_2re             / "springflow"           -- May-June mean flow
#   m_inund_i_l        / "inund_low"            -- Raw inundation (lower boundary)
#   m_inund_i_l_combined / "inund_low_combined" -- Combined inundation (lower)
#   m_lcc_i_l          / "lcc_low"              -- Hydraulic LCC (lower boundary)
#   m_lcc_i_l_combined / "lcc_low_combined"     -- Combined hydraulic LCC (lower)
#
# Key input files:
#   output/input_data.RData        -- In-sample monitoring and flow data
#   output/oos_data_new.RData      -- Out-of-sample monitoring data (2019-2024)
#   data/ecoval2d.csv              -- 2D hydraulic model inundation breakpoints
#   code/int3re.stan               -- Stan model code (IPM with 3 reaches)
#   output/full/*.rds              -- Saved fitted Stan model objects
#   output/ee_list.RData           -- Expert elicitation values
#
# Key output files:
#   output/full/*_mcmc_sub.csv     -- Subsampled MCMC posteriors for each model
#   output/full/*_summ.csv         -- Summary statistics for each model
#   output/inund_*_list*.RData     -- Covariate lists for centering in script 3
#   plots/fig3_oos_results_*.jpeg  -- Out-of-sample comparison figure
#   plots/r2_mods.jpeg             -- R2 heatmap across models and reaches
# =============================================================================

library(tidyverse)
library(rstan)            # Bayesian model fitting via Stan
library(zoo)              # rollapplyr() for rolling window calculations
library(hrbrthemes)       # theme_ipsum() for plot styling
library(dataRetrieval)    # readNWISdv() to fetch USGS gage data via API


# =============================================================================
# SECTION 1: LOAD IN-SAMPLE DATA
# =============================================================================

# input_data.RData contains all data objects required to fit the Stan IPM,
# including monitoring catch data (mon0, mon1, mon01, monV, R0, R1),
# daily flow matrices (aQ [214 x 17], sQ [214 x 17]), expert elicitation
# (ee list), and all dimensional constants (Nyears, Nstrata, etc.).
# These objects are loaded directly into the global environment.
load(file = "output/input_data.RData")


# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# calcsig()
# Purpose: Recover the implied standard deviation from an expert-elicited
#   probability statement. Given a central estimate (mean), confidence interval
#   bounds (down, up), and the stated probability that the true value falls
#   within those bounds (q), finds the normal distribution SD consistent with
#   all three inputs via numerical optimization of the normal CDF.
#
# Used in all covariate calculations to reconstruct elicited parameter
# distributions from the summary statistics stored in ee_list.RData.
#
# Args:
#   up, down -- upper and lower bounds of the elicited confidence interval
#   mean     -- elicited central estimate
#   q        -- stated probability of containment (e.g., 0.9 for 90% CI)
#   INT      -- optimization search range for SD; default c(0, 1000)
#
# Returns: numeric scalar -- the recovered standard deviation
# -----------------------------------------------------------------------------
calcsig <- function(up, down, mean, q, INT = c(0, 1000)) {
  sig <- function(x) {
    abs(pnorm(up, mean, x) - pnorm(down, mean, x) - q)
  }
  optimize(sig, interval = INT)$minimum
}

# -----------------------------------------------------------------------------
# lin_int()
# Purpose: Linear interpolation between two known points. Core building block
#   for all flow-habitat relationships -- used to estimate habitat at arbitrary
#   discharge values between hydraulic model breakpoints.
#
# Args:
#   x           -- discharge at which to interpolate
#   xlow, xhi   -- discharge bounds of the interval
#   ylow, yhi   -- habitat values at xlow and xhi
#
# Returns: numeric -- interpolated habitat value at x
# -----------------------------------------------------------------------------
lin_int <- function(x, xlow, xhi, ylow, yhi) {
  (yhi - ylow) * (x - xlow) / (xhi - xlow) + ylow
}


# =============================================================================
# SECTION 3: FLOW-HABITAT INTERPOLATION FUNCTIONS
#
# Three versions of q2hab() exist, reflecting the evolution of the hydraulic
# data used across analyses. Each uses a different set of discharge breakpoints
# derived from different sources:
#   q2hab        -- original breakpoints from Yackulic et al. 2022 (expert)
#   q2hab_2d     -- breakpoints from 2D HEC-RAS hydraulic model
#   q2hab_combined -- combined set merging both expert and hydraulic breakpoints
# =============================================================================

# -----------------------------------------------------------------------------
# q2hab()
# Purpose: Original flow-habitat function from Yackulic et al. 2022.
#   Uses expert-elicited breakpoints at fine resolution at low flows (5, 50,
#   100, 150, 200, 250 cfs) and coarser resolution at high flows (1000-7000
#   cfs). Habitat is held flat above 7,000 cfs.
#
#   Includes a 0 cfs anchor (p = 0, qs = 0) so the function is defined at
#   zero flow. Unlike the 2D variants, pars here represents elicited habitat
#   values at each breakpoint rather than hydraulic model output.
#
# Args:
#   q    -- discharge in cfs
#   pars -- 16-element vector of elicited habitat values at each breakpoint
#
# Returns: numeric -- interpolated habitat value
# -----------------------------------------------------------------------------
q2hab <- function(q, pars) {
  # Prepend a zero-habitat anchor at q = 0 and extend flat above 7,000 cfs
  p  <- c(0, pars, pars[length(pars)])
  qs <- c(0, 5, 50, 100, 150, 200, 250, 1000, 1500, 2000, 2500, 3000, 4000,
          5000, 6000, 7000, 2 * 10 ^ 4)
  t1 <- findInterval(q, qs)
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}

# -----------------------------------------------------------------------------
# q2hab_2d()
# Purpose: Flow-habitat function using breakpoints from the 2D HEC-RAS
#   hydraulic model. Breakpoints are at the seven modeled discharge levels
#   (700, 2000, 3000, 4000, 5000, 6000, 7000 cfs), plus a 0 cfs anchor
#   and a 20,000 cfs ceiling to prevent extrapolation.
#
#   Unlike q2hab(), this version does NOT apply a log transform -- the raw
#   habitat values from the hydraulic model are interpolated directly.
#   The log transform is applied in the covariate functions (calc_cov_2d,
#   calc_cov_combined) after this function returns.
#
#   NOTE: In script 3_scenarios.R, the breakpoints are supplied explicitly
#   as a function argument (the "input" parameter), allowing any set of
#   discharge breakpoints to be used. Here they are hard-coded to match
#   the ecoval2d.csv breakpoints used during model fitting.
#
# Args:
#   q    -- discharge in cfs
#   pars -- 8-element vector of habitat values at each hydraulic breakpoint
#
# Returns: numeric -- interpolated habitat in acres
# -----------------------------------------------------------------------------
q2hab_2d <- function(q, pars) {
  p  <- c(pars, pars[length(pars)])  # extend flat above 7,000 cfs
  qs <- c(0, 700, 2000, 3000, 4000, 5000, 6000, 7000, 2 * 10 ^ 4)
  t1 <- findInterval(q, qs)
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}

# -----------------------------------------------------------------------------
# q2hab_combined()
# Purpose: Flow-habitat function using a combined set of breakpoints that
#   merges the low-flow expert elicitation breakpoints (5, 50, 100, 150,
#   200 cfs) with the high-flow 2D hydraulic model breakpoints (1000-7000 cfs).
#
#   This hybrid approach captures the expert-elicited relationship at low flows
#   (where habitat changes rapidly with small flow increments and is difficult
#   to model hydraulically) while using physically grounded hydraulic model
#   output at higher flows where inundation dynamics are better characterized.
#   The resulting q_inund_curve_combined table is saved as q_inund_combined.csv
#   for use in 3_scenarios.R.
#
# Args:
#   q    -- discharge in cfs
#   pars -- vector of habitat values at each combined breakpoint
#
# Returns: numeric -- interpolated habitat value
# -----------------------------------------------------------------------------
q2hab_combined <- function(q, pars) {
  p  <- c(pars, pars[length(pars)])
  qs <- c(0, 5, 50, 100, 150, 200, 1000, 1500, 2000, 3000, 4000, 5000, 6000,
          7000, 2 * 10 ^ 4)
  t1 <- findInterval(q, qs)
  lin_int(q, qs[t1], qs[(t1 + 1)], p[t1], p[(t1 + 1)])
}


# =============================================================================
# SECTION 4: LCC COVARIATE FUNCTIONS
#
# Three versions of calc_cov() compute the larval carrying capacity (LCC)
# covariate from a daily flow time series, using the three q2hab variants above.
# All three share the same LCC logic (see calc_cov_combined in 3_scenarios.R
# for full documentation); they differ only in which flow-habitat function
# they call and whether a log transform is applied inside the function.
# =============================================================================

# prop()
# Compute the daily proportional egg-laying distribution across the 132-day
# spawning window (March 1 = day 1 to July 10 = day 132), interpolating
# between expert-elicited values at 14 10-day breakpoints and normalizing
# so the vector sums to 1. Identical across all three model variants.
# See 3_scenarios.R for full documentation.
prop <- function(pars) {
  t    <- seq(2, 132, 10)
  tout <- pars[1]
  for (i in 2:131) {
    t1      <- findInterval(i, t)
    tout[i] <- lin_int(i, t[t1], t[(t1 + 1)], pars[t1], pars[(t1 + 1)])
  }
  tout[132] <- pars[14]
  tout / sum(tout)
}

# -----------------------------------------------------------------------------
# calc_cov()
# Purpose: Original LCC covariate calculation from Yackulic et al. 2022.
#   Uses q2hab() with expert-elicited breakpoints. No log transform is
#   applied to habitat -- raw habitat values are used directly. This is the
#   key difference from calc_cov_2d and calc_cov_combined, which apply
#   log(hab + 1) before computing the duration minimum.
#
# Args:
#   Q       -- numeric vector [214] of daily discharge in cfs
#   q2fPARS -- habitat parameters at elicited breakpoints (t3 in the code)
#   kappa   -- flow cue threshold (cfs) triggering spawning onset
#   D       -- duration requirement (days); Qdur = 20 in deployed model
#   prPARS  -- 14-element vector of elicited daily spawning proportions
#
# Returns: numeric scalar -- the LCC covariate value for this year/reach
# -----------------------------------------------------------------------------
calc_cov <- function(Q, q2fPARS, kappa, D, prPARS) {
  thab  <- numeric()
  for (i in 1:214) {
    # Note: NO log transform here; raw habitat values from q2hab() are used
    thab[i] <- q2hab(Q[i], pars = q2fPARS)
  }
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

# -----------------------------------------------------------------------------
# calc_cov_2d()
# Purpose: LCC covariate using 2D hydraulic model breakpoints (q2hab_2d).
#   Applies log(hab + 1) before computing the duration minimum, which
#   compresses the habitat scale and stabilizes the covariate across the
#   wide range of inundation values produced by the hydraulic model.
#
# Args: same as calc_cov() but q2fPARS contains hydraulic model habitat values
# Returns: numeric scalar -- LCC covariate value
# -----------------------------------------------------------------------------
calc_cov_2d <- function(Q, q2fPARS, kappa, D, prPARS) {
  thab <- numeric()
  for (i in 1:214) {
    thab[i] <- log(q2hab_2d(Q[i], pars = q2fPARS) + 1)  # log transform applied here
  }
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

# -----------------------------------------------------------------------------
# calc_cov_combined()
# Purpose: LCC covariate using the combined expert + hydraulic breakpoints
#   (q2hab_combined). Same log(hab + 1) transform as calc_cov_2d. This is
#   the version used by the deployed "lcc_low_combined" model (m_lcc_i_l_combined)
#   and carried forward into 3_scenarios.R.
#
# Args: same as calc_cov() but q2fPARS contains combined breakpoint values
# Returns: numeric scalar -- LCC covariate value
# -----------------------------------------------------------------------------
calc_cov_combined <- function(Q, q2fPARS, kappa, D, prPARS) {
  thab <- numeric()
  for (i in 1:214) {
    thab[i] <- log(q2hab_combined(Q[i], pars = q2fPARS) + 1)
  }
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


# =============================================================================
# SECTION 5: CALCULATE HABITAT COVARIATES FOR MODEL FITTING
#
# Six covariate types are computed for the 17 in-sample years (2002-2018):
#   preds2             -- original LCC (expert elicitation only, q2hab)
#   simpflow           -- May-June mean flow in thousands of cfs
#   inund_cov_list     -- peak Qdur-day rolling min of 2D hydraulic inundation
#   inund_cov_combined_list -- same using combined curve breakpoints
#   inund_lcc_list     -- LCC covariate using 2D hydraulic breakpoints
#   inund_lcc_list_combined -- LCC covariate using combined breakpoints
#
# All covariates are later centered (subtracting mean) before passing to Stan.
# =============================================================================

# ── 5A. Original LCC covariate (preds2) ─────────────────────────────────────
# Replicates the covariate from Yackulic et al. 2022 using expert-elicited
# habitat values at reach-specific discharge breakpoints.
#
# t1 = prPARS: 14-element spawning phenology vector (EE question type "1")
# t2 = kappa:  flow cue threshold in cfs (EE question type "2")
# t3 = q2fPARS: reach-specific habitat values at elicited breakpoints
#               (EE question type coded by riversegment: "3e","3c","3a"
#               for Angostura, Isleta, San Acacia respectively)
# t4 = D:      duration requirement in days (EE question type "4")
#
# NOTE: t3 construction here differs from t1/t2/t4. For question type "1"
# (spawning phenology), SD is subtracted from lower breakpoints and added
# to upper breakpoints. For question type "3" (reach habitat), it is
# REVERSED: SD is added to lower breakpoints (t3[t] = central + se for t=1:6)
# and subtracted from upper breakpoints (t3[t] = central - se for t=8:15).
# This reversal reflects the directional interpretation of habitat uncertainty
# in the original Yackulic et al. model.
preds2       <- array(NA, dim = c(17, 3))
riversegment <- c("3e", "3c", "3a")  # Angostura, Isleta, San Acacia

# Reconstruct t1 (spawning phenology) from EE question type "1", expert 3
temp <- subset(ee[[3]], ee[[3]][, 1] == "1" & is.na(ee[[3]][, 5]) == FALSE)
t1   <- numeric()
for (t in 1:6) {
  se      <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
  t1[t]   <- temp[t, 6] - se      # lower breakpoints: central - SD
  se      <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100)
  t1[(t + 7)] <- temp[(t + 6), 6] + se  # upper breakpoints: central + SD
}
t1[7]  <- 1  # peak relative spawning = 1 (fixed)
t1[14] <- 0  # end of window = 0 (fixed)

# Reconstruct t2 (kappa, flow cue) from EE question type "2"
temp <- subset(ee[[3]], ee[[3]][, 1] == "2")
se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
t2   <- temp[1, 6]  # central value used directly

# Reconstruct t4 (D, duration requirement) from EE question type "4"
temp <- subset(ee[[3]], ee[[3]][, 1] == "4")
se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
t4   <- temp[1, 6] - se  # conservative (lower) estimate

# Compute original LCC covariate for each reach and year
for (r in 1:3) {
  # Reconstruct t3 (reach-specific habitat) from EE question type riversegment[r]
  temp <- subset(ee[[3]], ee[[3]][, 1] == riversegment[r] & is.na(ee[[3]][, 5]) == FALSE)
  t3   <- numeric()
  for (t in 1:6) {
    se    <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t3[t] <- temp[t, 6] + se  # lower flow breakpoints: central + SD (opposite sign from t1)
  }
  for (t in 8:15) {
    se    <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t3[t] <- temp[t, 6] - se  # upper flow breakpoints: central - SD (opposite sign from t1)
  }
  t3[7] <- temp[7, 6]  # peak breakpoint uses central value directly
  
  for (j in 1:17) {
    # Reach-gage assignment: Angostura (r==3) uses ABQ gage (aQ); others use SA gage (sQ)
    if (r == 3) {q <- aQ[, j]} else {q <- sQ[, j]}
    preds2[j, r] <- calc_cov(q, t3, t2, t4, t1)
  }
}

# ── 5B. May-June mean flow covariate (simpflow) ──────────────────────────────
# Simple alternative covariate: mean daily discharge during May-June at each
# gage, expressed in thousands of cfs (dividing by 1000 to place on a similar
# scale as the LCC covariate for model comparison).
# San Acacia gage is used for both San Acacia and Isleta (columns 1 and 2);
# Albuquerque gage is used for Angostura (column 3). Years 2002-2018 only.
mjflow_sana <- subset(sanaQ, (sanaQ$month == 5 | sanaQ$month == 6) &
                        sanaQ$year > 2001 & sanaQ$year < 2019)
mjflow_ang  <- subset(angQ,  (angQ$month  == 5 | angQ$month  == 6) &
                        angQ$year  > 2001 & angQ$year  < 2019)

# tapply computes the mean cfs for each year; cbind combines gages into a
# matrix [17 years x 3 reaches] matching the shape of other covariates.
simpflow <- cbind(tapply(mjflow_sana$cfs, mjflow_sana$year, mean),  # San Acacia
                  tapply(mjflow_sana$cfs, mjflow_sana$year, mean),  # Isleta (same gage)
                  tapply(mjflow_ang$cfs,  mjflow_ang$year,  mean)) / 1000  # Angostura / 1000

# ── 5C. Raw inundation covariate (inund_cov_list) ────────────────────────────
# Computes the peak Qdur-day rolling minimum of inundated area (acres) within
# the April-June spawning window for each year and reach, using the 2D hydraulic
# model lookup table (q_hab_lookup). Four habitat variants are computed (i=1:4):
#   i=1: inund_lower (conservative channel boundary) -- used in deployed model
#   i=2: inund_upper (wider channel boundary)
#   i=3: wua_lower (weighted usable area, conservative)
#   i=4: wua_upper (weighted usable area, wider boundary)
# Only inund_cov_list[[1]] (inund_lower) is used for model fitting; all four
# are saved for potential sensitivity analysis.

# Load 2D hydraulic model breakpoints (7 discharge levels per reach)
q_inund_curve <- read_csv("data/ecoval2d.csv") %>%
  select(q,
         inund_lower = tot_outside_chan_lower,  # floodplain area, conservative boundary
         inund_upper = tot_outside_chan_upper,  # floodplain area, wider boundary
         wua_lower, wua_upper,                  # weighted usable area variants
         reach, reach_num)

reaches <- c("San Acacia", "Isleta", "Angostura")

# Build flow-habitat lookup table at 1-cfs resolution (1 to 7,000 cfs).
# NOTE: this lookup extends only to 7,000 cfs (matching the highest modeled
# discharge), unlike the 10,000 cfs range used in 3_scenarios.R. A separate
# combined lookup (q_hab_lookup_combined) extending to 10,000 cfs is built
# later in Section 5D.
q_hab_lookup <- data.frame()

for (i in 1:length(reaches)) {
  q_inund_i   <- filter(q_inund_curve, reach == reaches[i])
  
  # Interpolate all four habitat metrics to 1-cfs resolution using q2hab_2d()
  hab_lookup_i <- tibble(cfs   = 1:7e3,
                         reach = rep(reaches[i], 7e3)) %>%
    mutate(inund_lower = q2hab_2d(cfs, q_inund_i$inund_lower),
           inund_upper = q2hab_2d(cfs, q_inund_i$inund_upper),
           wua_lower   = q2hab_2d(cfs, q_inund_i$wua_lower),
           wua_upper   = q2hab_2d(cfs, q_inund_i$wua_upper))
  
  q_hab_lookup <- bind_rows(q_hab_lookup, hab_lookup_i)
}

# Qdur = 20 days: the duration requirement elicited from expert 3 (EE question
# type "4"). Used as the rolling window width in rollapplyr() below.
Qdur <- 20

inund_cov_list <- list()

for (i in 1:4) {
  # maxmins [17 years x 3 reaches]: peak Qdur-day rolling min of habitat
  # within the April-June window for each year and reach
  maxmins <- array(NA, dim = c(17, 3))
  
  for (j in 1:3) {
    # Reach-gage assignment: San Acacia (j==1) uses SA gage; others use ABQ gage
    if (j == 1) {
      q <- filter(sanaQ, year %in% 2002:2018)
    } else {
      q <- filter(angQ, year %in% 2002:2018)
    }
    
    # Select the habitat metric for this iteration (columns 3-6 = inund_lower,
    # inund_upper, wua_lower, wua_upper; hence names(q_hab_lookup)[i+2])
    hab_lookup_ij <- filter(q_hab_lookup, reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup)[i + 2])
    
    # Join daily habitat to flow data, compute Qdur-day rolling minimum,
    # then find the maximum of that rolling min within April-June (months 4-6).
    # This "max of rolling mins" represents the best sustained inundation event.
    # NOTE: here rollapplyr uses min (not mean as in 3_scenarios.R); both
    # represent the "bottleneck" habitat concept but this version is more
    # conservative (minimum over window vs. mean over window).
    maxmins[, j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
  }
  inund_cov_list[[i]] <- maxmins
}

# ── 5D. Combined inundation covariate (inund_cov_combined_list) ──────────────
# Constructs a hybrid flow-inundation curve merging expert-elicited low-flow
# habitat values with 2D hydraulic model high-flow values, then computes the
# same rolling-min covariate as Section 5C.
#
# Construction logic:
#   1. Extract expert-elicited habitat values at low-flow breakpoints (5-200 cfs)
#      from ee_list.RData (riversegment-coded question types "3e","3c","3a")
#   2. Remove the 700 cfs row from the 2D hydraulic curve (it falls between
#      the expert elicitation range and the hydraulic model range)
#   3. Join expert values to hydraulic values to create a combined breakpoint table
#   4. Where hydraulic data exist, use them; where only expert data exist (< 700
#      cfs), scale using the ratio of hydraulic to expert values at 2,000 cfs
#   5. Interpolate to 1-cfs resolution using q2hab_combined()

# Extract expert-elicited habitat at low-flow breakpoints for all three reaches
elicited_vals <- data.frame()

for (r in 1:3) {
  # t3 here uses the REVERSED sign convention: SD added for lower breakpoints,
  # subtracted for upper (same as preds2 construction above)
  temp <- subset(ee[[3]], ee[[3]][, 1] == riversegment[r] & is.na(ee[[3]][, 5]) == FALSE)
  t3   <- numeric()
  for (t in 1:6) {
    se    <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t3[t] <- temp[t, 6] + se   # lower flows: central + SD
  }
  for (t in 8:15) {
    se    <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t3[t] <- temp[t, 6] - se   # upper flows: central - SD
  }
  t3[7] <- temp[7, 6]  # peak: central value directly
  
  temp$best  <- t3     # store best estimate
  temp$reach <- reaches[r]
  elicited_vals <- bind_rows(elicited_vals, temp)
}

# Reformat elicited values: select the flow breakpoints (X column) and best
# estimates, dropping flows not needed for merging (250 and 2500 cfs rows that
# appear in the EE data but don't have corresponding hydraulic breakpoints)
ee_vals_sub <- select(elicited_vals, q = X, best, reach) %>%
  filter(!q %in% c(250, 2500)) %>%
  mutate(q = as.numeric(q))

# Merge 2D hydraulic curve with elicited values.
# Remove 700 cfs row from hydraulic data: this discharge falls in the gap
# between expert elicitation range and hydraulic model range and would
# create a discontinuity in the combined curve.
q_inund_curve_mod <- q_inund_curve %>%
  select(-reach_num) %>%
  filter(!q %in% c(700)) %>%
  full_join(ee_vals_sub) %>%
  arrange(reach, q) %>%
  pivot_longer(inund_lower:wua_upper,
               names_to  = "inund_type",
               values_to = "hab_orig")

# Compute scaling ratio: ratio of hydraulic habitat to expert habitat at 2,000
# cfs, where both sources overlap. Used to fill expert-elicited habitat at
# low flows proportionally to the hydraulic model's scale.
inund_ratio <- filter(q_inund_curve_mod, q == 2000) %>%
  select(-c(q)) %>%
  mutate(ratio = hab_orig / best) %>%
  select(reach, inund_type, ratio)

# Build combined curve:
#   - Where hydraulic data exist (hab_orig is not NA): use hydraulic value
#   - Where only expert data exist (hab_orig is NA): scale using ratio * best
# The ratio ensures a smooth transition between the two data sources.
q_inund_curve_combined <- q_inund_curve_mod %>%
  right_join(inund_ratio) %>%
  mutate(
    log_hab_orig = case_when(is.na(hab_orig) ~ NA,
                             TRUE             ~ log(hab_orig)),
    hab          = case_when(is.na(hab_orig) ~ ratio * best,  # use scaled expert value
                             TRUE            ~ hab_orig)       # use hydraulic value
  ) %>%
  arrange(reach, inund_type, q) %>%
  select(q, inund_type, hab, reach) %>%
  pivot_wider(names_from = inund_type, values_from = hab)

# Interpolate combined curve to 1-cfs resolution (1 to 10,000 cfs) using
# q2hab_combined(). Extended to 10,000 cfs (vs. 7,000 in q_hab_lookup) to
# accommodate the wider range of the combined breakpoints.
q_hab_lookup_combined <- data.frame()

for (i in 1:length(reaches)) {
  q_inund_i    <- filter(q_inund_curve_combined, reach == reaches[i])
  
  hab_lookup_i <- tibble(cfs   = 1:10e3,
                         reach = rep(reaches[i], 10e3)) %>%
    mutate(inund_lower = q2hab_combined(cfs, q_inund_i$inund_lower),
           inund_upper = q2hab_combined(cfs, q_inund_i$inund_upper),
           wua_lower   = q2hab_combined(cfs, q_inund_i$wua_lower),
           wua_upper   = q2hab_combined(cfs, q_inund_i$wua_upper))
  
  q_hab_lookup_combined <- bind_rows(q_hab_lookup_combined, hab_lookup_i)
}

# Save for use in 3_scenarios.R (used as the baseline lookup for restoration runs)
# write_csv(q_hab_lookup_combined, "data/q_hab_lookup_combined.csv")

# Compute rolling-min covariate from combined curve (same logic as inund_cov_list)
inund_cov_combined_list <- list()

for (i in 1:4) {
  maxmins <- array(NA, dim = c(17, 3))
  
  for (j in 1:3) {
    if (j == 1) {
      q <- filter(sanaQ, year %in% 2002:2018)
    } else {
      q <- filter(angQ, year %in% 2002:2018)
    }
    
    hab_lookup_ij <- filter(q_hab_lookup_combined, reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup_combined)[i + 2])
    
    maxmins[, j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
  }
  inund_cov_combined_list[[i]] <- maxmins
}

# ── 5E. Hydraulic LCC covariate (inund_lcc_list) ─────────────────────────────
# Computes the full LCC covariate (integrating spawning phenology, flow cues,
# and duration requirements) using the 2D hydraulic model's habitat breakpoints
# via calc_cov_2d(). This replaces the expert-elicited habitat values in t3
# with the hydraulic model's inundation values pulled from q_inund_curve.
#
# Loop structure: outer loop over 4 habitat metrics (i), inner loops over 3
# reaches (r) and 17 years (t). For each combination, pull() extracts the
# habitat values at the hydraulic breakpoints for that metric and reach.
# pull(i+1) offsets by 1 because column 1 of q_inund_curve is the discharge
# column (q), so columns 2-5 are the four habitat metrics.

inund_lcc_list <- list()

for (i in 1:4) {
  preds <- array(NA, dim = c(17, 3))
  
  # Reconstruct expert elicitation parameters (t1, t2, t4) -- same as preds2
  temp <- subset(ee[[3]], ee[[3]][, 1] == "1" & is.na(ee[[3]][, 5]) == FALSE)
  t1   <- numeric()
  for (t in 1:6) {
    se         <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t1[t]      <- temp[t, 6] - se
    se         <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100)
    t1[(t + 7)] <- temp[(t + 6), 6] + se
  }
  t1[7]  <- 1
  t1[14] <- 0
  temp <- subset(ee[[3]], ee[[3]][, 1] == "2")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t2   <- temp[1, 6]
  temp <- subset(ee[[3]], ee[[3]][, 1] == "4")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t4   <- temp[1, 6] - se
  
  for (r in 1:3) {
    # t3 here comes from the hydraulic model, not expert elicitation.
    # pull(i+1) selects the i-th habitat metric (columns 2-5 of q_inund_curve).
    # reach_num is used instead of reach name to match the integer indexing.
    t3 <- filter(q_inund_curve, reach_num == r) %>%
      pull(i + 1)  # columns: 1=q, 2=inund_lower, 3=inund_upper, 4=wua_lower, 5=wua_upper
    
    for (t in 1:17) {
      if (r == 3) {q <- aQ[, t]} else {q <- sQ[, t]}
      preds[t, r] <- calc_cov_2d(q, t3, t2, t4, t1)
    }
  }
  
  inund_lcc_list[[i]] <- preds
}

# ── 5F. Combined hydraulic + expert LCC covariate (inund_lcc_list_combined) ──
# Same as Section 5E but uses the combined breakpoint curve (q_inund_combined_format)
# and calc_cov_combined() instead of calc_cov_2d(). This is the covariate used
# in the deployed model (m_lcc_i_l_combined / "lcc_low_combined") and carried
# forward into 3_scenarios.R as inund_lcc_list_combined.

# Reformat combined curve to match the column structure of q_inund_curve,
# with reach in column 5 and reach_num as an integer (3=Angostura, 2=Isleta,
# 1=San Acacia -- note reversed order reflects data arrangement in the source).
q_inund_combined_format <- q_inund_curve_combined %>%
  relocate(reach, .after = wua_upper) %>%
  mutate(reach_num = rep(c(3, 2, 1), each = 14))

# Save for use in 3_scenarios.R berm scenario (as q_inund_combined.csv)
# write_csv(q_inund_combined_format, "data/q_inund_combined.csv")

inund_lcc_list_combined <- list()

for (i in 1:4) {
  preds <- array(NA, dim = c(17, 3))
  
  # Expert elicitation reconstruction: identical to inund_lcc_list
  temp <- subset(ee[[3]], ee[[3]][, 1] == "1" & is.na(ee[[3]][, 5]) == FALSE)
  t1   <- numeric()
  for (t in 1:6) {
    se         <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t1[t]      <- temp[t, 6] - se
    se         <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100)
    t1[(t + 7)] <- temp[(t + 6), 6] + se
  }
  t1[7]  <- 1
  t1[14] <- 0
  temp <- subset(ee[[3]], ee[[3]][, 1] == "2")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t2   <- temp[1, 6]
  temp <- subset(ee[[3]], ee[[3]][, 1] == "4")
  se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
  t4   <- temp[1, 6] - se
  
  for (r in 1:3) {
    # t3 from combined curve (merged expert + hydraulic breakpoints)
    t3 <- filter(q_inund_combined_format, reach_num == r) %>%
      pull(i + 1)  # columns: 1=q, 2=inund_lower, 3=inund_upper, 4=wua_lower, 5=wua_upper
    
    for (t in 1:17) {
      if (r == 3) {q <- aQ[, t]} else {q <- sQ[, t]}
      preds[t, r] <- calc_cov_combined(q, t3, t2, t4, t1)
    }
  }
  
  inund_lcc_list_combined[[i]] <- preds
}


# =============================================================================
# SECTION 6: FIT STAN MODELS
#
# Six models are fit using int3re.stan, an integrated population model (IPM)
# with three reaches and random effects for year-to-year variation. The models
# share the same structure and differ only in the covariate matrix X:
#
#   rgsm.data$X = ct(covariate)
#
# ct() centers the covariate by subtracting the mean. The Stan data list
# includes all monitoring data objects loaded from input_data.RData.
# Each model is fit with 3 chains and 10,000 iterations (5,000 warmup default).
#
# Models are computationally intensive. On a modern workstation with
# parallel::detectCores() chains, expect several hours per model.
# Fitted model objects are saved immediately as .rds for recovery if needed.
# =============================================================================

# Enable parallel chain execution and pre-compile Stan model
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ct(): center a covariate by subtracting its mean. Applied to X before
# passing to Stan so that the intercept (mu_lbeta) represents recruitment
# at mean historical covariate conditions, not at covariate = 0.
ct <- function(x) {x - mean(x)}

# --- Model 1: Original LCC (from Yackulic et al. 2022) ----------------------
# X = centered preds2 (expert-elicited LCC covariate)
rgsm.data <- list(
  Nyears = Nyears, Nstrata = Nstrata, NAVsamps = NAVsamps,
  Nobs_mon1 = Nobs_mon1, Nobs_mon0 = Nobs_mon0, Nobs_mon01 = Nobs_mon01,
  Nobs_monV = Nobs_monV, ewidths = ewidths, Nexperts = Nexperts,
  emove = emove, refQ = refQ,
  X = ct(preds2),  # <-- original LCC covariate
  mon1 = mon1, mon0 = mon0, mon01 = mon01, monV = monV,
  StrataLen = StrataLen, lNz = lNz, lCz = lCz, R0 = R0, R1 = R1,
  cum_nd = cum_nd, Ntotjul = Ntotjul, cum_phiR = cum_phiR,
  NR0 = NR0, NR1 = NR1, w_um = w_um, mAV = mAV, nmons = nmons,
  w_m = w_m, mon1_effQ = mon1_effQ, mon01_effQ = mon01_effQ,
  mon0_effQ = mon0_effQ, monV_effQ = monV_effQ, prop_nd = prop_nd)

M2_1re <- stan("code/int3re.stan", data = rgsm.data, chains = 3, iter = 10000)

# --- Model 2: May-June mean flow (springflow) --------------------------------
# X = centered simpflow (mean May-June discharge / 1000)
rgsm.data <- list(
  Nyears = Nyears, Nstrata = Nstrata, NAVsamps = NAVsamps,
  Nobs_mon1 = Nobs_mon1, Nobs_mon0 = Nobs_mon0, Nobs_mon01 = Nobs_mon01,
  Nobs_monV = Nobs_monV, ewidths = ewidths, Nexperts = Nexperts,
  emove = emove, refQ = refQ,
  X = ct(simpflow),  # <-- May-June flow covariate
  mon1 = mon1, mon0 = mon0, mon01 = mon01, monV = monV,
  StrataLen = StrataLen, lNz = lNz, lCz = lCz, R0 = R0, R1 = R1,
  cum_nd = cum_nd, Ntotjul = Ntotjul, cum_phiR = cum_phiR,
  NR0 = NR0, NR1 = NR1, w_um = w_um, mAV = mAV, nmons = nmons,
  w_m = w_m, mon1_effQ = mon1_effQ, mon01_effQ = mon01_effQ,
  mon0_effQ = mon0_effQ, monV_effQ = monV_effQ, prop_nd = prop_nd)

M2_2re <- stan("code/int3re.stan", data = rgsm.data, chains = 3, iter = 10000)

# --- Model 3: Raw inundation lower boundary (inund_low) ----------------------
# X = log(inund_cov_list[[1]] + 1) centered
# [[1]] selects inund_lower; log transform before centering compresses the
# wide range of inundation values onto an approximately normal scale.
rgsm.data <- list(
  Nyears = Nyears, Nstrata = Nstrata, NAVsamps = NAVsamps,
  Nobs_mon1 = Nobs_mon1, Nobs_mon0 = Nobs_mon0, Nobs_mon01 = Nobs_mon01,
  Nobs_monV = Nobs_monV, ewidths = ewidths, Nexperts = Nexperts,
  emove = emove, refQ = refQ,
  X = ct(log(inund_cov_list[[1]] + 1)),  # <-- log-transformed raw inundation
  mon1 = mon1, mon0 = mon0, mon01 = mon01, monV = monV,
  StrataLen = StrataLen, lNz = lNz, lCz = lCz, R0 = R0, R1 = R1,
  cum_nd = cum_nd, Ntotjul = Ntotjul, cum_phiR = cum_phiR,
  NR0 = NR0, NR1 = NR1, w_um = w_um, mAV = mAV, nmons = nmons,
  w_m = w_m, mon1_effQ = mon1_effQ, mon01_effQ = mon01_effQ,
  mon0_effQ = mon0_effQ, monV_effQ = monV_effQ, prop_nd = prop_nd)

m_inund_i_l <- stan("code/int3re.stan", data = rgsm.data, chains = 3, iter = 10000)

# --- Model 4: Combined inundation lower boundary (inund_low_combined) --------
# X = log(inund_cov_combined_list[[1]] + 1) centered
rgsm.data <- list(
  Nyears = Nyears, Nstrata = Nstrata, NAVsamps = NAVsamps,
  Nobs_mon1 = Nobs_mon1, Nobs_mon0 = Nobs_mon0, Nobs_mon01 = Nobs_mon01,
  Nobs_monV = Nobs_monV, ewidths = ewidths, Nexperts = Nexperts,
  emove = emove, refQ = refQ,
  X = ct(log(inund_cov_combined_list[[1]] + 1)),  # <-- combined inundation
  mon1 = mon1, mon0 = mon0, mon01 = mon01, monV = monV,
  StrataLen = StrataLen, lNz = lNz, lCz = lCz, R0 = R0, R1 = R1,
  cum_nd = cum_nd, Ntotjul = Ntotjul, cum_phiR = cum_phiR,
  NR0 = NR0, NR1 = NR1, w_um = w_um, mAV = mAV, nmons = nmons,
  w_m = w_m, mon1_effQ = mon1_effQ, mon01_effQ = mon01_effQ,
  mon0_effQ = mon0_effQ, monV_effQ = monV_effQ, prop_nd = prop_nd)

m_inund_i_l_combined <- stan("code/int3re.stan", data = rgsm.data, chains = 3, iter = 10000)

# --- Model 5: Hydraulic LCC lower boundary (lcc_low) ------------------------
# X = inund_lcc_list[[1]] centered (no additional log transform; the log is
# already applied inside calc_cov_2d() during covariate construction)
rgsm.data <- list(
  Nyears = Nyears, Nstrata = Nstrata, NAVsamps = NAVsamps,
  Nobs_mon1 = Nobs_mon1, Nobs_mon0 = Nobs_mon0, Nobs_mon01 = Nobs_mon01,
  Nobs_monV = Nobs_monV, ewidths = ewidths, Nexperts = Nexperts,
  emove = emove, refQ = refQ,
  X = ct(inund_lcc_list[[1]]),  # <-- hydraulic LCC covariate
  mon1 = mon1, mon0 = mon0, mon01 = mon01, monV = monV,
  StrataLen = StrataLen, lNz = lNz, lCz = lCz, R0 = R0, R1 = R1,
  cum_nd = cum_nd, Ntotjul = Ntotjul, cum_phiR = cum_phiR,
  NR0 = NR0, NR1 = NR1, w_um = w_um, mAV = mAV, nmons = nmons,
  w_m = w_m, mon1_effQ = mon1_effQ, mon01_effQ = mon01_effQ,
  mon0_effQ = mon0_effQ, monV_effQ = monV_effQ, prop_nd = prop_nd)

m_lcc_i_l <- stan("code/int3re.stan", data = rgsm.data, chains = 3, iter = 10000)

# --- Model 6: Combined hydraulic LCC lower (lcc_low_combined) ---------------
# X = inund_lcc_list_combined[[1]] centered
# This is the DEPLOYED MODEL used in the Shiny app and 3_scenarios.R.
rgsm.data <- list(
  Nyears = Nyears, Nstrata = Nstrata, NAVsamps = NAVsamps,
  Nobs_mon1 = Nobs_mon1, Nobs_mon0 = Nobs_mon0, Nobs_mon01 = Nobs_mon01,
  Nobs_monV = Nobs_monV, ewidths = ewidths, Nexperts = Nexperts,
  emove = emove, refQ = refQ,
  X = ct(inund_lcc_list_combined[[1]]),  # <-- combined hydraulic LCC; DEPLOYED MODEL
  mon1 = mon1, mon0 = mon0, mon01 = mon01, monV = monV,
  StrataLen = StrataLen, lNz = lNz, lCz = lCz, R0 = R0, R1 = R1,
  cum_nd = cum_nd, Ntotjul = Ntotjul, cum_phiR = cum_phiR,
  NR0 = NR0, NR1 = NR1, w_um = w_um, mAV = mAV, nmons = nmons,
  w_m = w_m, mon1_effQ = mon1_effQ, mon01_effQ = mon01_effQ,
  mon0_effQ = mon0_effQ, monV_effQ = monV_effQ, prop_nd = prop_nd)

m_lcc_i_l_combined <- stan("code/int3re.stan", data = rgsm.data, chains = 3, iter = 10000)


# =============================================================================
# SECTION 7: SAVE MODEL OUTPUTS
# =============================================================================

# Load previously fitted models if re-running analysis without refitting.
# These .rds files are the full Stan model objects saved by output_fun() below.
M2_1re               <- readRDS("output/full/orig_lcc.rds")
M2_2re               <- readRDS("output/full/springflow.rds")
m_inund_i_l          <- readRDS("output/full/inund_low.rds")
m_inund_i_l_combined <- readRDS("output/full/inund_low_combined.rds")
m_lcc_i_l            <- readRDS("output/full/lcc_low.rds")
m_lcc_i_l_combined   <- readRDS("output/full/lcc_low_combined.rds")

# Model metadata dataframe: maps short model names to display names and
# pre-allocates R2 columns (filled in Section 8)
model_df <- data.frame(
  model  = c("orig_lcc", "springflow", "inund_low", "inund_low_combined",
             "lcc_low", "lcc_low_combined"),
  name   = c("LCC IPM", "May-June IPM", "Inundation IPM",
             "Composite Inundation IPM", "Hydraulic LCC IPM",
             "Composite Hydraulic LCC IPM"),
  r2_sa  = NA, r2_isl = NA, r2_ang = NA)

# -----------------------------------------------------------------------------
# output_fun()
# Purpose: Save a fitted Stan model and extract key outputs for downstream use.
#
# Saves three outputs per model:
#   (1) Full .rds object: the complete Stan model for diagnostics and reloading
#   (2) *_mcmc_sub.csv: subsampled posterior for the key parameters only
#       (a, mu_lbeta, sd_lbeta, B_lbeta, effS). This is the file loaded in
#       3_scenarios.R as sim_mod for recruitment forecasting.
#   (3) *_summ.csv: summary statistics (mean, SD, quantiles, Rhat, n_eff)
#       for all parameters, used for convergence checking.
#
# Args:
#   mod      -- fitted Stan model object
#   mod_name -- short model name string (from model_df$model)
# -----------------------------------------------------------------------------
output_fun <- function(mod, mod_name) {
  print(mod_name)
  
  # Save full model object
  saveRDS(mod, paste0("output/full/", mod_name, ".rds"))
  
  # Extract subsampled MCMC posterior for key parameters only.
  # pars = c("a", "mu_lbeta", "sd_lbeta", "B_lbeta", "effS") matches the
  # columns expected by forecast_recruit() in 3_scenarios.R.
  # This subsample contains all post-warmup iterations (typically 15,000).
  mcmc_sub <- as.data.frame(mod, pars = c("a", "mu_lbeta", "sd_lbeta",
                                          "B_lbeta", "effS"))
  
  # Summary statistics for all parameters (convergence diagnostics)
  mod_summ <- as.data.frame(summary(mod)$summary)
  
  # Save outputs
  write_csv(mcmc_sub, paste0("output/", mod_name, "_mcmc_sub.csv"))
  write_csv(mod_summ, paste0("output/full/", mod_name, "_summ.csv"))
}

# Ordered list of model objects matching model_df row order
model_list <- list(M2_1re, M2_2re, m_inund_i_l, m_inund_i_l_combined,
                   m_lcc_i_l, m_lcc_i_l_combined)

# Run output_fun() for all models
for (i in 1:length(model_list)) {
  output_fun(model_list[[i]], model_df$model[i])
}

# Save all covariate lists used during model fitting. These are loaded in
# 3_scenarios.R to center predictions using mean(cov_list):
#   inund_cov_list[[1]]          --> used with m_inund_i_l
#   inund_cov_combined_list[[1]] --> used with m_inund_i_l_combined
#   inund_lcc_list[[1]]          --> used with m_lcc_i_l
#   inund_lcc_list_combined[[1]] --> used with m_lcc_i_l_combined (DEPLOYED)
saveRDS(inund_cov_list,          "output/inund_cov_list.RData")
saveRDS(inund_cov_combined_list, "output/inund_cov_combined_list.RData")
saveRDS(inund_lcc_list,          "output/inund_lcc_list.RData")
saveRDS(inund_lcc_list_combined, "output/inund_lcc_list_combined.RData")


# =============================================================================
# SECTION 8: MODEL DIAGNOSTICS
# =============================================================================

# ── 8A. R² calculation ───────────────────────────────────────────────────────
# calc_R2() computes a Bayesian R² for each reach by comparing the variance of
# the model's predictions to the total variance (predictions + residuals).
#
# For reach k:
#   R² = 1 - Var(residuals) / Var(predictions)
#
# where Var is computed across years for each MCMC iteration, then averaged
# across iterations to give a posterior mean R².
#
# eps  = lbeta_eps[iteration, reach, year]: year-level random effect residuals
#        (deviations from the covariate-based prediction on the log-lbeta scale)
# Rf   = Rf[iteration, year, reach]: predicted larval carrying capacity
#        (lbeta on the natural scale); log(Rf) gives the log-scale prediction
#
# A higher R² indicates more of the interannual variation in carrying capacity
# is explained by the covariate, rather than the year random effect.
calc_R2 <- function(mod) {
  eps1   <- apply(rstan::extract(mod, "lbeta_eps")[[1]][, 1, ], 1, var)
  lpred1 <- apply(log(rstan::extract(mod, "Rf")[[1]][, , 1]), 1, var)
  eps2   <- apply(rstan::extract(mod, "lbeta_eps")[[1]][, 2, ], 1, var)
  lpred2 <- apply(log(rstan::extract(mod, "Rf")[[1]][, , 2]), 1, var)
  eps3   <- apply(rstan::extract(mod, "lbeta_eps")[[1]][, 3, ], 1, var)
  lpred3 <- apply(log(rstan::extract(mod, "Rf")[[1]][, , 3]), 1, var)
  
  r2_1 <- 1 - eps1 / lpred1  # San Acacia
  r2_2 <- 1 - eps2 / lpred2  # Isleta
  r2_3 <- 1 - eps3 / lpred3  # Angostura
  
  return(c(mean(r2_1), mean(r2_2), mean(r2_3)))
}

# Compute R² for all models and store in model_df
for (i in 1:length(model_list)) {
  model_df[i, 3:5] <- round(calc_R2(model_list[[i]]), 2)
}

# Format R² table for display and figure
r2_table <- model_df %>%
  mutate(mean_r2 = (r2_sa + r2_isl + r2_ang) / 3) %>%
  select("Model"       = name,
         "San Acacia"  = r2_sa,
         "Isleta"      = r2_isl,
         "Angostura"   = r2_ang,
         "Mean"        = mean_r2)

knitr::kable(r2_table)

# R² heatmap figure: each cell shows the R² value colored by magnitude.
# Horizontal line at y = 4.5 separates the three simpler models (LCC, May-June,
# raw inundation) from the three hydraulic/combined models.
r2s_long <- r2_table %>%
  mutate(Model = factor(Model,
                        levels = rev(c("LCC IPM", "May-June IPM", "Inundation IPM",
                                       "Composite Inundation IPM", "Hydraulic LCC IPM",
                                       "Composite Hydraulic LCC IPM")))) %>%
  pivot_longer(cols = 2:5, names_to = "fit_type", values_to = "r2") %>%
  mutate(fit_type = factor(fit_type,
                           levels = c("San Acacia", "Isleta", "Angostura", "Mean")),
         r2_label = sprintf("%.2f", r2))

sprintf("%02d", label_val)  # format helper for cell labels

ggplot(r2s_long, aes(fit_type, Model)) +
  geom_tile(aes(fill = r2)) +
  geom_hline(aes(yintercept = 4.5)) +  # divides simple from hydraulic models
  geom_text(aes(label = r2_label)) +
  scale_fill_viridis_c() +
  labs(x = "Reach", fill = expression(R ^ 2)) +
  theme_ipsum() +
  theme(plot.margin  = unit(c(1, 0, 1, 1), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12))

#ggsave("plots/r2_mods.jpeg", height = 2, width = 6, units = "in")

# ── 8B. MCMC traceplots ───────────────────────────────────────────────────────
# Visually inspect chain convergence for key parameters. Replace model and
# parameter names as needed. Well-mixed chains (no trends or divergences)
# indicate adequate convergence. B_lbeta is the most important parameter
# for the recruitment tool as it governs the effect of inundation on carrying
# capacity; poor convergence here would invalidate scenario projections.
traceplot(m_inund_i_l_combined, pars = "B_lbeta", inc_warmup = TRUE)

# ── 8C. Posterior predictive checks (from Yackulic et al. 2022) ──────────────
# mod_check_fun() generates 8-panel posterior predictive check figure (Fig S1
# in the manuscript). Checks whether the model's predicted catch distribution
# is consistent with observed catch across monitoring and rescue datasets.
#
# Panels A-B: Observed vs. predicted catch on log scale (R² shown)
# Panels C-D: Distribution of proportion of catches equal to zero (% zeros test)
#             Red line = observed proportion; histogram = simulated distribution
#             A well-fitting model's observed % zeros should fall within the
#             simulated distribution, not in the tails.
# Panels E-F: QQ plot of observed quantiles vs. uniform expected quantiles
#             Under a well-calibrated model, points should fall on the 1:1 line.
# Panels G-H: Observed quantile vs. discharge (G) and Julian date (H)
#             Systematic patterns (e.g., low quantiles at high flow) would
#             indicate a flow-related model misspecification.
#
# The qf() function computes the probability integral transform: the fraction
# of simulated values below the observed value. Under a well-calibrated model,
# these should be uniformly distributed between 0 and 1.
mod_check_fun <- function(model) {
  
  # Extract predicted log-catch and log-rescue densities from posterior
  lp0  <- rstan::extract(model, "lpmon0")[[1]]   # monitoring type 0 (undetected)
  lp1  <- rstan::extract(model, "lpmon1")[[1]]   # monitoring type 1 (detected)
  lp01 <- rstan::extract(model, "lpmon01")[[1]]  # monitoring type 01 (mixed)
  lpV  <- rstan::extract(model, "lpmonV")[[1]]   # vegetation monitoring
  lR0  <- rstan::extract(model, "lpR0")[[1]]     # rescue event type 0
  lR1  <- rstan::extract(model, "lpR1")[[1]]     # rescue event type 1
  sz   <- rstan::extract(model, "sz")[[1]]       # negative binomial size (monitoring)
  rsz  <- rstan::extract(model, "rsz")[[1]]      # negative binomial size (rescue)
  
  pc0  <- lp0; pc1 <- lp1; pc01 <- lp01; pcV <- lpV
  pr0  <- lR0; pr1 <- lR1
  
  # Quantile function (probability integral transform):
  # returns the fraction of simulated values <= observed value.
  # Handles ties (observed == simulated) by sampling uniformly from tied ranks.
  # Handles out-of-range values (below all simulated) by returning 1/n.
  q_pc0  <- numeric(); q_pc1 <- numeric(); q_pc01 <- numeric()
  q_pcV  <- numeric(); q_pr0 <- numeric(); q_pr1  <- numeric()
  
  qf <- function(obs, pred) {
    temp <- which(obs == sort(pred))
    if (length(temp) > 0) {
      sample(temp, 1) / length(pred)
    } else {
      temp2 <- findInterval(obs, sort(pred))
      if (temp2 == 0) {1 / length(pred)} else {temp2 / length(pred)}
    }
  }
  
  # Simulate predicted catches from the negative binomial model for each
  # observation, then compute the probability integral transform
  for (i in 1:(dim(lp0)[2]))  {pc0[, i]  <- rnbinom(15000, mu = exp(lp0[, i]),  size = sz);  q_pc0[i]  <- qf(mon0[i, 5],  pc0[, i])}
  for (i in 1:(dim(lp1)[2]))  {pc1[, i]  <- rnbinom(15000, mu = exp(lp1[, i]),  size = sz);  q_pc1[i]  <- qf(mon1[i, 5],  pc1[, i])}
  for (i in 1:(dim(lp01)[2])) {pc01[, i] <- rnbinom(15000, mu = exp(lp01[, i]), size = sz);  q_pc01[i] <- qf(mon01[i, 5], pc01[, i])}
  for (i in 1:(dim(lpV)[2]))  {pcV[, i]  <- rnbinom(15000, mu = exp(lpV[, i]),  size = sz);  q_pcV[i]  <- qf(monV[i, 5],  pcV[, i])}
  for (i in 1:(dim(lR0)[2]))  {pr0[, i]  <- rnbinom(15000, mu = exp(lR0[, i]),  size = rsz); q_pr0[i]  <- qf(R0[i, 4],   pr0[, i])}
  for (i in 1:(dim(lR1)[2]))  {pr1[, i]  <- rnbinom(15000, mu = exp(lR1[, i]),  size = rsz); q_pr1[i]  <- qf(R1[i, 4],   pr1[, i])}
  
  # Combine monitoring data types for plotting
  q_mon    <- c(q_pc0, q_pc1, q_pc01, q_pcV)
  pc_mon   <- cbind(pc0, pc1, pc01, pcV)
  obsc_mon <- c(mon0[, 5], mon1[, 5], mon01[, 5], monV[, 5])
  Q_mon    <- c(mon0_effQ[, 2], mon1_effQ[, 2], mon01_effQ[, 2], monV_effQ[, 2])
  q_res    <- c(q_pr0, q_pr1)
  pc_res   <- cbind(pr0, pr1)
  obsc_res <- c(R0[, 4], R1[, 4])
  jul_res  <- c(R0[, 2], R1[, 2])
  
  # 8-panel diagnostic figure
  par(mfrow = c(4, 2), mar = c(4, 4, 1, 1))
  
  # Panels A-B: Observed vs. predicted mean catch (log scale)
  xymax <- max(c(obsc_mon, apply(pc_mon, 2, mean)))
  plot(1 + obsc_mon, 1 + apply(pc_mon, 2, mean),
       ylim = c(1, 1 + xymax), xlim = c(1, 1 + xymax), axes = FALSE,
       xlab = "", ylab = "", pch = 19, col = rgb(0, 0, 0, .1),
       main = "Monitoring data", log = "xy")
  axis(1, at = c(1, 2, 11, 101, 1001, 10001), labels = c(0, 1, 10, 100, 1000, 10000), pos = 1)
  axis(2, at = c(1, 2, 11, 101, 1001, 10001), labels = c(0, 1, 10, 100, 1000, 10000), pos = 1, las = TRUE)
  mtext("Observed catch", 1, 2, cex = 0.7)
  mtext("Predicted catch", 2, 2, cex = 0.7)
  text(1000, 3, "R2 = 0.76")
  text(.3, 3000, "A)", xpd = TRUE)
  
  xymax <- max(c(obsc_res, apply(pc_res, 2, mean)))
  plot(1 + obsc_res, 1 + apply(pc_res, 2, mean),
       ylim = c(1, 1 + xymax), xlim = c(1, 1 + xymax), axes = FALSE,
       xlab = "", ylab = "", pch = 19, col = rgb(0, 0, 0, .1),
       main = "Rescue data", log = "xy")
  axis(1, at = c(1, 2, 11, 101, 1001, 10001), labels = c(0, 1, 10, 100, 1000, 10000), pos = 1)
  axis(2, at = c(1, 2, 11, 101, 1001, 10001), labels = c(0, 1, 10, 100, 1000, 10000), pos = 1, las = TRUE)
  mtext("Observed catch", 1, 2, cex = 0.7)
  mtext("Predicted catch", 2, 2, cex = 0.7)
  text(8000, 3, "R2 = 0.77")
  text(.15, 30000, "B)", xpd = TRUE)
  
  # Panels C-D: % zeros check -- histogram of simulated % zeros, red line = observed
  iszero <- function(x) {ifelse(x == 0, 1, 0)}
  hist(apply(iszero(pc_mon), 1, mean), axes = FALSE, xlab = "", ylab = "", main = "")
  abline(v = mean(iszero(obsc_mon)), col = "red", lty = 2, lwd = 2)
  axis(1, pos = 0)
  axis(2, pos = min(apply(iszero(pc_mon), 1, mean)),
       at = c(0, 750, 1500, 2250), labels = seq(0, 0.15, .05), las = TRUE)
  mtext("Proportion of catch equal to zero", 1, 2, cex = 0.7)
  mtext("% of simulations", 2, 2, cex = 0.7)
  text(mean(iszero(obsc_mon)) + .01, 1500, "Observed", col = "red")
  text(0.52, 2600, "C)", xpd = TRUE)
  
  hist(apply(iszero(pc_res), 1, mean), axes = FALSE, xlab = "", ylab = "", main = "")
  abline(v = mean(iszero(obsc_res)), col = "red", lty = 2, lwd = 2)
  axis(1, pos = 0)
  axis(2, pos = min(apply(iszero(pc_res), 1, mean)),
       at = c(0, 1500, 3000, 4500), labels = seq(0, 0.3, .1), las = TRUE)
  mtext("Proportion of catch equal to zero", 1, 2, cex = 0.7)
  mtext("% of simulations", 2, 2, cex = 0.7)
  text(mean(iszero(obsc_res)) - 0.03, 4000, "Observed", col = "red")
  text(0.14, 4500, "D)", xpd = TRUE)
  
  # Panels E-F: QQ plots of observed quantile vs. expected uniform quantile
  plot(sort(q_mon), c(1:length(q_mon)) / length(q_mon),
       pch = 19, col = rgb(0, 0, 0, .01), main = "", xlab = "", ylab = "", axes = FALSE)
  curve(1 * x, add = TRUE, lty = 2, lwd = 2)  # reference 1:1 line
  axis(1, pos = 0); axis(2, pos = 0, las = TRUE)
  mtext("Observed quantile", 1, 2, cex = 0.7)
  mtext("Expected quantile", 2, 2, cex = 0.7)
  text(-0.16, 1.1, "E)", xpd = TRUE)
  
  plot(sort(q_res), c(1:length(q_res)) / length(q_res),
       pch = 19, col = rgb(0, 0, 0, .1), main = "", xlab = "", ylab = "", axes = FALSE)
  curve(1 * x, add = TRUE, lty = 2, lwd = 2)
  axis(1, pos = 0); axis(2, pos = 0, las = TRUE)
  mtext("Observed quantile", 1, 2, cex = 0.7)
  mtext("Expected quantile", 2, 2, cex = 0.7)
  text(-0.16, 1.1, "F)", xpd = TRUE)
  
  # Panels G-H: Observed quantile vs. discharge and Julian date
  # Systematic patterns here would indicate residual flow or timing effects
  # not captured by the covariate
  plot(Q_mon, q_mon, pch = 19, col = rgb(0, 0, 0, .1), main = "", xlab = "", ylab = "", axes = FALSE)
  axis(1, at = c(0, 0.5, 1), labels = c(0, 500, 1000), pos = 0)  # Q stored as cfs/1000
  axis(2, pos = 0, las = TRUE)
  mtext("Discharge (cfs)", 1, 2, cex = 0.7)
  mtext("Observed quantile", 2, 2, cex = 0.7)
  text(-0.16, 1.1, "G)", xpd = TRUE)
  
  plot(jul_res, q_res, pch = 19, col = rgb(0, 0, 0, .1), main = "", xlab = "", ylab = "", axes = FALSE)
  axis(1, pos = 0)
  axis(2, pos = 0, las = TRUE)
  mtext("Days after April 1st", 1, 2, cex = 0.7)
  mtext("Observed quantile", 2, 2, cex = 0.7)
  text(-27, 1.1, "H)", xpd = TRUE)
}

mod_check_fun(M2_2re)  # run diagnostic figure for the springflow model

# ── 8D. Prior vs. posterior checks ───────────────────────────────────────────
# checkpp_rep(): compare posterior distributions (red) to prior distributions
# (blue dashed) for three key recruitment parameters:
#   a (density-independent survival): prior N(675, 135)
#   beta_stk (stocking effect multiplier): prior N(1, 0.1)
#   beta_2 (age-2+ fish vs. age-1 spawner effectiveness): prior N(2, 0.1)
# Large posterior-prior differences indicate the data are informative about
# these parameters; near-identical distributions indicate non-identification.
checkpp_rep <- function(mod) {
  par(mfrow = c(1, 3))
  
  plot(density(rstan::extract(mod, "a")[[1]]),
       xlim = c(300, 1500), axes = FALSE, main = "", xlab = "a", col = "red")
  curve(dnorm(x, 675, 135), col = "blue", add = TRUE, lty = 2)
  axis(1, pos = 0); axis(2, pos = 300, las = TRUE)
  
  plot(density(rstan::extract(mod, "beta_stk")[[1]]),
       xlim = c(0.5, 3), axes = FALSE, main = "", xlab = "beta_stk", col = "red")
  curve(dnorm(x, 1, 0.1), col = "blue", add = TRUE, lty = 2)
  axis(1, pos = 0); axis(2, pos = 0.5, las = TRUE)
  
  plot(density(rstan::extract(mod, "beta_2")[[1]]),
       xlim = c(0.5, 3), axes = FALSE, main = "", xlab = "beta_2", col = "red")
  curve(dnorm(x, 2, 0.1), col = "blue", add = TRUE, lty = 2)
  axis(1, pos = 0); axis(2, pos = 0.5, las = TRUE)
}

checkpp_rep(M2_2re)     # springflow model
checkpp_rep(m_inund_i_l) # raw inundation model

# checkpp_sd(): compare posteriors to priors for mortality and detection
# parameters. All have Uniform(lower, upper) priors, so the prior is flat
# across the plotted range. Parameters:
#   mu_lM0/mu_lM1/mu_lMw: log-scale mortality rates for age-0, age-1+, winter
#                          (3 reaches each); prior Uniform(-8, -3)
#   irphi, p0, p1, rp0, rp1: detection probabilities; prior Uniform(0, 1)
#   sz, rsz: negative binomial overdispersion; prior Uniform(0, 1)
checkpp_sd <- function(mod) {
  par(mfrow = c(4, 4))
  for (k in 1:3) {
    plot(density(rstan::extract(mod, "mu_lM0")[[1]][, k]),
         xlim = c(-8, -3), axes = FALSE, main = "", xlab = paste("mu_lM0", k), col = "red")
    curve(dunif(x, -8, -3), col = "blue", add = TRUE, lty = 2)
    axis(1, pos = 0); axis(2, pos = -8, las = TRUE)
  }
  for (k in 1:3) {
    plot(density(rstan::extract(mod, "mu_lM1")[[1]][, k]),
         xlim = c(-8, -3), axes = FALSE, main = "", xlab = paste("mu_lM1", k), col = "red")
    curve(dunif(x, -8, -3), col = "blue", add = TRUE, lty = 2)
    axis(1, pos = 0); axis(2, pos = -8, las = TRUE)
  }
  for (k in 1:3) {
    plot(density(rstan::extract(mod, "mu_lMw")[[1]][, k]),
         xlim = c(-8, -3), axes = FALSE, main = "", xlab = paste("mu_lMw", k), col = "red")
    curve(dunif(x, -8, -3), col = "blue", add = TRUE, lty = 2)
    axis(1, pos = 0); axis(2, pos = -8, las = TRUE)
  }
  for (param in c("irphi", "p0", "p1", "rp0", "rp1")) {
    plot(density(rstan::extract(mod, param)[[1]]),
         xlim = c(0, 1), axes = FALSE, main = "", xlab = param, col = "red")
    curve(dunif(x, 0, 1), col = "blue", add = TRUE, lty = 2)
    axis(1, pos = 0); axis(2, pos = 0, las = TRUE)
  }
  for (param in c("sz", "rsz")) {
    plot(density(rstan::extract(mod, param)[[1]]),
         main = "", xlab = param, col = "red", axes = FALSE, xlim = c(0.25, 0.6))
    curve(dunif(x, 0, 1), col = "blue", add = TRUE, lty = 2)
    axis(1, pos = 0); axis(2, pos = 0.25, las = TRUE)
  }
}

checkpp_sd(M2_1re)  # run for original LCC model


# =============================================================================
# SECTION 9: OUT-OF-SAMPLE VALIDATION (2019-2024)
#
# Tests whether models fit on 2002-2018 data accurately predict RGSM catch
# during the validation period (2019-2024). Six model variants plus the
# Biological Opinion linear model are compared against observed annual catch.
#
# Approach:
#   1. Load 2019-2024 RGSM monitoring catch data (oos_data_new.RData)
#   2. Download 2021-2024 USGS flow data and append to historical records
#   3. Compute covariates for 2019-2024 using the same methods as fitting
#   4. Run forecast_oos_re() for each model to simulate 2019-2024 catch
#   5. Compare predicted vs. observed total annual catch (Figure 3)
# =============================================================================

load(file = "output/oos_data_new.RData")

# Download 2021-2024 USGS flow data for both gages using the dataRetrieval package.
# USGS parameter code 00060 = instantaneous discharge; readNWISdv fetches daily values.
# Reformatted to match the column structure of the historical angQ/sanaQ objects.
start.date <- "2021-01-01"
end.date   <- "2024-12-31"
siteAng    <- "08330000"  # Albuquerque gage (for Angostura reach)
siteSan    <- "08354900"  # San Acacia gage
pCode      <- "00060"     # USGS parameter code for discharge

new_ang_data <- readNWISdv(siteNumbers = siteAng,
                           parameterCd = pCode,
                           startDate   = start.date,
                           endDate     = end.date) %>%
  mutate(year  = as.numeric(format(Date, "%Y")),
         month = as.numeric(format(Date, "%m")),
         day   = as.numeric(format(Date, "%d")),
         Date  = paste(month, day, year, sep = "/")) %>%
  rename(cfs = X_00060_00003) %>%
  select(names(angQ))  # enforce column order matching historical data

new_SanA_data <- readNWISdv(siteNumbers = siteSan,
                            parameterCd = pCode,
                            startDate   = start.date,
                            endDate     = end.date) %>%
  mutate(year  = as.numeric(format(Date, "%Y")),
         month = as.numeric(format(Date, "%m")),
         day   = as.numeric(format(Date, "%d")),
         Date  = paste(month, day, year, sep = "/")) %>%
  rename(cfs = X_00060_00003) %>%
  select(names(sanaQ))

# Combine historical (2002-2018 from input_data.RData) with new (2021-2024)
# to create full-period gage records for out-of-sample covariate calculation
angQ_all  <- bind_rows(angQ,  new_ang_data)
sanaQ_all <- bind_rows(sanaQ, new_SanA_data)

# -----------------------------------------------------------------------------
# forecast_oos_re()
# Purpose: Project RGSM population dynamics and catch for 2019-2024, using
#   the fitted Stan model parameters and out-of-sample covariate values.
#
# This function implements the full integrated population model (IPM) forward
# simulation for 6 years, carrying population state (SpN) from year to year.
# It mirrors the Stan model's structure but is implemented in R for flexibility.
#
# Key parameters extracted from the Stan model:
#   Sp_N0:    initial spawner abundance in 2019 by age class and reach
#             (taken from the year-18 posterior, i.e., the last fitted year)
#   a:        Beverton-Holt density-independent survival
#   mu_lbeta, sd_lbeta, B_lbeta: larval carrying capacity parameters
#   mu_lM0/M1/Mw: age-class-specific mortality rate parameters
#   move, rp0, rp1: movement and detection probability parameters
#   A0/AtQ_perpool, sl_width, bankfull: habitat-catch relationship parameters
#   alpha0/1_max/int: detection efficiency parameters
#   p0, p1: sampling probabilities
#   sz: negative binomial overdispersion
#   irphi: stocking survival factor
#
# Population state transitions (year 1 to 6):
#   SpN[age-0]: newly recruited fish, predicted by Beverton-Holt from effS and tRf
#   SpN[age-1+]: surviving fish from previous year FN (post-July)
#   SpN[stocked]: stocked fish scaled by irphi and winter mortality (Mw)
#   effS: effective spawner abundance = age-0 + age-1 * beta_2 + stocked * beta_stk
#   tRf:  larval carrying capacity = exp(mu_lbeta + sd_lbeta noise + B_lbeta * Xout)
#   tN0:  age-0 recruits = Beverton-Holt(a, effS, tRf)
#
# Catch prediction:
#   pC[i]: predicted catch probability for each out-of-sample monitoring event i,
#          accounting for pool/run habitat composition, detection efficiency
#          by habitat type (talpha0/1), and sampling effort.
#   C[i]:  simulated catch drawn from NegBinom(mu = pC[i], size = sz).
#
# Args:
#   mod  -- fitted Stan model object
#   Xout -- centered LCC covariate matrix [6 years x 3 reaches] for 2019-2024
#
# Returns: list with:
#   ef_6yr:     total electrofishing effort by year
#   pc_6yr:     predicted catch by year [iter x 6]
#   c_6yr:      simulated catch by year [iter x 6]
#   SpN:        spawner abundance [iter x 6 x 3 age classes x 3 reaches]
#   FN:         July-1 abundance [iter x 6 x 2 size classes x 3 reaches]
#   effS:       effective spawner abundance [iter x 6 x 3 reaches]
#   jul1_age0:  age-0 recruits on July 1 [iter x 6 x 3 reaches]
#   tRf_arr:    larval carrying capacity [iter x 6 x 3 reaches]
#   c_SA/Isl/Ang: simulated catch by reach and year [iter x 6 each]
# -----------------------------------------------------------------------------
forecast_oos_re <- function(mod, Xout) {
  
  # Extract all required parameter vectors from the posterior
  Sp_N0       <- rstan::extract(mod, "Sp_N")[[1]][, 18, , ]  # year 18 = 2018 (last fitted year)
  irphi       <- rstan::extract(mod, "irphi")[[1]]           # stocking survival factor
  beta_2      <- rstan::extract(mod, "beta_2")[[1]]          # age-2+ spawner effectiveness
  beta_stk    <- rstan::extract(mod, "beta_stk")[[1]]        # stocked fish effectiveness
  a           <- rstan::extract(mod, "a")[[1]]               # Beverton-Holt a parameter
  mu_lM0      <- rstan::extract(mod, "mu_lM0")[[1]]          # log mortality rate, age-0 [iter x 3]
  mu_lM1      <- rstan::extract(mod, "mu_lM1")[[1]]          # log mortality rate, age-1+ [iter x 3]
  mu_lMw      <- rstan::extract(mod, "mu_lMw")[[1]]          # log mortality rate, winter [iter x 3]
  sd_lM       <- rstan::extract(mod, "sd_lM")[[1]]           # among-year SD in log mortality
  move        <- rstan::extract(mod, "move")[[1]]            # movement probability
  rp0         <- rstan::extract(mod, "rp0")[[1]]             # rescue detection probability (age-0)
  rp1         <- rstan::extract(mod, "rp1")[[1]]             # rescue detection probability (age-1+)
  mu_lbeta    <- rstan::extract(mod, "mu_lbeta")[[1]]        # log LCC intercept [iter x 3]
  sd_lbeta    <- rstan::extract(mod, "sd_lbeta")[[1]]        # among-year SD in log LCC
  B_lbeta     <- rstan::extract(mod, "B_lbeta")[[1]]         # LCC covariate coefficient
  A0_perpool  <- rstan::extract(mod, "A0_perpool")[[1]]      # pool density-area relationship intercept
  AtQ_perpool <- rstan::extract(mod, "AtQ_perpool")[[1]]     # pool area-discharge coefficient
  sl_width    <- rstan::extract(mod, "sl_width")[[1]]        # stream width at low flow [iter x 3]
  bankfull    <- rstan::extract(mod, "bankfull")[[1]]        # bankfull discharge [iter x 3]
  alpha0_max  <- rstan::extract(mod, "alpha0_max")[[1]]      # max detection in pools (age-0)
  alpha1_max  <- rstan::extract(mod, "alpha1_max")[[1]]      # max detection in runs (age-1+)
  alpha0_int  <- rstan::extract(mod, "alpha0_int")[[1]]      # detection-discharge coefficient (age-0)
  alpha1_int  <- rstan::extract(mod, "alpha1_int")[[1]]      # detection-discharge coefficient (age-1+)
  p0          <- rstan::extract(mod, "p0")[[1]]              # base sampling probability (age-0)
  p1          <- rstan::extract(mod, "p1")[[1]]              # base sampling probability (age-1+)
  sz          <- rstan::extract(mod, "sz")[[1]]              # NegBinom overdispersion
  sd_lbeta    <- rstan::extract(mod, "sd_lbeta")[[1]]
  iter        <- length(p0)                                   # number of MCMC samples
  
  # Pre-allocate mortality arrays [iter x 3 reaches x 6 years]
  Mw <- array(NA, dim = c(iter, 3, 6))
  M0 <- array(NA, dim = c(iter, 3, 6))
  M1 <- array(NA, dim = c(iter, 3, 6))
  
  # Draw mortality rates for each year from the posterior predictive distribution.
  # Mortality rates are drawn fresh for each out-of-sample year (no pooling across
  # years) because these years were not in the fitting data.
  for (k in 1:Nstrata) {
    for (yr in 1:6) {
      Mw[, k, yr] <- exp(rnorm(iter, mu_lMw[, k], sd_lM))
      M0[, k, yr] <- exp(rnorm(iter, mu_lM0[, k], sd_lM))
      M1[, k, yr] <- exp(rnorm(iter, mu_lM1[, k], sd_lM))
    }
  }
  
  # Pre-allocate population state and catch arrays
  FN        <- array(NA, dim = c(iter, 6, 2, 3))    # July-1 abundance [iter x yr x age x reach]
  pC        <- matrix(NA, ncol = Nobs_oosC, nrow = iter)  # predicted catch probability
  C         <- matrix(NA, ncol = Nobs_oosC, nrow = iter)  # simulated catch
  SpN       <- array(NA, dim = c(iter, 6, 3, 3))    # spawner abundance [iter x yr x age x reach]
  effS      <- array(NA, dim = c(iter, 6, 3))        # effective spawners [iter x yr x reach]
  jul1_age0 <- array(NA, dim = c(iter, 6, 3))        # age-0 on July 1
  tRf_arr   <- array(NA, dim = c(iter, 6, 3))        # larval carrying capacity
  
  for (k in 1:Nstrata) {
    
    # Year 1 (2019): initialize from model's year-18 posterior (2018 fish)
    SpN[, 1, 1, k] <- Sp_N0[, k, 1]  # age-0 cohort from 2018
    SpN[, 1, 2, k] <- Sp_N0[, k, 2]  # age-1+ cohort from 2018
    SpN[, 1, 3, k] <- (w_m_oos[1, k, 1] * exp(-60  * Mw[, k, 1]) +
                         w_m_oos[1, k, 2] * exp(-151 * Mw[, k, 1])) * irphi  # stocked fish
    
    # Effective spawners: weighted sum of age classes
    effS[, 1, k] <- SpN[, 1, 1, k] + SpN[, 1, 2, k] * beta_2 + SpN[, 1, 3, k] * beta_stk
    tN1          <- rowSums(SpN[, 1, 1:2, k])  # total age-1+ spawners
    
    # Larval carrying capacity: exp(mu_lbeta + noise + B_lbeta * Xout)
    tRf              <- exp(rnorm(iter, mu_lbeta[, k], sd_lbeta) + B_lbeta * Xout[1, k])
    tRf_arr[, 1, k]  <- tRf
    tN0              <- a * (effS[, 1, k]) / (1 + a * (effS[, 1, k]) / tRf)  # Beverton-Holt
    jul1_age0[, 1, k] <- tN0
    
    # July-1 abundance after survival to sampling date
    # FN[age-0]: age-0 fish survival accounting for movement and downstream drift (cum_nd)
    # FN[age-1+]: age-1+ fish survival to end of year
    FN[, 1, 1, k] <- tN0 * exp(-M0[, k, 1] * 124) *
      ((1 - cum_nd_oos[1, 91, k]) + (move - 1) * (cum_nd_oos[1, 215, k] - cum_nd_oos[1, 91, k]) +
         (1 - move) * rp0 * (cum_phiR_oos[1, 215, k] - cum_phiR_oos[1, 91, k]))
    FN[, 1, 2, k] <- tN1 * exp(-M1[, k, 1] * 215) *
      (1 + (move - 1) * cum_nd_oos[1, 215, k] + (1 - move) * rp1 * cum_phiR_oos[1, 215, k])
    
    # Years 2-6 (2020-2024): carry state forward from previous year
    for (t in 2:6) {
      SpN[, t, 1:3, k] <- cbind(
        FN[, t - 1, 1, k] * exp(-Mw[, k, t] * 150),              # age-0 overwinter survival
        FN[, t - 1, 2, k] * exp(-Mw[, k, t] * 150),              # age-1+ overwinter survival
        (w_m_oos[t, k, 1] * exp(-60 * Mw[, k, t]) +
           w_m_oos[t, k, 2] * exp(-151 * Mw[, k, t])) * irphi)   # stocked fish in year t
      
      effS[, t, k] <- SpN[, t, 1, k] + SpN[, t, 2, k] * beta_2 + SpN[, t, 3, k] * beta_stk
      tN1          <- rowSums(SpN[, t, 1:2, k])
      tRf          <- exp(rnorm(iter, mu_lbeta[, k], sd_lbeta) + B_lbeta * Xout[t, k])
      tRf_arr[, t, k] <- tRf
      tN0          <- a * (effS[, t, k]) / (1 + a * (effS[, t, k]) / tRf)
      jul1_age0[, t, k] <- tN0
      
      FN[, t, 1, k] <- tN0 * exp(-M0[, k, t] * 124) *
        ((1 - cum_nd_oos[t, 91, k]) + (move - 1) * (cum_nd_oos[t, 215, k] - cum_nd_oos[t, 91, k]) +
           (1 - move) * rp0 * (cum_phiR_oos[t, 215, k] - cum_phiR_oos[t, 91, k]))
      FN[, t, 2, k] <- tN1 * exp(-M1[, k, t] * 215) *
        (1 + (move - 1) * cum_nd_oos[t, 215, k] + (1 - move) * rp1 * cum_phiR_oos[t, 215, k])
    }
  }
  
  # Predict catch for each out-of-sample monitoring event.
  # The catch model accounts for:
  #   - habitat composition at each sample (pool vs. run, via A0/AtQ/sl_width/bankfull)
  #   - detection efficiency by habitat type and fish size class (talpha0/1)
  #   - sampling effort (oos_C$effort)
  #   - fish abundance available at sampling time (FN[y, age, reach])
  #   - survival from July 1 to sampling date
  for (i in 1:Nobs_oosC) {
    q       <- oos_C$cQ[i] / 1000  # discharge at sampling event, scaled to match model
    totpool <- exp(A0_perpool + AtQ_perpool * q) * 200 * sl_width[oos_C[i, 5]] * q /
      (1 + sl_width[oos_C[i, 5]] * q / bankfull[oos_C[i, 5]])
    totrun  <- (1 - exp(A0_perpool + AtQ_perpool * q)) * 200 * sl_width[oos_C[i, 5]] * q /
      (1 + sl_width[oos_C[i, 5]] * q / bankfull[oos_C[i, 5]])
    
    # Detection efficiency by habitat type (pool = type 1, run = type 2)
    talpha0 <- c(alpha0_int * q / (1 + alpha0_int * q / alpha0_max), 1)
    talpha1 <- c(alpha1_int * q / (1 + alpha1_int * q / alpha1_max), 1)
    
    y <- oos_C[i, 3] - 2018  # year index: 2019 = 1, 2020 = 2, etc.
    
    # Predicted catch probability: sum of age-0 and age-1+ contributions,
    # each weighted by detection efficiency, effort, and habitat composition
    pC[, i] <- (p0 * oos_C$effort[i] * talpha0[oos_C$type[i]] * .2 *
                  FN[, y, 1, oos_C[i, 5]] *
                  exp(-M0[, oos_C[i, 5], (oos_C[i, 3] - 2018)] * (oos_C[i, 2] - 91)) *
                  ((1 - cum_nd_oos[y, 91, oos_C[i, 5]]) +
                     (move - 1) * (cum_nd_oos[y, oos_C[i, 2], oos_C[i, 5]] - cum_nd_oos[y, 91, oos_C[i, 5]]) +
                     (1 - move) * rp0 * (cum_phiR_oos[y, oos_C[i, 2], oos_C[i, 5]] - cum_phiR_oos[y, 91, oos_C[i, 5]]))) /
      (StrataLen_oos[y, oos_C[i, 2], oos_C[i, 5]] * (talpha0[1] * totpool + totrun)) +
      (p1 * oos_C$effort[i] * talpha1[oos_C$type[i]] * .2 *
         FN[, y, 2, oos_C[i, 5]] *
         exp(-M1[, oos_C[i, 5], (oos_C[i, 3] - 2018)] * oos_C[i, 2]) *
         (1 + (move - 1) * cum_nd_oos[y, oos_C[i, 2], oos_C[i, 5]] +
            (1 - move) * rp1 * cum_phiR_oos[y, oos_C[i, 2], oos_C[i, 5]])) /
      (StrataLen_oos[y, oos_C[i, 2], oos_C[i, 5]] * (talpha1[1] * totpool + totrun))
    
    C[, i] <- rnbinom(iter, mu = pC[, i], size = sz)
  }
  
  # Index rows of oos_C by year and reach for annual catch aggregation
  w19 <- which(oos_C$year == 2019); w20 <- which(oos_C$year == 2020)
  w21 <- which(oos_C$year == 2021); w22 <- which(oos_C$year == 2022)
  w23 <- which(oos_C$year == 2023); w24 <- which(oos_C$year == 2024)
  
  wSA19  <- which(oos_C$year == 2019 & oos_C$Cstrata == 1)
  wIsl19 <- which(oos_C$year == 2019 & oos_C$Cstrata == 2)
  wAng19 <- which(oos_C$year == 2019 & oos_C$Cstrata == 3)
  wSA20  <- which(oos_C$year == 2020 & oos_C$Cstrata == 1)
  wIsl20 <- which(oos_C$year == 2020 & oos_C$Cstrata == 2)
  wAng20 <- which(oos_C$year == 2020 & oos_C$Cstrata == 3)
  wSA21  <- which(oos_C$year == 2021 & oos_C$Cstrata == 1)
  wIsl21 <- which(oos_C$year == 2021 & oos_C$Cstrata == 2)
  wAng21 <- which(oos_C$year == 2021 & oos_C$Cstrata == 3)
  wSA22  <- which(oos_C$year == 2022 & oos_C$Cstrata == 1)
  wIsl22 <- which(oos_C$year == 2022 & oos_C$Cstrata == 2)
  wAng22 <- which(oos_C$year == 2022 & oos_C$Cstrata == 3)
  wSA23  <- which(oos_C$year == 2023 & oos_C$Cstrata == 1)
  wIsl23 <- which(oos_C$year == 2023 & oos_C$Cstrata == 2)
  wAng23 <- which(oos_C$year == 2023 & oos_C$Cstrata == 3)
  wSA24  <- which(oos_C$year == 2024 & oos_C$Cstrata == 1)
  wIsl24 <- which(oos_C$year == 2024 & oos_C$Cstrata == 2)
  wAng24 <- which(oos_C$year == 2024 & oos_C$Cstrata == 3)
  
  # Aggregate simulated catch to annual totals across all monitoring events per year
  pc_6yr <- data.frame(pc19 = rowSums(pC[, w19]), pc20 = rowSums(pC[, w20]),
                       pc21 = rowSums(pC[, w21]), pc22 = rowSums(pC[, w22]),
                       pc23 = rowSums(pC[, w23]), pc24 = rowSums(pC[, w24]))
  
  # Total electrofishing effort (person-seconds) per year
  ef_6yr <- c(sum(oos_C$effort[w19]), sum(oos_C$effort[w20]), sum(oos_C$effort[w21]),
              sum(oos_C$effort[w22]), sum(oos_C$effort[w23]), sum(oos_C$effort[w24]))
  
  c_6yr  <- data.frame(c19 = rowSums(C[, w19]), c20 = rowSums(C[, w20]),
                       c21 = rowSums(C[, w21]), c22 = rowSums(C[, w22]),
                       c23 = rowSums(C[, w23]), c24 = rowSums(C[, w24]))
  
  # Reach-specific catch for diagnostic decomposition
  c_SA  <- data.frame(c19 = rowSums(C[, wSA19]),  c20 = rowSums(C[, wSA20]),
                      c21 = rowSums(C[, wSA21]),  c22 = rowSums(C[, wSA22]),
                      c23 = rowSums(C[, wSA23]),  c24 = rowSums(C[, wSA24]))
  c_Ang <- data.frame(c19 = rowSums(C[, wAng19]), c20 = rowSums(C[, wAng20]),
                      c21 = rowSums(C[, wAng21]), c22 = rowSums(C[, wAng22]),
                      c23 = rowSums(C[, wAng23]), c24 = rowSums(C[, wAng24]))
  c_Isl <- data.frame(c19 = rowSums(C[, wIsl19]), c20 = rowSums(C[, wIsl20]),
                      c21 = rowSums(C[, wIsl21]), c22 = rowSums(C[, wIsl22]),
                      c23 = rowSums(C[, wIsl23]), c24 = rowSums(C[, wIsl24]))
  
  return(list(ef_6yr = ef_6yr, pc_6yr = pc_6yr, c_6yr = c_6yr,
              SpN = SpN, FN = FN, effS = effS,
              jul1_age0 = jul1_age0, tRf_arr = tRf_arr,
              c_SA = c_SA, c_Isl = c_Isl, c_Ang = c_Ang))
}

# =============================================================================
# SECTION 10: COMPUTE OUT-OF-SAMPLE COVARIATES AND RUN FORECASTS
# =============================================================================

# Flow matrices for 2019-2024, April-September (months 4-9).
# Structured as [days x 6 years] matrices matching the format of aQ/sQ
# used in the in-sample analysis. April-September covers the full spawning
# and summer rearing period relevant to covariate calculation.
aQ_out <- angQ_all %>%
  filter(year %in% 2019:2024, month %in% 4:9) %>%
  select(year, cfs, month, day) %>%
  pivot_wider(names_from = year, values_from = cfs) %>%
  select(-month, -day) %>%
  as.matrix

sQ_out <- sanaQ_all %>%
  filter(year %in% 2019:2024, month %in% 4:9) %>%
  select(year, cfs, month, day) %>%
  pivot_wider(names_from = year, values_from = cfs) %>%
  select(-month, -day) %>%
  as.matrix

# ── OOS covariates for Model 1 (original LCC) ────────────────────────────────
# Recompute original LCC covariate for 2019-2024 using the same expert
# elicitation parameters (t1, t2, t3, t4) and calc_cov() as the in-sample fit.
# Expert parameters are identical to preds2 construction above.
preds_out <- array(NA, dim = c(6, 3))
temp <- subset(ee[[3]], ee[[3]][, 1] == "1" & is.na(ee[[3]][, 5]) == FALSE)
t1   <- numeric()
for (t in 1:6) {
  se         <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
  t1[t]      <- temp[t, 6] - se
  se         <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100)
  t1[(t + 7)] <- temp[(t + 6), 6] + se
}
t1[7]  <- 1; t1[14] <- 0
temp <- subset(ee[[3]], ee[[3]][, 1] == "2")
se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
t2   <- temp[1, 6]
temp <- subset(ee[[3]], ee[[3]][, 1] == "4")
se   <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100)
t4   <- temp[1, 6] - se

for (r in 1:3) {
  temp <- subset(ee[[3]], ee[[3]][, 1] == riversegment[r] & is.na(ee[[3]][, 5]) == FALSE)
  t3   <- numeric()
  for (t in 1:6)  {se <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100); t3[t]  <- temp[t, 6] + se}
  for (t in 8:15) {se <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100); t3[t]  <- temp[t, 6] - se}
  t3[7] <- temp[7, 6]
  for (j in 1:6) {
    if (r == 3) {q <- aQ_out[, j]} else {q <- sQ_out[, j]}
    preds_out[j, r] <- calc_cov(q, t3, t2, t4, t1)
  }
}

# Run OOS forecasts for all six models plus the BioOp linear model.
# Each forecast uses mean(covariate from fitting period) to center predictions.

# Model 1: Original LCC (M2_1re)
f1 <- forecast_oos_re(M2_1re, preds_out - mean(preds2))

# Model 2: May-June flow (M2_2re)
# OOS covariate: mean May-June flow for 2019-2024 at each gage / 1000
mjflow_sana_out <- subset(sanaQ_all, (sanaQ_all$month == 5 | sanaQ_all$month == 6) &
                            sanaQ_all$year > 2018)
mjflow_ang_out  <- subset(angQ_all,  (angQ_all$month  == 5 | angQ_all$month  == 6) &
                            angQ_all$year  > 2018)
simpflow_out <- cbind(tapply(mjflow_sana_out$cfs, mjflow_sana_out$year, mean),
                      tapply(mjflow_sana_out$cfs, mjflow_sana_out$year, mean),
                      tapply(mjflow_ang_out$cfs,  mjflow_ang_out$year,  mean)) / 1000

f2 <- forecast_oos_re(M2_2re, simpflow_out - mean(simpflow))

# Models 3-4: Raw and combined inundation covariates
# Loop over 4 habitat metrics; only [[1]] (inund_lower) is used for forecasting.
inund_list_out <- list()
for (i in 1:4) {
  maxmins <- array(NA, dim = c(6, 3))
  for (j in 1:3) {
    if (j == 1) {q <- filter(sanaQ_all, year %in% 2019:2024)} else {q <- filter(angQ_all, year %in% 2019:2024)}
    hab_lookup_ij <- filter(q_hab_lookup, reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup)[i + 2])
    maxmins[, j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
  }
  inund_list_out[[i]] <- maxmins
}

f3 <- forecast_oos_re(m_inund_i_l,
                      log(inund_list_out[[1]] + 1) - mean(log(inund_cov_list[[1]] + 1)))

inund_combined_list_out <- list()
for (i in 1:4) {
  maxmins <- array(NA, dim = c(6, 3))
  for (j in 1:3) {
    if (j == 1) {q <- filter(sanaQ_all, year %in% 2019:2024)} else {q <- filter(angQ_all, year %in% 2019:2024)}
    hab_lookup_ij <- filter(q_hab_lookup_combined, reach == reaches[j]) %>%
      select(cfs, hab = names(q_hab_lookup_combined)[i + 2])
    maxmins[, j] <- left_join(q, hab_lookup_ij) %>%
      mutate(roll_min = rollapplyr(hab, Qdur, min, fill = NA)) %>%
      filter(month %in% c(4:6)) %>%
      group_by(year) %>%
      slice_max(roll_min, with_ties = FALSE) %>%
      pull(roll_min)
  }
  inund_combined_list_out[[i]] <- maxmins
}

f4 <- forecast_oos_re(m_inund_i_l_combined,
                      log(inund_combined_list_out[[1]] + 1) - mean(log(inund_cov_combined_list[[1]] + 1)))

# Models 5-6: Hydraulic LCC covariates for 2019-2024
# Recompute LCC covariate using 2D hydraulic model breakpoints and the same
# expert elicitation parameters (t1, t2, t4) as the in-sample analysis.
lcc_list_out <- list()
for (i in 1:4) {
  preds <- array(NA, dim = c(6, 3))
  # (Expert parameter reconstruction omitted for brevity; identical to Section 5E)
  temp <- subset(ee[[3]], ee[[3]][, 1] == "1" & is.na(ee[[3]][, 5]) == FALSE)
  t1   <- numeric()
  for (t in 1:6) {
    se <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100); t1[t] <- temp[t, 6] - se
    se <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100); t1[(t + 7)] <- temp[(t + 6), 6] + se
  }
  t1[7] <- 1; t1[14] <- 0
  temp <- subset(ee[[3]], ee[[3]][, 1] == "2"); se <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100); t2 <- temp[1, 6]
  temp <- subset(ee[[3]], ee[[3]][, 1] == "4"); se <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100); t4 <- temp[1, 6] - se
  for (r in 1:3) {
    t3 <- filter(q_inund_curve, reach_num == r) %>% pull(i + 1)
    for (t in 1:6) {
      if (r == 3) {q <- aQ_out[, t]} else {q <- sQ_out[, t]}
      preds[t, r] <- calc_cov_2d(q, t3, t2, t4, t1)
    }
  }
  lcc_list_out[[i]] <- preds
}

f5 <- forecast_oos_re(m_lcc_i_l, lcc_list_out[[1]] - mean(inund_lcc_list[[1]]))

lcc_list_out_combined <- list()
for (i in 1:4) {
  preds <- array(NA, dim = c(6, 3))
  temp <- subset(ee[[3]], ee[[3]][, 1] == "1" & is.na(ee[[3]][, 5]) == FALSE)
  t1   <- numeric()
  for (t in 1:6) {
    se <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100); t1[t] <- temp[t, 6] - se
    se <- calcsig(temp[(t + 6), 4], temp[(t + 6), 3], temp[(t + 6), 6], temp[(t + 6), 5] / 100); t1[(t + 7)] <- temp[(t + 6), 6] + se
  }
  t1[7] <- 1; t1[14] <- 0
  temp <- subset(ee[[3]], ee[[3]][, 1] == "2"); se <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100); t2 <- temp[1, 6]
  temp <- subset(ee[[3]], ee[[3]][, 1] == "4"); se <- calcsig(temp[1, 4], temp[1, 3], temp[1, 6], temp[1, 5] / 100); t4 <- temp[1, 6] - se
  for (r in 1:3) {
    t3 <- filter(q_inund_combined_format, reach_num == r) %>% pull(i + 1)
    for (t in 1:6) {
      if (r == 3) {q <- aQ_out[, t]} else {q <- sQ_out[, t]}
      preds[t, r] <- calc_cov_combined(q, t3, t2, t4, t1)
    }
  }
  lcc_list_out_combined[[i]] <- preds
}

f6 <- forecast_oos_re(m_lcc_i_l_combined,
                      lcc_list_out_combined[[1]] - mean(inund_lcc_list_combined[[1]]))

# ── Biological Opinion linear model (reference benchmark) ────────────────────
# The BioOp uses a simple empirical linear model relating May-June mean flow
# at Albuquerque to RGSM catch, developed in the 2016 Biological Opinion.
# Formula: log10(catch + 1) ~ N(-0.1477 + 0.0004*Q - 0.000000014284*Q^2, 0.1)
# Prediction is scaled by electrofishing effort (ef_6yr from model f1).
obs_catch <- tapply(oos_C$hybama, oos_C$year, sum)  # observed annual catch

biop_covs <- angQ_all %>%
  filter(month %in% 5:6, year %in% 2019:2024) %>%
  group_by(year) %>%
  summarize(tbQ = mean(cfs)) %>%  # mean May-June discharge
  mutate(effort = f1$ef_6yr)      # total effort from model f1

biop_preds <- matrix(NA, nrow = 1000000, ncol = 6)
for (i in 1:nrow(biop_covs)) {
  # Simulate 1,000,000 predictions from the BioOp model and scale by effort
  biop_preds[, i] <- (10 ^ (rnorm(1000000,
                                  -0.1477 + (0.0004 * biop_covs$tbQ[i]) -
                                    (biop_covs$tbQ[i] ^ 2) * 0.000000014284,
                                  .1)) - 1) * biop_covs$effort[i] / 100
}


# =============================================================================
# SECTION 11: OUT-OF-SAMPLE COMPARISON FIGURE (Figure 3)
# =============================================================================

# Quantile functions for 80% and 95% credible intervals
q025 <- function(x) {quantile(x, .025)}
q975 <- function(x) {quantile(x, .975)}
q10  <- function(x) {quantile(x, .1)}
q90  <- function(x) {quantile(x, .9)}

# Compile all model predictions into a single list for plotting.
# c_6yr contains the posterior predictive distribution of total annual catch
# (not just predicted means) for each model.
mod_list_oos <- list(
  "LCC IPM"             = f1$c_6yr,
  "May-June IPM"        = f2$c_6yr,
  "Inund. IPM"          = f3$c_6yr,
  "Compos. Inund. IPM"  = f4$c_6yr,
  "Hydraul. LCC"        = f5$c_6yr,
  "Compos. Hydraul. LCC" = f6$c_6yr,
  "BioOp Lin. Mod."     = as.data.frame(biop_preds)
)

mod_out_oos <- data.frame()

for (i in 1:length(mod_list_oos)) {
  output <- mod_list_oos[[i]]
  
  # Summarize the posterior predictive distribution for each year/model:
  # mean, median, and 80% / 95% credible intervals
  df <- data.frame(
    model  = rep(names(mod_list_oos[i]), 6),
    year   = 2019:2024,
    mean   = colMeans(output),
    median = sapply(output, median),
    min    = sapply(output, min),
    l025   = sapply(output, q025),
    l10    = sapply(output, q10),
    u90    = sapply(output, q90),
    u975   = sapply(output, q975))
  
  mod_out_oos <- bind_rows(mod_out_oos, df)
}

obs_catch_df <- data.frame(year = 2019:2024, catch = obs_catch)

# Join observed catch and set factor order for plotting (simpler models first,
# then hydraulic models, to facilitate visual comparison)
mod_out_catch <- left_join(mod_out_oos, obs_catch_df) %>%
  mutate(model = factor(model,
                        levels = c("LCC IPM", "May-June IPM", "BioOp Lin. Mod.",
                                   "Hydraul. LCC", "Compos. Hydraul. LCC",
                                   "Inund. IPM", "Compos. Inund. IPM")))

# Figure 3: Faceted comparison of predicted vs. observed annual RGSM catch.
# Each facet = one year. For each model:
#   Black point = posterior mean
#   Thin red bars = 95% credible interval
#   Thick black bars = 80% credible interval
#   Dashed line = observed total annual catch
# The light blue rectangle highlights the four hydraulic/combined models,
# distinguishing them from the three simpler reference models.
mod_out_catch %>%
  ggplot(aes(x = model)) +
  geom_point(aes(y = mean), color = "black") +
  geom_errorbar(aes(ymin = l025, ymax = u975), width = 0, color = "red") +   # 95% CI
  geom_errorbar(aes(ymin = l10, ymax = u90), width = 0, color = "black") +   # 80% CI
  geom_hline(data = obs_catch_df, aes(yintercept = catch), lty = 2) +        # observed
  geom_rect(aes(xmin = 3.5, xmax = 7.5, ymin = 0, ymax = Inf),              # highlight box
            fill = "lightblue", alpha = .1) +
  facet_wrap(~ year, scales = "free_y", nrow = 2) +
  labs(x = "Model", y = "RGSM Catch") +
  theme_ipsum(axis_title_just  = 1,
              base_size         = 14,
              strip_text_size   = 18,
              axis_title_size   = 16) +
  theme(axis.text.x    = element_text(angle = 90),
        panel.spacing.y = unit(0.1, "lines"))

ggsave(filename = "plots/fig3_oos_results_new_mod_order.jpeg",
       bg = "white", width = 8, height = 7, units = "in")
