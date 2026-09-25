library(shiny)
library(bslib)
library(leaflet)
library(htmltools)
library(sf)
library(tidyverse)
library(zoo)
library(patchwork)
library(hrbrthemes)
library(scales)

# New functions ----
## Helper functions ----
calcsig <- function(up, down, mean, q, INT = c(0, 1000)) {
  sig <- function(x) {
    abs(pnorm(up, mean, x) - pnorm(down, mean, x) - q)
  }
  optimize(sig, interval = INT)$minimum
}

lin_int <- function(x, xlow, xhi, ylow, yhi) {
  (yhi - ylow) * (x - xlow) / (xhi - xlow) + ylow
}

q025 <- function(x) {quantile(x, .025)}
q975 <- function(x) {quantile(x, .975)}
q10  <- function(x) {quantile(x, .1)}
q90  <- function(x) {quantile(x, .9)}

## Selecting hydrographs for simulation ----
q_sim <- function(flow_input) {
  q_all <- array(NA, dim = c(214, length(flow_input), 2))
  
  if (length(flow_input[flow_input < 1]) > 0) {
    for (i in 1:length(flow_input)) {
      aQ_quant  <- quantile(aQ_med_spr$tot_vol, flow_input[i])
      a_year_sel <- aQ_med_spr$year[which.min(abs(aQ_quant - aQ_med_spr$tot_vol))]
      
      sQ_quant  <- quantile(sQ_med_spr$tot_vol, flow_input[i])
      s_year_sel <- sQ_med_spr$year[which.min(abs(sQ_quant - sQ_med_spr$tot_vol))]
      
      aQ_temp <- filter(angQ,  year == a_year_sel, month > 2, month < 10)
      sQ_temp <- filter(sanaQ, year == s_year_sel, month > 2, month < 10)
      
      q_all[, i, 1] <- aQ_temp$cfs
      q_all[, i, 2] <- sQ_temp$cfs
    }
  } else {
    for (i in 1:length(flow_input)) {
      aQ_temp <- filter(angQ,  year == flow_input[i], month > 2, month < 10)
      sQ_temp <- filter(sanaQ, year == flow_input[i], month > 2, month < 10)
      
      q_all[, i, 1] <- aQ_temp$cfs
      q_all[, i, 2] <- sQ_temp$cfs
    }
  }
  return(q_all)
}

## Functions for calculating the larval carrying capacity index covariate ----
q2hab_2d <- function(q, pars, input) {
  p  <- c(pars, pars[length(pars)])
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

preds_calc <- function(exp, flows, hab_lookup) {
  preds <- array(NA, dim = c(dim(flows)[2], 3))
  
  temp <- subset(ee[[exp]], ee[[exp]][, 1] == "1" & is.na(ee[[exp]][, 5]) == FALSE)
  t1 <- numeric()
  for (t in 1:6) {
    se    <- calcsig(temp[t, 4], temp[t, 3], temp[t, 6], temp[t, 5] / 100)
    t1[t] <- temp[t, 6] - se
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
    t3 <- filter(hab_lookup, reach_num == r)
    for (j in 1:dim(flows)[2]) {
      if (r == 3) {q <- flows[, j, 1]} else {q <- flows[, j, 2]}
      preds[j, r] <- calc_cov_combined(q, t3, t2, t4, t1)
    }
  }
  
  return(preds)
}

## Functions for calculating inundation covariate ----
hab_lookup_fun <- function(inund_curve, reaches) {
  hab_lookup <- tibble()
  
  for (i in 1:length(reaches)) {
    q_hab_i <- filter(inund_curve, reach == reaches[i])
    sites   <- unique(q_hab_i$site)
    hab_lookup_i <- tibble()
    
    for (j in 1:length(sites)) {
      q_hab_ij <- if ("site" %in% colnames(q_hab_i)) {
        q_hab_i %>% filter(site %in% sites[j])
      } else {
        q_hab_i
      }
      
      hab_lookup_ij <- tibble(cfs = 1:10e3) %>%
        mutate(hab = q2hab_2d(cfs, q_hab_ij$hab, q_hab_ij$cfs))
      
      hab_lookup_i <- bind_rows(hab_lookup_i, hab_lookup_ij)
    }
    
    hab_lookup_i <- group_by(hab_lookup_i, cfs) %>%
      summarize(hab = sum(hab)) %>%
      mutate(reach = reaches[i])
    
    hab_lookup <- bind_rows(hab_lookup, hab_lookup_i)
  }
  return(hab_lookup)
}

maxmins_fun <- function(lookup, sim_q, reaches) {
  maxmins_hab <- array(NA, dim = c(dim(sim_q)[2], 3))
  maxmins_cfs <- array(NA, dim = c(dim(sim_q)[2], 3))
  
  for (i in 1:length(reaches)) {
    hab_lookup_i <- filter(lookup, reach == reaches[i])
    
    if (i == 1) {
      q <- round(sim_q[, , 2], 0) %>% as_tibble() %>% set_names("cfs")
    } else {
      q <- round(sim_q[, , 1], 0) %>% as_tibble() %>% set_names("cfs")
    }
    
    for (t in 1:dim(sim_q)[2]) {
      x <- q[, t]
      maxmins_hab[t, i] <- max(rollapplyr(hab_lookup_i$hab[match(x$cfs, hab_lookup_i$cfs)],
                                          Qdur, mean, fill = NA)[32:122], na.rm = T)
      maxmins_cfs[t, i] <- max(rollapplyr(x$cfs,
                                          Qdur, mean, fill = NA)[32:122], na.rm = T)
    }
  }
  
  return(list(maxmins_hab = maxmins_hab, maxmins_cfs = maxmins_cfs))
}

## Flow augmentation function ----
aug_q_fun <- function(st_date, en_date, vol_AcFt, q_array) {
  simq_dates <- data.frame(
    date    = seq(as.Date("2000-03-01"), as.Date("2000-09-30"), 1),
    day_ind = 1:214
  ) %>%
    mutate(date = as.character(format(date, "%m-%d")))
  
  start_ind  <- simq_dates$day_ind[simq_dates$date == st_date]
  end_ind    <- simq_dates$day_ind[simq_dates$date == en_date]
  aug_vec    <- rep(0, 214)
  ndays_aug  <- length(start_ind:end_ind)
  augCFS     <- round(vol_AcFt * 43559.935 / ndays_aug / (24 * 60 * 60), 0)
  
  aug_vec[start_ind:end_ind] <- rep(augCFS, ndays_aug)
  aug_q_array <- q_array + aug_vec
  return(aug_q_array)
}

## Forecast recruit function ----
forecast_recruit <- function(mod,
                             Xout,
                             years_sim,
                             maxminsBase,
                             maxminsAug,
                             maxminsRest,
                             reaches) {
  a          <- pull(mod, a)
  mu_lbeta   <- select(mod, starts_with("mu_lbeta")) %>% as.matrix()
  sd_lbeta   <- pull(mod, sd_lbeta)
  B_lbeta    <- pull(mod, B_lbeta)
  effS_all   <- select(mod, starts_with("effS"))
  iter       <- length(a)
  
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
  
  age0_base  <- array(NA, dim = c(iter, nrow(Xout), 3))
  age0_rest  <- array(NA, dim = c(iter, nrow(Xout), 3))
  age0_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  diff_rest  <- array(NA, dim = c(iter, nrow(Xout), 3))
  prop_rest  <- array(NA, dim = c(iter, nrow(Xout), 3))
  diff_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  prop_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  ratio_rest <- array(NA, dim = c(iter, nrow(Xout), 3))
  ratio_hydro <- array(NA, dim = c(iter, nrow(Xout), 3))
  effS_sim   <- array(NA, dim = c(iter, 3))
  
  for (k in 1:3) {
    effS <- filter(effS_long, year == effS_ind[k], reach == k) %>% pull(effS)
    
    for (i in 1:nrow(Xout)) {
      tRf0       <- exp(rnorm(iter, mu_lbeta[, k], sd_lbeta))
      tRf_base   <- exp(B_lbeta * Xout[i, k, 1]) * tRf0
      tRf_rest   <- exp(B_lbeta * Xout[i, k, 2]) * tRf0
      tRf_hydro  <- exp(B_lbeta * Xout[i, k, 3]) * tRf0
      
      age0_base[, i, k]  <- a * (effS) / (1 + a * (effS) / tRf_base)
      age0_rest[, i, k]  <- a * (effS) / (1 + a * (effS) / tRf_rest)
      age0_hydro[, i, k] <- a * (effS) / (1 + a * (effS) / tRf_hydro)
      diff_rest[, i, k]  <- age0_rest[, i, k]  - age0_base[, i, k]
      diff_hydro[, i, k] <- age0_hydro[, i, k] - age0_base[, i, k]
      prop_hydro[, i, k] <- diff_hydro[, i, k] / age0_base[, i, k]
      prop_rest[, i, k]  <- diff_rest[, i, k]  / age0_base[, i, k]
      ratio_hydro[, i, k] <- age0_hydro[, i, k] / age0_base[, i, k]
      ratio_rest[, i, k]  <- age0_rest[, i, k]  / age0_base[, i, k]
    }
    
    effS_sim[, k] <- effS
  }
  
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
  
  age0_dat <- tibble(
    reach                = rep(reaches, each = nrow(Xout)),
    flow_years           = rep(years_sim, 3),
    cfs_baseline         = as.vector(maxminsBase$maxmins_cfs),
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
    
    ratio_rest_med       = ratio_rest_med,
    lower_ratio_rest     = lower_ratio_rest,
    upper_ratio_rest     = upper_ratio_rest,
    
    ratio_hydro_med      = ratio_hydro_med,
    lower_ratio_hydro    = lower_ratio_hydro,
    upper_ratio_hydro    = upper_ratio_hydro
  )
  
  return(age0_dat)
}

## Plotting functions ----
recruit_plot <- function(recruit_dat, rest_reach, rest_acres, aug_vol) {
  
  col_values_1 <- c("black", "red", "blue")
  names(col_values_1) <- c(
    "Baseline",
    paste0("Floodplain\nrestoration\n(", round(rest_acres, 0), " acres)"),
    paste0("Add\n", format(aug_vol, big.mark = ","), "\nacre-ft")
  )
  
  p1 <-
    recruit_dat %>%
    filter(reach %in% rest_reach) %>%
    ggplot(aes(x = cfs_baseline)) +
    geom_point(aes(y = log(recruits_baseline), color = names(col_values_1)[1]),
               cex = 2.5, stroke = 2, pch = 21, position = position_nudge(x = -15)) +
    geom_point(aes(y = log(recruits_rest),     color = names(col_values_1)[2]),
               size = 2.5, stroke = 2, pch = 22, position = position_nudge(x = 15)) +
    geom_point(aes(y = log(recruits_hydro),    color = names(col_values_1)[3]),
               size = 2.5, stroke = 2, pch = 24) +
    facet_wrap(~ reach) +
    labs(x = "Flow index", y = "Age0 recruits (log-scale)") +
    scale_color_manual(name = "Scenarios",
                       breaks = names(col_values_1),
                       values = col_values_1) +
    theme_ipsum() +
    theme(plot.margin   = unit(c(1, 0, 1, 1), "cm"),
          axis.title.x  = element_text(size = 14, hjust = 0.5),
          axis.title.y  = element_text(size = 16, hjust = 0.5),
          plot.title    = element_text(size = 16, face = "plain"),
          legend.title  = element_text(size = 16),
          legend.text   = element_text(size = 14))
  
  nudge_size_b <- max(c(max(recruit_dat$diff_recruits_hydro),
                        max(recruit_dat$diff_recruits_rest))) / (1000 * 7)
  
  col_values_2 <- c("red", "blue")
  names(col_values_2) <- c(
    paste0("Floodplain\nrestoration\n(", round(rest_acres, 0), " acres)"),
    paste0("Add\n", format(aug_vol, big.mark = ","), "\nacre-ft")
  )
  
  p2 <-
    recruit_dat %>%
    filter(reach %in% rest_reach) %>%
    ggplot(aes(x = cfs_baseline)) +
    geom_point(aes(y = diff_recruits_rest  / 1000, color = names(col_values_2)[1]),
               size = 2.5, stroke = 2, pch = 22) +
    geom_point(aes(y = diff_recruits_hydro / 1000, color = names(col_values_2)[2]),
               size = 2.5, stroke = 2, pch = 24) +
    geom_text(aes(y = nudge_size_b + (diff_recruits_hydro / 1000), label = flow_years)) +
    facet_wrap(~ reach) +
    labs(x = "Flow index", y = "Recruits over baseline (in thousands)") +
    scale_color_manual(name = "Scenarios",
                       breaks = names(col_values_2),
                       values = col_values_2) +
    theme_ipsum() +
    theme(plot.margin          = unit(c(1, 0, 1, 1), "cm"),
          axis.title.x         = element_text(size = 14, hjust = 0.5),
          axis.title.y         = element_text(size = 16, hjust = 0.5),
          plot.title           = element_text(size = 16, face = "plain"),
          legend.title         = element_text(size = 16),
          legend.text          = element_text(size = 14),
          legend.key.spacing.y = unit(.5, "lines"))
  
  p2 / p1 +
    plot_annotation(title = rest_reach) &
    theme(strip.text = element_blank())
}

flow_plot <- function(flow_input, aug_q, reach) {
  sim_q <- q_sim(flow_input)
  
  if (reach == "San Acacia") {
    q          <- round(sim_q[, , 2], 0)
    aug_reach  <- round(aug_q[, , 2], 0)
    med_flow   <- sQ_med_dly
  } else {
    q          <- round(sim_q[, , 1], 0)
    aug_reach  <- round(aug_q[, , 1], 0)
    med_flow   <- aQ_med_dly
  }
  
  aug_df <- as_tibble(aug_reach) %>%
    set_names(flow_input) %>%
    bind_cols(med_flow[, c("month", "day")]) %>%
    pivot_longer(cols = as.character(flow_input),
                 names_to  = "input_val",
                 values_to = "aug_cfs")
  
  q %>%
    as_tibble() %>%
    set_names(flow_input) %>%
    bind_cols(med_flow) %>%
    pivot_longer(cols = as.character(flow_input),
                 names_to  = "input_val",
                 values_to = "cfs") %>%
    left_join(aug_df) %>%
    arrange(input_val) %>%
    mutate(day = as.Date(paste(month, day, sep = "-"), format = "%m-%d")) %>%
    filter(month < 8) %>%
    ggplot() +
    geom_line(aes(day, aug_cfs, group = input_val, color = "E-flow",  linetype = "E-flow")) +
    geom_line(aes(day, cfs,     group = input_val, color = "Flow",    linetype = "Flow")) +
    geom_line(aes(day, med_cfs, group = input_val, color = "Median",  linetype = "Median")) +
    scale_color_manual(breaks = c("Flow", "Median", "E-flow"),
                       values = c("Flow" = "salmon", "Median" = "black", "E-flow" = "blue")) +
    scale_linetype_manual(breaks = c("Flow", "Median", "E-flow"),
                          values = c("Flow" = "solid", "Median" = "dashed", "E-flow" = "solid")) +
    labs(color = "", linetype = "", y = "CFS", x = "Date") +
    facet_wrap(~ input_val, nrow = length(flow_input), scales = "free") +
    ggtitle(reach) +
    theme_ipsum() +
    theme(plot.margin  = unit(c(0.5, 0, 1, 1), "cm"),
          axis.title.x = element_text(size = 14, hjust = 0.5),
          axis.title.y = element_text(size = 16, hjust = 0.5),
          plot.title   = element_text(size = 16, face = "plain"),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 12))
}

flow_rest_plot <- function(flow_input, reach_sel, aug_q, lookup, lookup_rest) {
  sim_q <- q_sim(flow_input)
  
  hab_lookup <- lookup_rest %>%
    rename(rest_hab = hab) %>%
    left_join(lookup) %>%
    mutate(tot_hab = hab + rest_hab) %>%
    filter(reach == reach_sel)
  
  if (reach_sel == "San Acacia") {
    q          <- round(sim_q[, , 2], 0)
    aug_reach  <- round(aug_q[, , 2], 0)
    med_flow   <- sQ_med_dly
  } else {
    q          <- round(sim_q[, , 1], 0)
    aug_reach  <- round(aug_q[, , 1], 0)
    med_flow   <- aQ_med_dly
  }
  
  aug_df <- as_tibble(aug_reach) %>%
    set_names(flow_input) %>%
    bind_cols(med_flow[, c("month", "day")]) %>%
    pivot_longer(cols = as.character(flow_input),
                 names_to  = "input_val",
                 values_to = "cfs")
  
  qhab1 <- as_tibble(q) %>%
    set_names(flow_input) %>%
    bind_cols(med_flow) %>%
    mutate(date = as.Date(paste(month, day, sep = "-"), format = "%m-%d")) %>%
    pivot_longer(cols = -c("date", "month", "day", "med_cfs"),
                 names_to  = "input_val",
                 values_to = "cfs") %>%
    left_join(hab_lookup)
  
  qhab2 <- select(qhab1, date, month, day, input_val) %>%
    left_join(aug_df) %>%
    left_join(hab_lookup[, c("cfs", "hab")]) %>%
    rename(aug_cfs = cfs, aug_hab = hab)
  
  qhab_all <- left_join(qhab1, qhab2) %>%
    filter(month < 8) %>%
    mutate(hab_diff_aug = aug_hab - hab)
  
  p2 <- ggplot(qhab_all) +
    geom_line(aes(date, tot_hab, group = input_val, color = "With\nRestoration")) +
    geom_line(aes(date, aug_hab, group = input_val, color = "Augmented\nFlow")) +
    geom_line(aes(date, hab,     group = input_val, color = "Baseline")) +
    scale_color_manual(breaks = c("With\nRestoration", "Baseline", "Augmented\nFlow"),
                       values = c("With\nRestoration" = "salmon", "Baseline" = "black",
                                  "Augmented\nFlow" = "blue")) +
    labs(color = "", y = "Acres", x = "Date") +
    facet_wrap(~ input_val, nrow = length(flow_input), scales = "free_y") +
    ggtitle("Inundated habitat") +
    theme_ipsum() +
    theme(plot.margin  = unit(c(0.5, 1, 1, 0), "cm"),
          axis.title.x = element_text(size = 14, hjust = 0.5),
          axis.title.y = element_text(size = 16, hjust = 0.5),
          plot.title   = element_text(size = 16, face = "plain"),
          legend.title = element_text(size = 14),
          legend.text  = element_text(size = 12))
  
  p3 <- ggplot(qhab_all) +
    geom_line(aes(date, rest_hab,     group = input_val, color = "With\nRestoration",
                  linetype = "With\nRestoration")) +
    geom_line(aes(date, hab_diff_aug, group = input_val, color = "Augmented\nFlow",
                  linetype = "Augmented\nFlow")) +
    scale_color_manual(breaks = c("With\nRestoration", "Augmented\nFlow"),
                       values = c("With\nRestoration" = "salmon", "Augmented\nFlow" = "blue")) +
    scale_linetype_manual(breaks = c("With\nRestoration", "Augmented\nFlow"),
                          values = c("With\nRestoration" = "solid", "Augmented\nFlow" = "solid")) +
    labs(color = "", linetype = "", y = "Acres", x = "Date") +
    facet_wrap(~ input_val, nrow = length(flow_input), scales = "free_y") +
    ggtitle("Habitat above baseline") +
    theme_ipsum() +
    theme(plot.margin   = unit(c(0.5, 1, 1, 0), "cm"),
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(size = 16, hjust = 0.5),
          plot.title    = element_text(size = 16, face = "plain"),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 12),
          legend.position = "none")
  
  p2 + p3 + plot_annotation(title = reach_sel) &
    theme(plot.title = element_text(face = "bold", size = 16))
}

# Select model attributes ---------------------------------------------------
spawn_rank <- 9
mod_name   <- "lcc_low_combined"
hab_cov    <- "inund_lower"
Qdur       <- 20

sim_mod <- read_csv(paste0("output/", mod_name, "_mcmc_sub.csv"))
inund_lcc_list_combined <- readRDS("output/inund_lcc_list_combined.RData")
cov_list <- inund_lcc_list_combined[[1]]

# Prepare data ----------------------------------------------------------
## Reach and map data ---
reaches <- c("San Acacia", "Isleta", "Angostura")

flow_perm <- st_read("data/flow_perm/flow_perm.shp") %>%
  filter(month == "July") %>%
  slice(-1)

line_segments <- lapply(1:(nrow(flow_perm) - 1), function(i) {
  coords <- st_coordinates(flow_perm)[i:(i + 1), ]
  line   <- st_linestring(coords)
  st_sf(geometry  = st_sfc(line, crs = 4326),
        flow_perm = flow_perm$value[i])
})
lines_sf <- do.call(rbind, line_segments) %>%
  mutate(flow_perm = round(flow_perm, 2))

## Flow data ----
angQ  <- read.csv("data/yackulic2022_data/abq_gage_08330000.csv") %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))
sanaQ <- read.csv("data/yackulic2022_data/SanAcacia_gage_08354900.csv") %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

aQ_med_dly <- angQ %>%
  filter(month > 2, month < 10) %>%
  group_by(month, day) %>%
  summarise(med_cfs = median(cfs))

sQ_med_dly <- sanaQ %>%
  filter(month > 2, month < 10) %>%
  group_by(month, day) %>%
  summarise(med_cfs = median(cfs))

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

ee          <- readRDS("output/ee_list.RData")
expert_num  <- 3

## Hydraulic data ----
q_inund_curve <- read_csv("data/ecoval2d.csv") %>%
  select(cfs = q, inund_lower = tot_outside_chan_lower,
         inund_upper = tot_outside_chan_upper,
         wua_lower, wua_upper, reach, reach_num) %>%
  pivot_longer(cols = inund_lower:wua_upper,
               names_to  = "hab_type",
               values_to = "hab")

q_hab_lookup_combined <- read_csv("data/q_hab_lookup_combined.csv")

reach_df <- data.frame(reach = reaches, reach_num = c(1, 2, 3))

no_flow_hab <- reach_df %>% mutate(hab = 0, cfs = 0)

hab_lookup_baseline <- q_hab_lookup_combined %>%
  select(cfs, reach, hab = any_of(hab_cov)) %>%
  left_join(reach_df) %>%
  bind_rows(no_flow_hab) %>%
  arrange(reach, cfs)

## Restoration lookup (used for generic acre-based scaling) ----
q_hab_rest_raw <- read_csv("data/san_acacia_rest_ecovals.csv") %>%
  filter(terrain == "AB") %>%
  select(-terr_inund_notes, raw_hab = hab)

terr_ratio <- q_hab_rest_raw %>%
  filter(cfs == 7000) %>%
  mutate(ratio = raw_hab / terr_ac) %>%
  select(site, ratio)

q_hab_rest_norm <- q_hab_rest_raw %>%
  left_join(terr_ratio) %>%
  mutate(norm_terr_hab = case_when(ratio == Inf ~ terr_ac * 1,
                                   TRUE ~ terr_ac * ratio),
         hab = raw_hab - norm_terr_hab)

q_hab_rest <- select(q_hab_rest_norm, cfs, reach, hab, site)

sana_site_props_2 <- q_hab_rest_raw %>%
  filter(cfs == 7000) %>%
  select(site, size, raw_hab) %>%
  mutate(ratio = raw_hab / 6284)

no_rest_reaches    <- reaches[c(2, 3)]
q_hab_curves_non_rest <- tibble(cfs   = rep(c(0, 7000), length(no_rest_reaches)),
                                reach = rep(no_rest_reaches, each = 2),
                                hab   = 0)

hab_lookup_rest_all <- tibble()
rest_sites   <- unique(q_hab_rest$site)
restReach    <- "San Acacia"

init_inund <- q_hab_rest_raw %>%
  select(site, cfs, raw_hab) %>%
  filter(cfs != 0) %>%
  group_by(site) %>%
  slice_min(raw_hab) %>%
  slice_max(cfs) %>%
  select(site, cfs)

for (i in 1:length(rest_sites)) {
  reach_prop_i  <- filter(sana_site_props_2, site == rest_sites[i]) %>% pull(ratio)
  init_inund_i  <- filter(init_inund, site == rest_sites[i]) %>% pull()
  
  baseline_hab_prop_i <- hab_lookup_baseline %>%
    mutate(hab_site_base = hab * reach_prop_i)
  
  q_hab_i  <- filter(q_hab_rest_raw, site == rest_sites[i]) %>%
    select(cfs, reach, hab = raw_hab, site)
  curves_i <- bind_rows(q_hab_curves_non_rest, q_hab_i)
  
  # Matches analysis script: does NOT add channel habitat; clamps negatives to zero
  lookup_i <- hab_lookup_fun(curves_i, reaches) %>%
    mutate(site = rest_sites[i]) %>%
    rename(raw_hab = hab) %>%
    left_join(baseline_hab_prop_i) %>%
    mutate(rest_hab_norm = raw_hab - hab_site_base) %>%
    mutate(rest_hab_norm = case_when(rest_hab_norm < 0 ~ 0,
                                     TRUE ~ rest_hab_norm))
  
  hab_lookup_rest_all <- bind_rows(hab_lookup_rest_all, lookup_i)
}

hab_lookup_rest_global <- hab_lookup_rest_all %>%
  group_by(cfs, reach) %>%
  summarize(hab = sum(rest_hab_norm)) %>%
  left_join(reach_df) %>%
  bind_rows(no_flow_hab) %>%
  arrange(reach, cfs) %>%
  mutate(hab = case_when(!reach %in% restReach ~ 0, TRUE ~ hab))

floodplain_site_tot <- sum(sana_site_props_2$raw_hab)

lookup_per_acre <- hab_lookup_rest_global %>%
  mutate(hab = hab / floodplain_site_tot)

lookup_per_acre_rest_reach <- lookup_per_acre %>%
  filter(reach == "San Acacia") %>%
  select(cfs, hab)

hab_lookup_rest_total <- bind_rows(hab_lookup_baseline, hab_lookup_rest_global) %>%
  group_by(reach, cfs, reach_num) %>%
  summarize(hab = sum(hab))

# UI function ----------------------------------------------------------
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      
      h4("Hydrology inputs", style = "font-style: italic;"),
      
      textInput("integerInput", "Years or numbers from 0 - 1",
                value = "2003, 2005, 2006, 2007, 2008"),
      helpText("Enter years between 1993 and 2020 separated by commas, or enter 
            values between 0 and 1 to select flow quantiles (e.g., 0.1 = dry 
            year, 0.9 = wet year). For site comparison, using all flow years 
            is strongly recommended."),
      
      checkboxInput("usePreset", "Use all flow years", value = FALSE),
      
      radioButtons("reachInput", "Reach of interest", reaches),
      helpText("Select the reach where your restoration site is located."),
      
      h4("Environmental flow inputs", style = "font-style: italic;"),
      
      numericInput("augInput", "E-flow volume (ac-ft)", value = 10000),
      helpText("Set to 0 if you are evaluating the restoration scenario only 
            and do not want to include a flow augmentation scenario."),
      
      h5("E-flow dates", style = "font-style: italic;"),
      fluidRow(
        column(6, textInput("startDate", "Start (MM-DD):", value = "05-20")),
        column(6, textInput("endDate",   "End (MM-DD):",   value = "06-14"))
      ),
      helpText("The specified volume will be divided evenly across all days 
            in this date range and added to the simulated hydrograph."),
      
      h4("Restoration inputs", style = "font-style: italic;"),
      
      radioButtons(
        "restInputMode",
        "Restoration input method",
        choices  = c("Generic (scale by acres)" = "generic",
                     "Custom flow-habitat curve" = "custom"),
        selected = "generic"
      ),
      
      conditionalPanel(
        condition = "input.restInputMode == 'generic'",
        numericInput("restAcres", "Restoration area (acres)", value = 50),
        helpText("If you do not have site-specific hydraulic data, enter your 
              best estimate of the restorable area in acres. Results can be 
              re-run with a custom curve once hydraulic modeling is complete.")
      ),
      
      conditionalPanel(
        condition = "input.restInputMode == 'custom'",
        fileInput(
          "restCurveFile",
          label   = "Upload CSV with columns 'cfs' and 'hab'",
          accept  = c(".csv", "text/csv")
        ),
        helpText("'hab' should be inundated habitat in acres at your restoration 
              site for each discharge value in 'cfs'. Values will be linearly 
              interpolated across the 1–7,000 cfs range. Flows outside the 
              uploaded range will be clamped to the nearest endpoint value.")
      ),
      
      h5("Restoration location", style = "font-style: italic;"),
      helpText("Optional. Enter decimal-degree coordinates to visualize the 
            location of your restoration site relative to summer flow 
            permanence on the Map tab."),
      fluidRow(
        column(6, numericInput("lat", "Latitude:",  value = 34.225, step = 0.01)),
        column(6, numericInput("lng", "Longitude:", value = -106.899, step = 0.01))
      ),
      
      downloadButton("downloadData", "Download output"),
      helpText("Downloads a CSV with one row per flow year containing recruitment 
            estimates under baseline, restoration, and flow augmentation 
            scenarios. Run and download separately for each site, renaming 
            the file each time (e.g., site_A_output.csv).")
    ),
    
    mainPanel(
      style = "height: 1200px; overflow-y: auto;",
      tabsetPanel(
        tabPanel(
          "Home",
          h3("Rio Grande Silvery Minnow Recruitment Tool"),
          
          h4("What does this tool do?"),
          p("This tool estimates the number of Rio Grande Silvery Minnow (RGSM) age-0 recruits 
    that would be produced under different floodplain restoration and flow management 
    scenarios at sites along the Middle Rio Grande, New Mexico. It is designed for 
    restoration practitioners, agency planners, and resource managers who want to compare 
    the expected recruitment benefits of different restoration sites or management actions 
    before committing to detailed design or planning work."),
          p("A typical workflow involves selecting a set of historical flow years to simulate, 
    specifying one or more restoration sites using either an estimated acreage or an 
    uploaded flow-habitat curve, reviewing the Habitat and Recruitment tabs to confirm 
    that results look reasonable, and downloading the output CSV for each site. Downloaded 
    files can then be combined in a spreadsheet to rank or compare sites based on their 
    projected recruitment benefit across the range of historical flow conditions."),
          p("For a detailed step-by-step guide to this workflow, see the ",
            a(strong("User Walkthrough"), 
              href   = "https://usace-wrises.github.io/rgsm_hydraulic/walkthrough",
              target = "_blank",
              rel    = "noopener noreferrer"),
            " linked from this tab."),
          
          h4("Instructions"),
          p("Adjust the inputs in the sidebar to configure your scenario. 
   The plots and download output automatically update based on your selections."),
          tags$ul(
            tags$li(
              strong("Flow years"),
              tags$ul(
                tags$li("Enter specific years between 1993 and 2020, separated by commas 
               (e.g., 2003, 2005, 2010). The tool will simulate hydrology from 
               each of those years and show recruitment estimates for each."),
                tags$li("You can also enter values between 0 and 1 to select flow quantiles 
               rather than specific years. For example, entering 0.1 selects a year 
               representative of the driest 10% of the record, while 0.9 selects a 
               year representative of the wettest 10%. This is useful for 
               understanding how a site performs across a stylized range of 
               conditions rather than specific historical years."),
                tags$li(strong("For site comparison, checking 'Use all flow years' is 
               strongly recommended."), " This runs the full 1993–2020 historical 
               record and captures a large range of wet and dry conditions. 
               Comparing sites using only a few years risks selecting years that 
               favor one site over another for reasons unrelated to site quality. 
               Averaging or summarizing the downloaded output across all years 
               gives the most robust basis for comparison.")
              )
            ),
            tags$li(
              strong("Reach of interest"),
              tags$ul(
                tags$li("Select the reach where your restoration site is located. This 
               controls which reach is highlighted in the Recruitment and Habitat 
               plots.")
              )
            ),
            tags$li(
              strong("Environmental flow inputs"),
              tags$ul(
                tags$li("Enter the volume of water in acre-feet that would be added as an 
               environmental flow release. This volume is distributed evenly across 
               the date range you specify."),
                tags$li("If you are interested only in the restoration scenario and not in 
               evaluating a flow augmentation scenario, set the volume to 0. The 
               baseline and restoration results are not affected by this input."),
                tags$li("The start and end dates control the window over which the 
               environmental flow is added to the hydrograph. Earlier windows 
               (e.g., May) tend to coincide with peak spawning and generally 
               produce larger recruitment responses than later windows, though 
               this varies by year.")
              )
            ),
            tags$li(
              strong("Restoration inputs"),
              tags$ul(
                tags$li(HTML("<b>Generic (scale by acres):</b> Use this option when you have 
               an estimate of the total restorable area at your site but have not 
               yet conducted site-specific hydraulic modeling. The tool scales a 
               representative per-acre flow-inundation curve by the acreage you 
               enter. This is appropriate for early-stage screening and site 
               prioritization.")),
                tags$li(HTML("<b>Custom flow-habitat curve:</b> Use this option when you 
               have site-specific hydraulic modeling results. Upload a CSV with 
               two columns — <code>cfs</code> (discharge) and <code>hab</code> 
               (inundated habitat in acres) — and the tool will use that 
               relationship directly. This produces more accurate recruitment 
               projections and is recommended for later-stage analyses and 
               formal planning documents.")),
                tags$li("You can re-run the tool for the same site using first the generic 
               option and then a custom curve to see how sensitive the results are 
               to the assumed flow-habitat relationship.")
              )
            ),
            tags$li(
              strong("Downloading and comparing results"),
              tags$ul(
                tags$li("Click 'Download output' to save a CSV for the current site 
               configuration. Rename each file before running the next site so 
               you do not overwrite it (e.g., site_A_output.csv, site_B_output.csv)."),
                tags$li(HTML("To compare sites, open the CSV files in Excel or R and focus 
               on two columns: <code>diff_recruits_rest</code> (additional recruits 
               above baseline) and <code>ratio_rest_med</code> (ratio of restoration 
               to baseline recruitment). See the 'How to interpret results' section 
               below for guidance on which metric to use."))
              )
            )
          ),
          h4("How to interpret results"),
          p("The following points will help you make sense of what you see in the plots 
   and in the downloaded data."),
          tags$ul(
            tags$li(
              strong("The x-axis is a flow index, not a calendar year."),
              " The x-axis in the Recruitment and Habitat plots shows the seasonal maximum 
    of the 20-day rolling minimum discharge for that year, a flow metric designed to indicate whether 
    there were high flows of suitable duration (i.e., 20 days) for adequate RGSM spawning. Higher values indicate 
    wetter years with more sustained high flows. Points are labeled with the 
    corresponding year so you can identify specific years of interest."
            ),
            tags$li(
              strong("Absolute recruitment gains are larger in wet years."),
              " In wet years, both baseline recruitment and the restoration benefit are 
    larger in absolute terms because more of the floodplain is inundated for 
    longer. This does not necessarily mean that a restoration site is more 
    valuable in wet years — it may simply reflect the higher baseline. The ratio 
    column (", code("ratio_rest_med"), ") adjusts for this by expressing the 
    restoration benefit relative to what would have happened without restoration."
            ),
            tags$li(
              strong("Use the difference column to estimate total fish production; use the ratio 
    column to assess proportional benefit."),
              HTML(" <code>diff_recruits_rest</code> tells you how many additional fish the 
    restoration scenario is projected to produce compared to the baseline. This is 
    the most direct measure of restoration benefit in terms of fish numbers. 
    <code>ratio_rest_med</code> tells you what fraction of baseline recruitment is 
    added by the restoration — a value of 1.20 means 20% more recruits than the 
    baseline. The ratio is more useful for comparing sites that differ in size or 
    for understanding benefit in dry years when absolute numbers are small but 
    proportional gains may still be meaningful.")
            ),
            tags$li(
              strong("Uncertainty bounds reflect model parameter uncertainty, not future variability."),
              " The lower and upper bounds in the downloaded CSV (e.g., ", 
              code("lower_ratio_rest"), ", ", code("upper_ratio_rest"), 
              ") reflect the range of plausible parameter values in the Bayesian population 
    model. They do not represent a prediction interval over future years. A wide 
    interval means the model parameters are uncertain for that reach or flow 
    condition; a narrow interval means the estimate is more precisely constrained."
            ),
            tags$li(
              strong("Dry-year results warrant caution."),
              " In very dry years (low x-axis values), both baseline recruitment and the 
    restoration benefit can be near zero because flows are insufficient to inundate 
    meaningful floodplain area at most sites. Very small absolute differences in 
    these years can produce large or unstable ratio values. When comparing sites, 
    consider summarizing across years or focusing on years with non-trivial baseline 
    recruitment."
            ),
            tags$li(
              strong("The restoration scenario is additive to the baseline."),
              " The tool models the restoration site as providing habitat above and beyond 
    what exists in the baseline river channel. The baseline already reflects the 
    existing floodplain inundation at the natural channel geometry. The restoration 
    benefit you see in the plots and output represents the incremental contribution 
    of the proposed restoration area."
            )
          ),
          br(),
          p("Disclaimer: This Beta Version is still in development and should be used for 
    demo and research purposes only. The tool should not be used for project planning 
    without prior USACE certification. Neither the authors nor the U.S. Army Corps of 
    Engineers accepts responsibility or liability for the model's use by third parties.",
            style = "color: red; font-size: 12px; border-top: 1px solid #ccc; padding-top: 10px;")
          #)
        ), # End of "Home" Tab Panel
        
        tabPanel("Hydrology",
                 p(HTML("<b>Hydrographs for specific years selected by the user</b>.
                         The salmon line shows the given year, the blue line indicates the change in discharge
                         with user-defined environmental flow parameters, and the dashed line
                         indicates the median daily flow from 1993 through 2020.")),
                 plotOutput("hydroPlot")),
        
        tabPanel(
          "Habitat",
          p(HTML(
            "<b>Inundated floodplain habitat under baseline, restoration, and 
    environmental flow scenarios.</b> Each row of panels corresponds to one 
    of the simulated flow years. The left panels show total inundated habitat 
    in acres over the March–July period: the black line is the baseline 
    (existing conditions), the salmon line is the total habitat with the 
    restoration scenario added, and the blue line is the total habitat with 
    the environmental flow added to the baseline. The right panels show 
    habitat above the baseline level only, isolating the incremental 
    contribution of each scenario."
          )),
          p(HTML(
            "<b>What to look for:</b> In a typical restoration scenario, the salmon 
    line should track above the black line at flows sufficient to inundate 
    the restoration area, and the two lines should converge at low flows 
    when neither the baseline nor the restoration area is inundated. If the 
    salmon line does not rise above the black line at any flow level, check 
    that the restoration acreage is correctly specified or that the uploaded 
    CSV contains non-zero habitat values at the expected discharge levels. 
    The right-hand panels are particularly useful for understanding at what 
    point in the season and at what discharge level the restoration begins 
    to contribute meaningfully above the baseline."
          )),
          plotOutput("habPlot")
        ),
        
        tabPanel(
          "Recruitment",
          p(HTML(
            "<b>Estimated age-0 Rio Grande Silvery Minnow recruitment under baseline, 
    restoration, and environmental flow scenarios.</b> Each point represents 
    one simulated flow year. The x-axis shows the flow index for that year 
    (seasonal maximum of the 20-day rolling minimum discharge), so points 
    toward the right represent wetter years with more sustained high flows."
          )),
          p(HTML(
            "In the <b>upper panel</b>, circles show baseline recruitment (log scale), 
    squares show recruitment under the restoration scenario, and triangles show 
    recruitment under the environmental flow scenario. Points that cluster 
    closely together indicate years where the scenarios produce similar 
    outcomes to the baseline; widely separated points indicate years where 
    the scenarios produce substantially more recruits. Note that the log scale 
    compresses differences at high recruitment levels — refer to the lower 
    panel and the downloaded data for absolute differences."
          )),
          p(HTML(
            "In the <b>lower panel</b>, each point shows the difference in recruits 
    above the baseline for the restoration (red squares) and environmental 
    flow (blue triangles) scenarios, expressed in thousands of fish. Years 
    labeled near the points can be used to identify which specific flow years 
    drive the largest restoration benefits. If all points in the lower panel 
    are near zero, this may indicate that the restoration area is too small 
    to produce a detectable signal at the flow levels present in the 
    selected years, or that the flow-habitat relationship specified for the 
    site does not generate meaningful inundation above the baseline."
          )),
          plotOutput("recruitPlot")
        ),
        
        tabPanel("Map",
                 h4(HTML("<b>Change the coordinates in the sidepanel to visualize flow permanence data for the location of a potential restoration project.</b>")),
                 leafletOutput("Map", height = 500),
                 p("Flow permanence data comes from the River Eyes Program. The
                   value for each 0.5 stream mile segment gives the percentage of
                   years from 2007 to 2023 where that segment did not dry completely
                   during July. For example, a value of 80% means that that segment
                   did not dry out during 80% of the years.")),
        
        tabPanel(
          "About This Model",
          
          h3("About the RGSM Recruitment Tool"),
          
          # ── OVERVIEW ──────────────────────────────────────────────────────────────
          h4("Overview"),
          p("This tool estimates age-0 Rio Grande Silvery Minnow (RGSM) recruitment — 
    the number of juvenile fish produced during spring spawning — as a function of floodplain inundation during the spring spawning 
    period. It is designed to support ", strong("relative comparisons"), " of 
    floodplain restoration scenarios, helping practitioners and agency planners 
    identify which restoration sites or management actions are most likely to 
    increase RGSM recruitment across a range of historical flow conditions. It 
    is not designed to generate absolute population forecasts or to substitute 
    for species-level population viability analysis."),
          p("The tool covers three management reaches of the Middle Rio Grande — 
    San Acacia, Isleta, and Angostura — and simulates recruitment under 
    baseline, floodplain restoration, and environmental flow scenarios using 
    the full historical streamflow record from 1993 to 2020."),
          
          # ── SCIENTIFIC FOUNDATION ─────────────────────────────────────────────────
          h4("Scientific Foundation"),
          p("The tool builds on a state-of-the-art hierarchical Bayesian population 
    model developed by Yackulic et al. (2022) and published in ",
            em("Ecosphere"), ":"),
          tags$blockquote(
            style = "border-left: 4px solid #2E74B5; padding-left: 12px; 
             color: #333; margin: 10px 24px;",
            p("Yackulic, C.B., et al. (2022). Quantifying flow and nonflow management 
      impacts on an endangered fish by integrating data, research, and expert 
      opinion. ", em("Ecosphere"), ", 13(8), e4240. ",
              a("https://doi.org/10.1002/ecs2.4240",
                href = "https://doi.org/10.1002/ecs2.4240", target = "_blank"))
          ),
          p("That model integrates three streams of evidence within a single 
    framework:"),
          tags$ul(
            tags$li(strong("Long-term RGSM monitoring data"), " collected across the 
             Middle Rio Grande, providing the empirical basis for the 
             relationship between flow conditions and recruitment outcomes."),
            tags$li(strong("USGS streamflow records"), " from gages at Albuquerque 
             and San Acacia, used to characterize the hydrologic conditions 
             experienced by spawning fish in each year."),
            tags$li(strong("Structured expert elicitation"), " of spawning phenology 
             parameters — including the flow cues that trigger spawning and 
             the duration of habitat required for larval survival — that 
             cannot be estimated reliably from monitoring data alone.")
          ),
          p("For the Recruitment Tool, this framework was extended by the same 
    modeling team, including the principal developer of the original population 
    model, to incorporate site-specific estimates of floodplain inundation 
    extent derived from two-dimensional (2D) hydraulic modeling conducted at 
    seven discharge levels between 700 and 7,000 cfs. These inundation 
    estimates replace the abstract habitat covariates in the original model 
    with physically grounded, spatially explicit measures of floodplain area, 
    enabling the model to generate recruitment projections for specific 
    restoration sites and acreages. This extension is described in a manuscript currently in review 
    at ", em("Ecosphere"), "."),
          
          # ── VALIDATION ────────────────────────────────────────────────────────────
          h4("Model Validation"),
          p("The original Yackulic et al. (2022) model was validated against 
    independent RGSM recruitment observations not used in model fitting, 
    demonstrating that the model captures the observed relationship between 
    spring hydrology and RGSM recruitment across the full range of wet and 
    dry years in the historical record."),
          p("The extended model implemented in this tool was similarly validated 
    against out-of-sample RGSM monitoring data prior to deployment. This 
    validation confirmed that the incorporation of 2D hydraulic inundation 
    estimates did not degrade the model's predictive skill and that the 
    tool produces reliable recruitment projections across the range of 
    historical flow conditions. Full details of the validation analyses, 
    including the withheld data, comparison metrics, and diagnostic plots, 
    are provided in the in-review manuscript and the associated Zenodo 
    archive (see Data and Code Availability below)."),
          p("The tool is currently undergoing the model certification process by the U.S. Army Corps of 
    Engineers for use in planning studies that compare the recruitment 
    benefits of floodplain restoration sites along the Middle Rio Grande. 
    The USACE certification document provides a full summary of the 
    technical basis for certification and is available upon request."),
          
          # ── WHAT THE MODEL PREDICTS ───────────────────────────────────────────────
          h4("What the Model Predicts — and What It Does Not"),
          p("Understanding the scope of the model helps ensure it is applied 
    appropriately and that outputs are interpreted correctly."),
          p(strong("The model does predict:")),
          tags$ul(
            tags$li("Age-0 RGSM recruitment as a function of floodplain inundation 
             extent during the March–June simulation window."),
            tags$li("The incremental recruitment benefit of a restoration scenario 
             relative to baseline conditions, expressed as additional recruits 
             above baseline or as a ratio of restoration to baseline 
             recruitment."),
            tags$li("How the recruitment benefit of a restoration scenario varies 
             across wet and dry years, based on the 1993–2020 historical 
             streamflow record."),
            tags$li("Uncertainty in recruitment projections arising from parameter 
             uncertainty in the Bayesian population model, reported as 
             10th and 90th percentile bounds on all estimates.")
          ),
          p(strong("The model does not predict:")),
          tags$ul(
            tags$li("Absolute RGSM population size or population-level recovery 
             outcomes."),
            tags$li("Survival of age-0 fish beyond July 1, or recruitment to older age classes."),
            tags$li("The effects of predation, water temperature, water quality, 
             or habitat features other than inundation extent."),
            tags$li("Recruitment under future climate or flow conditions not 
             represented in the 1993–2020 historical record."),
            tags$li("Outcomes at sites outside the Middle Rio Grande reaches 
             for which the model was developed and validated.")
          ),
          
          # ── GENERIC VS CUSTOM ─────────────────────────────────────────────────────
          h4("Choosing Between the Generic and Custom Input Options"),
          p(HTML(
            "The tool offers two ways to specify the flow-habitat relationship 
    for a restoration site, suited to different stages of project 
    development:"
          )),
          p(HTML("<b>Generic (scale by acres)</b> — appropriate for early-stage 
    screening and site prioritization when site-specific hydraulic data 
    are not yet available. The tool applies a representative per-acre 
    flow-inundation curve derived from existing 2D hydraulic modeling 
    in the San Acacia reach and scales it linearly by the acreage you 
    enter. This approximation works reasonably well for sites with 
    inundation dynamics broadly similar to the reference restoration site curves, but 
    may over- or underestimate habitat at sites with unusual topography, 
    connectivity, or channel geometry. It is well suited to answering 
    the question 'which sites are worth investigating further?' rather 
    than 'exactly how much recruitment will this site produce?'")),
          p(HTML("<b>Custom flow-habitat curve</b> — appropriate for later-stage 
    analyses, formal planning documents, and any application where 
    accuracy matters more than convenience. Upload a two-column CSV 
    with discharge (<code>cfs</code>) and inundated habitat in acres 
    (<code>hab</code>) derived from site-specific hydraulic modeling. 
    The tool interpolates linearly between your supplied values and 
    clamps to the nearest endpoint outside the uploaded flow range. 
    This option is recommended whenever 1D or 2D hydraulic modeling 
    has been conducted for a site, and may be required for results intended 
    to appear in formal planning or regulatory documents.")),
          p("If you are uncertain which option best fits your site, you can 
    run the tool with both for the same site and compare the results. 
    A large difference between the two suggests that the per-acre 
    approximation may not be a good fit for your site's specific 
    inundation dynamics."),
          
          # ── UNCERTAINTY ───────────────────────────────────────────────────────────
          h4("Understanding Uncertainty Bounds"),
          p("All recruitment projections in the downloaded CSV include lower and 
    upper bounds — for example, ", code("lower_ratio_rest"), " and ", 
            code("upper_ratio_rest"), ". These are the 10th and 90th percentiles 
    of the posterior distribution of recruitment estimates from the 
    Bayesian population model and reflect uncertainty in the model 
    parameters, particularly:"),
          tags$ul(
            tags$li("The reach-specific productivity parameters estimated from 
             RGSM monitoring data."),
            tags$li("The spawning phenology parameters informed by expert 
             elicitation, which carry uncertainty arising from the 
             elicitation process itself.")
          ),
          p("These bounds are ", strong("not"), " prediction intervals over 
    future years and should not be interpreted as the range of 
    recruitment you would expect to observe in any given year. They 
    are a measure of how precisely the model parameters are estimated 
    given the available data."),
          p("For site comparison, it is good practice to examine whether the 
    uncertainty intervals of two sites overlap substantially. 
    Overlapping intervals indicate that the available data do not 
    clearly distinguish the recruitment benefit of one site from 
    another, and the comparison should be treated with appropriate 
    caution. Non-overlapping intervals provide stronger evidence that 
    one site is genuinely more beneficial than the other under the 
    modeled conditions."),
          
          # ── DATA AND CODE ─────────────────────────────────────────────────────────
          h4("Data and Code Availability"),
          p("All data and code underlying the RGSM Recruitment Tool are 
    publicly archived and freely accessible:"),
          tags$ul(
            tags$li(
              strong("Yackulic et al. (2022)"), " — population model data and code: ",
              a("https://doi.org/10.5281/zenodo.6842195",
                href   = "https://doi.org/10.5281/zenodo.6842195",
                target = "_blank")
            ),
            tags$li(
              strong("Recruitment Tool extension and validation"),
              " (manuscript in review at ", em("Ecosphere"), "): ",
              a("https://doi.org/10.5281/zenodo.17641889",
                href   = "https://doi.org/10.5281/zenodo.17641889",
                target = "_blank")
            )
          ),
          p("The Zenodo archives include the RGSM monitoring data, USGS 
    streamflow records, 2D hydraulic modeling outputs, all R code 
    used to fit and validate the model, and annotated scripts for 
    reproducing the figures and tables in the supporting manuscripts.")
        )
      )
    )
  )
)

# Server function ----------------------------------------------------------
server <- function(input, output, session) {
  
  predefined_values <- c(1993, 1995:2020)
  
  final_input <- reactive({
    if (input$usePreset) {
      predefined_values
    } else {
      as.numeric(unlist(strsplit(input$integerInput, ",")))
    }
  })
  
  observeEvent(input$integerInput, { if (!input$usePreset) final_input() })
  observeEvent(input$usePreset,    { final_input() })
  
  # ---------------------------------------------------------------------------
  # NEW REACTIVE: build the restoration hab lookup from whichever input mode
  # is selected.  Returns a tibble with columns: cfs, reach, hab
  # (hab = incremental restoration habitat in acres above baseline)
  # ---------------------------------------------------------------------------
  hab_lookup_rest_reactive <- reactive({
    
    restReach   <- input$reachInput
    otherReaches <- reaches[reaches != restReach]
    
    # Zero-habitat placeholder for non-restoration reaches
    nonRestLookup <- tibble(
      cfs   = rep(seq(1, 7000, 1), length(otherReaches)),
      reach = rep(otherReaches, each = 7000),
      hab   = 0
    )
    
    if (input$restInputMode == "custom") {
      
      # Require a file to have been uploaded
      validate(
        need(!is.null(input$restCurveFile),
             "Please upload a CSV file with columns 'cfs' and 'hab'.")
      )
      
      # Read and validate the uploaded CSV
      custom_curve <- tryCatch(
        read_csv(input$restCurveFile$datapath, show_col_types = FALSE),
        error = function(e) NULL
      )
      
      validate(
        need(!is.null(custom_curve),
             "Could not read the uploaded file. Please check it is a valid CSV."),
        need(all(c("cfs", "hab") %in% colnames(custom_curve)),
             "Uploaded CSV must contain columns named 'cfs' and 'hab'.")
      )
      
      # Interpolate to integer cfs values 1:7000
      # rule = 2 clamps values outside the uploaded range to nearest endpoint
      cfs_full <- tibble(cfs = 1:7000)
      custom_interp <- cfs_full %>%
        mutate(
          hab   = approx(custom_curve$cfs, custom_curve$hab,
                         xout = cfs, rule = 2)$y,
          reach = restReach
        )
      
      bind_rows(custom_interp, nonRestLookup) %>%
        left_join(reach_df) # Add reach numbers
      
    } else {
      # Generic mode: scale the pre-computed per-acre curve by restAcres
      lookup_per_acre_rest_reach %>%
        mutate(reach = restReach,
               hab   = hab * input$restAcres) %>%
        bind_rows(nonRestLookup) %>%
        left_join(reach_df) # Add reach numbers
    }
  })
  
  # ---------------------------------------------------------------------------
  # NEW REACTIVE: a human-readable "acres" label for plot annotations.
  # In generic mode this is simply restAcres.
  # In custom mode we report the maximum habitat value in the uploaded curve
  # as a proxy for project size.
  # ---------------------------------------------------------------------------
  rest_acres_label <- reactive({
    if (input$restInputMode == "custom" && !is.null(input$restCurveFile)) {
      hab_lookup_rest_reactive() %>%
        filter(reach == input$reachInput) %>%
        summarize(max_hab = max(hab, na.rm = TRUE)) %>%
        pull(max_hab) %>%
        round(0)
    } else {
      input$restAcres
    }
  })
  # ---------------------------------------------------------------------------
  
  #Assesses if something is wrong with predictors
  preds_reactive <- reactive({
    
    input_values <- final_input()
    simQ         <- q_sim(input_values)
    augQ         <- aug_q_fun(
      input$startDate,
      input$endDate,
      input$augInput,
      simQ
    )
    
    hab_lookup_rest <- hab_lookup_rest_reactive() %>%
      bind_rows(hab_lookup_baseline) %>%
      group_by(reach, cfs, reach_num) %>%
      summarize(hab = sum(hab), .groups = "drop")
    
    lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline)
    lcc_rest     <- preds_calc(expert_num, simQ, hab_lookup_rest)
    lcc_hydro    <- preds_calc(expert_num, augQ, hab_lookup_baseline)
    
    maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
    maxmins_rest     <- maxmins_fun(hab_lookup_rest, simQ, reaches)
    maxmins_hydro    <- maxmins_fun(hab_lookup_baseline, augQ, reaches)
    
    preds_array <- array(
      NA,
      dim = c(length(input_values), length(reaches), 3)
    )
    
    preds_array[, , 1] <- lcc_baseline - mean(cov_list)
    preds_array[, , 2] <- lcc_rest     - mean(cov_list)
    preds_array[, , 3] <- lcc_hydro    - mean(cov_list)
    
    preds_array
  })
  
  
  recruitment_results <- reactive({
    
    input_values <- final_input()
    simQ         <- q_sim(input_values)
    augQ         <- aug_q_fun(
      input$startDate,
      input$endDate,
      input$augInput,
      simQ
    )
    
    hab_lookup_rest <- hab_lookup_rest_reactive() %>%
      bind_rows(hab_lookup_baseline) %>%
      group_by(reach, cfs, reach_num) %>%
      summarize(hab = sum(hab), .groups = "drop")
    
    lcc_baseline <- preds_calc(expert_num, simQ, hab_lookup_baseline)
    lcc_rest     <- preds_calc(expert_num, simQ, hab_lookup_rest)
    lcc_hydro    <- preds_calc(expert_num, augQ, hab_lookup_baseline)
    
    maxmins_baseline <- maxmins_fun(hab_lookup_baseline, simQ, reaches)
    maxmins_rest     <- maxmins_fun(hab_lookup_rest, simQ, reaches)
    maxmins_hydro    <- maxmins_fun(hab_lookup_baseline, augQ, reaches)
    
    preds_array <- array(NA,
                         dim = c(length(input_values), length(reaches), 3)
    )
    
    preds_array[, , 1] <- lcc_baseline - mean(cov_list)
    preds_array[, , 2] <- lcc_rest     - mean(cov_list)
    preds_array[, , 3] <- lcc_hydro    - mean(cov_list)
    
    age0_dat <- forecast_recruit(
      sim_mod,
      preds_array,
      input_values,
      maxminsBase = maxmins_baseline,
      maxminsAug  = maxmins_hydro,
      maxminsRest = maxmins_rest,
      reaches
    )
    
    age0_dat
  })
  
  output$hydroPlot <- renderPlot({
    input_values <- final_input()
    simQ         <- q_sim(input_values)
    augQ         <- aug_q_fun(input$startDate, input$endDate, input$augInput, simQ)
    flow_plot(input_values, augQ, input$reachInput)
  },
  height = reactive({ length(final_input()) * 200 }))
  
  output$habPlot <- renderPlot({
    input_values <- final_input()
    simQ         <- q_sim(input_values)
    augQ         <- aug_q_fun(input$startDate, input$endDate, input$augInput, simQ)
    
    # Use the unified reactive instead of inline construction
    hab_lookup_rest <- hab_lookup_rest_reactive()
    
    flow_rest_plot(input_values, input$reachInput, augQ,
                   hab_lookup_baseline, hab_lookup_rest)
  },
  height = reactive({ length(final_input()) * 200 }))
  
  output$recruitPlot <- renderPlot({
    
    recruit_plot(
      recruitment_results(),
      rest_reach = input$reachInput,
      rest_acres = rest_acres_label(),
      aug_vol    = input$augInput
    )
    
  }, height = 800)
  
  output$Map <- renderLeaflet({
    pal <- colorNumeric(palette = "RdYlBu", domain = lines_sf$flow_perm)
    
    leaflet(lines_sf) %>%
      addProviderTiles(providers$Esri.WorldImagery,
                       options = providerTileOptions(opacity = 0.5)) %>%
      addTiles(
        urlTemplate = "https://server.arcgisonline.com/ArcGIS/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile/{z}/{y}/{x}",
        options = tileOptions(tileSize = 256)
      ) %>%
      addPolylines(color   = ~ pal(flow_perm),
                   weight  = 5,
                   opacity = 0.8,
                   popup   = ~ htmltools::htmlEscape(flow_perm)) %>%
      addLegend(pal = pal, values = ~ flow_perm,
                title = "Flow\nPermanence", opacity = 1) %>%
      addCircleMarkers(lng = input$lng, lat = input$lat,
                       color = "red", radius = 6, fillOpacity = 1,
                       layerId = "user_point")
  })
  
  observe({
    leafletProxy("Map") %>%
      clearGroup("user_point") %>%
      addCircleMarkers(lng = input$lng, lat = input$lat,
                       color = "red", radius = 6, fillOpacity = 1,
                       group = "user_point")
  })
  
  datasetInput <- reactive({
    
    recruitment_results() %>%
      filter(reach %in% input$reachInput) %>%
      mutate(
        input_hydros = paste(final_input(), collapse = ","),
        rest_acres   = rest_acres_label(),
        rest_mode    = input$restInputMode,
        aug_vol      = input$augInput,
        aug_start    = input$startDate,
        aug_end      = input$endDate
      ) %>%
      select(!matches("aug|hydro|prop")) %>%
      mutate(across(c(cfs_baseline, acres_inund_baseline,
                      acres_inund_rest, recruits_baseline, lower_age0, upper_age0,
                      recruits_rest, lower_age0_rest, upper_age0_rest,
                      diff_recruits_rest, lower_diff_rest, upper_diff_rest
      ), round, digits = 0),
      across(c(ratio_rest_med, lower_ratio_rest, upper_ratio_rest), round, digits = 2))
    
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste0("rgsm_output_", Sys.Date(), ".csv") },
    content  = function(file) { write.csv(datasetInput(), file, row.names = FALSE) }
  )
  
  observe({
    print(
      preds_reactive()[,,3] - preds_reactive()[,,1]
    )
  })
  
}

shinyApp(ui = ui, server = server)
