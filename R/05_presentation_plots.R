## summary info for USFWS bull trout SSA

#### setup ####

## required libraries
## library(here)
## library(readr)
## library(MARSS)
library(tidyr)
library(dplyr)
library(viridisLite)

## set directories
clean_data_dir <- here::here("data", "clean")
output_dir <- here::here("output")

## first & last years of data considered in model
yr_first <- 1991
yr_last <- 2020

## time index for plots
t_index <- seq(yr_first, yr_last)


#### get observed adult data ####

## read data
adult_data <- readr::read_csv(file = here(clean_data_dir,
                                          "bull_trout_SSA_data_all_states_adults.csv"))

## trim years, reshape to "wide" format & transform for MARSS
yy <- adult_data %>%
  ## extract abundance data only
  filter(metric == "abundance") %>%
  ## drop metric
  select(-metric) %>%
  ## filter appropriate years
  filter(year >= yr_first) %>%
  ## pivot to wide
  pivot_wider(names_from = year,
              values_from = value,
              names_prefix = "yr") %>%
  ## sort by state & core area
  arrange(state, recovery_unit, core_area) %>%
  ## get total years of data
  rowwise(state:source) %>%
  mutate(n_yrs = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  ## select only those local popns with >=10 yrs of data
  filter(n_yrs >= 10) %>%
  select(-n_yrs) 


#### entire time period ####

## get model fits
mod_fit_CI90 <- readRDS(file = here(output_dir, "adult_model_fits_CI90.rds"))

## extract core area names
core_areas <- MARSS:::coef.marssMLE(mod_fit_CI90, matrix)$U %>%
  rownames() %>%
  strsplit(": ") %>%
  do.call(what = rbind) %>%
  as.data.frame() %>%
  `colnames<-`(c("state", "recovery_unit", "core_area"))

## number of core areas
n_cores <- nrow(core_areas)

## write plots to pdf - paper
# pdf(file = here(output_dir, "bull_trout_SSA_adult_summary_plots.pdf"),
#     height = 6, width = 9)

bias_smry <- read.csv(file = here(output_dir, "bull_trout_SSA_all_states_adults_biases.csv"))

bias_smry$mean %>% hist(breaks = seq(-0.15, 0.15, 0.03), plot = FALSE)


## write plots to pdf - presentation
pdf(file = here(output_dir, "bull_trout_SSA_adult_summary_plots_talk.pdf"),
    height = 6, width = 9)

## loop over core areas by state
for(i in 1:n_cores) {
  
  ## extract obs data for core area
  tmp <- yy %>%
    filter(core_area == core_areas[i,"core_area"]) %>%
    ## drop ID cols
    select(-(state:source)) %>%
    ## convert to matrix
    as.matrix() %>%
    + 1 %>%
    ## log-transform
    log() %>%
    ## remove the mean
    MARSS:::zscore(mean.only = FALSE) %>%
    ## pivot to long
    t()
  
  ## number of local popns in the core area
  nn <- ncol(tmp)
  
  ## trim states to data
  tmp_na <- (!is.na(tmp)) %>%
    apply(1, any)
  tmp_sf <- seq(length(t_index))[tmp_na] %>%
    range()
  tmp_t <- seq(tmp_sf[1], tmp_sf[2])
  
  ## get estimated trend line
  tmp2 <- mod_fit_CI90$states[i, tmp_t]
  
  ## get SE of estimated trend line
  tmp3u <- tmp2 + 1.645*mod_fit_CI90$states.se[i, tmp_t]
  tmp3l <- tmp2 - 1.645*mod_fit_CI90$states.se[i, tmp_t]
  
  ## set colormap
  clr <- mako(nn, begin = 0.4, end = 0.9)
  
  ## plot abundance index
  par(mai = c(0.9, 1.4, 0.6, 0.1))
  matplot(t_index, tmp,
          type = "o", lty = "solid", pch = 16, col = clr,
          ylim = range(c(tmp, tmp2, tmp3u, tmp3l), na.rm = TRUE),
          # for talk
          cex.axis = 1.5, cex.lab = 1.5, cex = 1.5, lwd = 2,
          las = 1, xlab = "Year", ylab = "")
  mtext("Abundance index", side = 2, cex = 1.5, line = 4.5)
  mtext(paste0(core_areas[i,"state"], ": ", core_areas[i,"core_area"]),
        # for talk
        cex = 1.5,
        side = 3, line = 0.5, adj = 0)
  ## plot estimated trend
  lines(t_index[tmp_t], tmp2, lwd = 3)
  lines(t_index[tmp_t], tmp3u, lwd = 2, col = "gray")
  lines(t_index[tmp_t], tmp3l, lwd = 2, col = "gray")
  
}

dev.off()


#### late-period trends ####

## get model fits
mod_fit_late_CI90 <- readRDS(file = here(output_dir, "adult_model_fits_late_CI90.rds"))

## write plots to pdf
pdf(file = here(output_dir, "bull_trout_SSA_adult_summary_plots_late_period.pdf"),
    height = 6, width = 9)

## loop over core areas by state
for(i in 1:n_cores) {
  
  ## extract obs data for core area
  tmp <- yy %>%
    filter(core_area == core_areas[i,"core_area"]) %>%
    ## drop ID cols
    select(-(state:source)) %>%
    ## convert to matrix
    as.matrix() %>%
    + 1 %>%
    ## log-transform
    log() %>%
    ## remove the mean
    MARSS:::zscore(mean.only = FALSE) %>%
    ## pivot to long
    t()
  
  ## number of local popns in the core area
  nn <- ncol(tmp)
  
  ## trim states to data
  tmp_na <- (!is.na(tmp)) %>%
    apply(1, any)
  tmp_sf <- seq(length(t_index))[tmp_na] %>%
    range()
  tmp_t <- seq(tmp_sf[1], tmp_sf[2])
  
  ## get estimated trend line
  tmp2 <- mod_fit_CI90$states[i, tmp_t]
  
  ## get SE of estimated trend line
  tmp3u <- tmp2 + 1.645*mod_fit_CI90$states.se[i, tmp_t]
  tmp3l <- tmp2 - 1.645*mod_fit_CI90$states.se[i, tmp_t]
  
  ## set colormap
  clr <- mako(nn, begin = 0.4, end = 0.9)
  
  ## plot abundance index
  par(mai = c(0.9, 0.9, 0.6, 0.1))
  matplot(t_index, tmp,
          type = "o", lty = "solid", pch = 16, col = clr,
          ylim = range(c(tmp, tmp2, tmp3u, tmp3l), na.rm = TRUE),
          las = 1, xlab = "Year", ylab = "Abundance index")
  mtext(paste0(core_areas[i,"state"], ": ", core_areas[i,"core_area"]),
        side = 3, line = 0.5, adj = 0)
  ## plot estimated trend
  lines(t_index[tmp_t], tmp2, lwd = 3)
  lines(t_index[tmp_t], tmp3u, lwd = 2, col = "gray")
  lines(t_index[tmp_t], tmp3l, lwd = 2, col = "gray")
  
}

dev.off()



#### get observed juvenile data ####

## read data
juvie_data <- readr::read_csv(file = here(clean_data_dir,
                                          "bull_trout_SSA_data_all_states_juveniles.csv"))

## trim years, reshape to "wide" format & transform for MARSS
yy <- juvie_data %>%
  arrange(year) %>%
  ## extract abundance data only
  # filter(metric == "abundance") %>%
  ## drop metric
  select(-metric) %>%
  ## filter appropriate years
  filter(year >= yr_first) %>%
  ## pivot to wide
  pivot_wider(names_from = year,
              values_from = value,
              names_prefix = "yr") %>%
  ## sort by state & core area
  arrange(state, recovery_unit, core_area) %>%
  ## get total years of data
  rowwise(state:source) %>%
  mutate(n_yrs = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  ## select only those local popns with >=10 yrs of data
  filter(n_yrs >= 10) %>%
  select(-n_yrs) 


#### entire time period ####

## get model fits
mod_fit_CI90 <- readRDS((file = here(output_dir, "juv_model_fits_CI90.rds")))

## extract core area names
core_areas <- MARSS:::coef.marssMLE(mod_fit_CI90, matrix)$U %>%
  rownames() %>%
  strsplit(": ") %>%
  do.call(what = rbind) %>%
  as.data.frame() %>%
  `colnames<-`(c("state", "recovery_unit", "core_area"))

## number of core areas
n_cores <- nrow(core_areas)

## write plots to pdf
pdf(file = here(output_dir, "bull_trout_SSA_juvenile_summary_plots.pdf"),
    height = 6, width = 9)

## loop over core areas by state
for(i in 1:n_cores) {
  
  ## extract obs data for core area
  tmp <- yy %>%
    filter(core_area == core_areas[i,"core_area"]) %>%
    ## drop ID cols
    select(-(state:source)) %>%
    ## convert to matrix
    as.matrix() %>%
    + 1 %>%
    ## log-transform
    log() %>%
    ## remove the mean
    MARSS:::zscore(mean.only = FALSE) %>%
    ## pivot to long
    t()
  
  ## number of local popns in the core area
  nn <- ncol(tmp)
  
  ## trim states to data
  tmp_na <- (!is.na(tmp)) %>%
    apply(1, any)
  tmp_sf <- seq(length(t_index))[tmp_na] %>%
    range()
  tmp_t <- seq(tmp_sf[1], tmp_sf[2])
  
  ## get estimated trend line
  tmp2 <- mod_fit_CI90$states[i, tmp_t]
  
  ## get SE of estimated trend line
  tmp3u <- tmp2 + 1.645*mod_fit_CI90$states.se[i, tmp_t]
  tmp3l <- tmp2 - 1.645*mod_fit_CI90$states.se[i, tmp_t]
  
  ## set colormap
  clr <- mako(nn, begin = 0.4, end = 0.9)
  
  ## plot abundance index
  par(mai = c(0.9, 0.9, 0.6, 0.1))
  matplot(t_index, tmp,
          type = "o", lty = "solid", pch = 16, col = clr,
          ylim = range(c(tmp, tmp2, tmp3u, tmp3l), na.rm = TRUE),
          las = 1, xlab = "Year", ylab = "Abundance index")
  mtext(paste0(core_areas[i,"state"], ": ", core_areas[i,"core_area"]),
        side = 3, line = 0.5, adj = 0)
  ## plot estimated trend
  lines(t_index[tmp_t], tmp2, lwd = 3)
  lines(t_index[tmp_t], tmp3u, lwd = 2, col = "gray")
  lines(t_index[tmp_t], tmp3l, lwd = 2, col = "gray")
  
}

dev.off()


