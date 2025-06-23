## model fitting for USFWS bull trout SSA

#### setup ####

## load libraries
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(MARSS)

## set directories
output_dir <- here("output")

## first year of data to consider in model
yr_first <- 1991

## read data
adult_data <- read_csv(file = here("data", "clean",
                                   "bull_trout_SSA_data_all_states_adults.csv"))


#### data formatting ####

## trim years & reshape to "wide" format for MARSS
adult_smry <- adult_data %>%
  filter(metric == "abundance") %>%
  select(-metric) %>%
  filter(year >= yr_first) %>%
  pivot_wider(names_from = year,
              values_from = value,
              names_prefix = "yr") %>%
  arrange(state, recovery_unit, core_area) %>%
  rowwise(state:source) %>%
  mutate(n_yrs = sum(!is.na(c_across(everything())))) %>%
  ungroup() %>%
  filter(n_yrs >= 10)

## write data summary to file
adult_smry %>% 
  select(state, core_area, popn_stream, source, n_yrs) %>%
  mutate(source = stringr::str_replace(source, "_", " ")) %>%
  write_csv(file = here(output_dir, "bull_trout_SSA_data_summary_adults.csv"))

yy <- adult_smry %>%
  select(-n_yrs)

## number of core areas (processes, x)
nc <- length(unique(yy$core_area))
## number of locations (streams/rivers, y)
nr <- nrow(yy)

## number of sites per core area
core_tbl <- yy %>% 
  group_by(state, recovery_unit, core_area) %>%
  summarise(n = length(core_area))


#### MARSS setup ####

## observation eqn

## empty Z matrix for mapping obs to processes
ZZ <- matrix(0, nr, nc)

## loop over core areas to set cols of Z
for (jj in 1:nc) {
  ## seq for row indices
  ll <- seq(as.integer(core_tbl[jj, 4]))
  ## last row
  xx <- ifelse(jj == 1, 0, max(ii))
  ## seq for row indices
  ii <- ll + xx
  ## assign 1's to respective rows/obs
  ZZ[ii, jj] <- 1
}

## offsets for obs (a); data de-meaned so all 0's
AA <- matrix(0, nr, 1)

## covariance matrix for obs (R)
RR <- matrix(list(0), nr, nr)
diag(RR) <- yy$source

## process eqn

## interactions matrix (B); set to I for RW's
BB <- diag(nc)

## bias terms (u); each core area gets a unique bias term
UU <- core_tbl %>%
  select(-n) %>%
  unite("core_area", state:core_area, sep = ": ") %>%
  as.matrix(nrow = nc, ncol = 1)

## cov matrix for processes (Q)
QQ <- matrix(list(0), nc, nc)
## diagonal and unequal
# diag(QQ) <- core_tbl %>%
#   select(-n) %>%
#   unite("core_area", state:core_area, sep = ": ") %>%
#   unlist()
## diagonal and equal (IID)
diag(QQ) <- rep("q", nc)

## data for fitting

## reformat for MARSS
yy <- yy %>%
  ## drop ID cols
  select(-(state:source)) %>%
  ## convert to matrix
  as.matrix() %>%
  + 1 %>%
  ## log-transform
  log() %>%
  ## remove the mean
  zscore(mean.only = FALSE)

## remove upstream tibbles that are no longer necessary
rm(adult_data,
   core_tbl)


#### model fitting ####

## model list
mod_list <- list(
  B = BB,
  U = UU,
  Q = QQ,
  Z = ZZ,
  A = AA,
  R = RR
)

## control list
con_list <- list(
  maxit = 5000
)

## fit base model
mod_fit <- MARSS(yy, model = mod_list, control = con_list)

## save fitted model object
saveRDS(mod_fit, here(output_dir, "adult_model_fits.rds"))

#### bootstrapped CI's ####

## bootstrap parameters from the Hessian matrix
mod_fit_CI90 <- MARSSparamCIs(mod_fit, method = "hessian", alpha = 0.1, nboot = 1000)

## save bootstrapped model object
saveRDS(mod_fit_CI90, here(output_dir, "adult_model_fits_CI90.rds"))

## extract bias params
bias_mean <- mod_fit_CI90$parMean[grep("U.", names(mod_fit_CI90$parMean))]

## summary table of bias CI's
bias_smry <- cbind(mod_fit_CI90$par.lowCI$U,
                   bias_mean,
                   mod_fit_CI90$par.upCI$U) %>%
  round(4) %>%
  as.data.frame()
## better row names
rownames(bias_smry) <- gsub("(U.)(.*)", "\\2", names(bias_mean))
## better col names
colnames(bias_smry) <- c("CI_lo_90", "mean", "CI_up_90")

## summarize trends
## negative
neg <- bias_smry %>%
  apply(1, function(x) x < 0) %>%
  t() %>%
  apply(1, all)
## positive
pos <- bias_smry %>%
  apply(1, function(x) x > 0) %>%
  t() %>%
  apply(1, all)
## add trend col
bias_smry$trend = "0"
bias_smry$trend[neg] <- "-"
bias_smry$trend[pos] <- "+"

## write bias summary to file
bias_smry %>% 
  write.csv(file = here(output_dir, "bull_trout_SSA_all_states_adults_biases.csv"))


#### late-period trend ####

## covariate for trend
cc <- ((seq(yr_first, 2020) >= 2008) * 1) %>%
  matrix(nrow = 1)

## update model list
mod_list$U <- matrix(0, ncol = 1, nrow = nc)

mod_list$C <- UU

mod_list$c <- cc

## fit late-period model
mod_fit_late <- MARSS(yy, model = mod_list, control = con_list)

## save fitted model object
saveRDS(mod_fit_late, here(output_dir, "adult_model_fits_late.rds"))

#### late-period bootstrapped CI's ####

## bootstrap parameters from the Hessian matrix
mod_fit_late_CI90 <- MARSSparamCIs(mod_fit_late, method = "hessian", alpha = 0.1, nboot = 1000)

## save bootstrapped model object
saveRDS(mod_fit_late_CI90, here(output_dir, "adult_model_fits_late_CI90.rds"))

## extract bias params
bias_late_mean <- mod_fit_late_CI90$parMean[grep("U.", names(mod_fit_late_CI90$parMean))]

## summary table of bias CI's
bias_smry <- cbind(mod_fit_late_CI90$par.lowCI$U,
                   bias_late_mean,
                   mod_fit_late_CI90$par.upCI$U) %>%
  round(4) %>%
  as.data.frame()
## better col names
colnames(bias_smry) <- c("CI_lo_90", "mean", "CI_up_90")
## better row names
rownames(bias_smry) <- gsub("(U.)(.*)", "\\2", names(bias_late_mean))

## summarize trends
## negative
neg <- bias_smry %>%
  apply(1, function(x) x < 0) %>%
  t() %>%
  apply(1, all)
## positive
pos <- bias_smry %>%
  apply(1, function(x) x > 0) %>%
  t() %>%
  apply(1, all)
## add trend col
bias_smry$trend = "0"
bias_smry$trend[neg] <- "-"
bias_smry$trend[pos] <- "+"

## write bias summary to file
bias_smry %>% 
  write.csv(file = here(output_dir, "bull_trout_SSA_all_states_adults_late_biases.csv"))


