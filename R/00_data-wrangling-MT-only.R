## data munging for USFWS bull trout SSA

#### NOTES ####

## The data from MT cam in a format different than the tidy template I sent
## to state coordinators, so they need to be cleaned before being processed
## with the other states.

## There are some questionable data based on the metadata in the original
## file (~50 site/year combinations out of 8000+). We can exlcude these with
## the following `qaqc` set to TRUE

## flag to do QA/QC
qaqc <- TRUE

## metadata codes
mt_metadata_codes <- c("a", "b", "d", "e", "f", "g", "z")

## save file for Excel out
mt_file_name_out <- "USFWS_bull_trout_SSA_data_MT.xlsx"

#### setup ####

## load libraries
library(here)
library(readxl)
library(writexl)
library(readr)
library(dplyr)
library(tidyr)

## set directories
raw_data_dir <- here("data", "raw")
mt_data_dir <- here("data", "raw", "MT")

## expected column names
colnames_default <- c("dataset", "recovery unit", "core area", "popn/stream",
                      "metric", "method", "year", "value")

## better column names
colnames_nice <- c("dataset", "recovery_unit", "core_area", "popn_stream",
                   "metric", "source", "year", "value")


#### get metadata ####

## get all file names
all_files <- dir(raw_data_dir)

metadata_file <- grep("metadata", all_files, value = TRUE)

metadata <- read_csv(file.path(raw_data_dir, metadata_file))


#### read raw data ####

## get Excel file name
mt_file_name <- "FWS_MT_LP_All_Redd_Data_to_Scheuerell.xlsx"

## read Excel file
mt_data_all <- file.path(mt_data_dir, mt_file_name) %>%
  read_xlsx(col_type = "text", na = c("", "NA", "n/a", "na", ".", "-")) %>%
  ## add integer col for dataset ID to match other states
  ## specific number doesn't matter b/c will be linked to data dictionary
  mutate(dataset = seq(nrow(.)), .before = 1)

## clean raw data
mt_data_clean <- mt_data_all %>%
  ## drop rows with <10 years of data
  filter(TotalYearsCounted >= 10) %>%
  ## drop cols that don't match template
  select(!c(`Local Population 2015`, LifeHistory, Kovach_ID, MT_ID, TotalYearsCounted)) %>%
  ## years are cols spread out wide; pivot to tidy long format
  pivot_longer(cols = `1979`:`2020`, names_to = "year", values_to = "redds") %>%
  ## need to extract embedded metadata codes from redd counts
  ## create cols for `n_redds` and `notes`
  mutate(n_redds = gsub(pattern = "\\D{1,2}", replacement  = "", x = redds),
         notes = gsub(pattern = "\\d{1,3}", replacement  = "", x = redds)) %>%
  ## convert empty strings from metadata into NA
  mutate(notes = ifelse(notes == "", NA, notes)) %>%
  ## convert year and n_redds to integer
  mutate(year = as.integer(year),
         n_redds = as.integer(n_redds)) %>%
  ## drop col with orig counts/codes
  select(-redds) %>%
  ## create abundance and method cols to match template
  mutate(metric = "abundance",
         method = "redd count", .before = year) %>%
  ## rename cols to match template
  `colnames<-`(c(colnames_default, "notes"))

## possibly exclude questionable data by setting values to NA
if(qaqc) {
  mt_data_clean <- mt_data_clean %>%
    mutate(value = ifelse(notes %in% mt_metadata_codes, NA, value))
}
  
## write revamped file to Excel
mt_data_clean %>%
  write_xlsx(path = file.path(raw_data_dir, mt_file_name_out))


#### metadata ####

## add metadata to existing .csv file
mt_data_clean %>%
  ## select the datasets and notes  
  select(dataset, notes) %>%
  ## contingency table
  table() %>%
  as_tibble() %>%
  ## pivot wider
  pivot_wider(names_from = notes,
              values_from = n,
              names_prefix = "note_") %>%
  ## add state ID
  mutate(dataset = as.integer(dataset),
         state = "MT", .before = dataset) %>%
  ## add life stage and data type
  mutate(lifestage = "A",
         DataType = "Redd Survey", .after = dataset) %>%
  ## methods (no info at present)
  mutate(`Escapement Methodology Description` = NA,
         `Biologist Methodology Description` = NA, .after = DataType) %>%
  ## drop contingency table
  select(state:`Biologist Methodology Description`) %>%
  ## combine with other metadata
  rbind.data.frame(metadata) %>%
  ## sort by state and dataset
  arrange(state, dataset) %>%
  ## write to csv
  write_csv(file = file.path(raw_data_dir, metadata_file))


#### lookup table for site data ####

mt_dataset_lut <- mt_data_all %>%
  select(dataset:MT_ID)




