rm(list=ls())

# Load libraries
library(tidyverse) # lots of useful functions

# Set working directory
setwd("C:/Users/User/OneDrive - M365 Universität Hamburg/4_DesmidFlowCam/analysis")

# Load data
raw_dat_meso <- read_csv("desmid_counts_meso.csv")
raw_dat_mo <- read_csv("desmid_counts_mo.csv")


# Tidy - pivot long
rawdat_long_meso <- raw_dat_meso %>% 
  pivot_longer(cols = 7:9,
               names_to = "species",
               values_to = "n")

rawdat_long_mo <- raw_dat_mo %>% 
  pivot_longer(cols = 7:9,
               names_to = "species",
               values_to = "n")

rawdat_long <- rbind(rawdat_long_meso, rawdat_long_mo)

# Calculations
dat <- rawdat_long %>% 
  # calc total desmids identified in sample
  group_by(day, community, media, temp, replicate) %>% 
  mutate(total_desmids = sum(n)) %>% 
  ungroup() %>% 
  # calc average species per triplicate
  group_by(day, community, media, temp, species) %>% 
  mutate(sp_mean = mean(n),
         total_desmids_mean = mean(total_desmids)) %>% 
  ungroup() %>% 
  # calc relative abundance of each species
  mutate(rel_abund = (sp_mean/total_desmids_mean)*100)


# sanity check that rel. abund = 100%
dat %>% 
  group_by(day, community, media, temp, replicate) %>% 
  mutate(pct_chk = sum(rel_abund)) %>% 
  ungroup() %>% 
  # select(day, sample, replicate, rel_abund, pct_chk) %>% 
  distinct(pct_chk, .keep_all = TRUE) %>% # only 100 -- good!
  select(1:7, 13) %>% 
  head()

head(dat, 20)

# Tidy
str(dat)

dat$day <- as.factor(dat$day)
dat$temp <- as.factor(dat$temp)
dat$replicate <- as.factor(dat$replicate)

glimpse(dat)
head(dat)

write.csv(dat, "flowcam_dat.csv", row.names = FALSE)
