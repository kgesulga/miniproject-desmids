# Set working directory
setwd("C:/Users/User/Desktop/Desmids/3_experiment/4_photometer")


# Load libraries
library(tidyverse)
library(cowplot) # to arrange plots
library(calecopal) # pretty palettes!


# List all .Rdata files 
my_files <- list.files(path = "data/", pattern = "*.RData", full.names = T)
# Load separate data frames into environment
load_data <- lapply(my_files, load, .GlobalEnv)
# Bind all rows into one main data frame
dat_all <- bind_rows(mget(unlist(load_data)))
# Clean up environment
rm(list=ls(pattern="T"))

# save(dat_all, file = "../../4_analyses/data/dat_all.RData")

dat_growthrate <- dat_all %>% 
  filter(Plate == "Fluorescence Ex 475",
         Wavelength == 690) %>% 
  filter(!is.na(community)) %>% 
  select(day, temp_treatment_C, community, media, replicate, value) %>% 
  rename("Day" = "day",
         "Temp" = "temp_treatment_C",
         "Nutrient" = "media",
         "BioRep" = "replicate",
         "PFU" = "value") # %>% 
  # mutate(Temp = case_when(Temp == 18 ~ "18°C",
  #                         Temp == 22 ~ "22°C"))

  
write.csv(dat_growthrate, file = "../../4_analyses/data/dat_growthrate_noedits.csv", row.names = FALSE)
