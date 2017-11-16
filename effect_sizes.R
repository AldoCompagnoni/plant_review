setwd("C:/cloud/MEGA/Projects/sApropos/results_plant_review")
library(tidyverse)

# file lists
airt_f    <- Filter(function(x) grepl("mod_summaries",x), list.files("airt") )
prec_f    <- Filter(function(x) grepl("mod_summaries",x), list.files("precip") )

# spp lists
airt_spp  <- gsub("mod_summaries_|.csv","",airt_f) 
prec_spp  <- gsub("mod_summaries_|.csv","",prec_f) 

# data frame
format_dfs <- function(x,y,clim_var){
  read.csv(paste0(clim_var,"/",x)) %>%
    mutate( species = y )
}

# list of data frames
airt_l  <- Map(format_dfs, airt_f, airt_spp, "airt")
prec_l  <- Map(format_dfs, prec_f, prec_spp, "precip")

# un-normalized data frames
airt_df <- Reduce(function(...) rbind(...), airt_l)
prec_df <- Reduce(function(...) rbind(...), prec_l)



# effect sizes species by specis ----------------------------------------------------------

# Astragalus_cremnophylax_var._cremnophylax 
ascr      <- subset(prec_df,
                    species == "Astragalus_cremnophylax_var._cremnophylax") %>%
                    dplyr::select(model, alpha_mean, beta_mean)

# Brassica_insularis
brin_prec <- subset(prec_df,
                species == "Brassica_insularis") %>%
                dplyr::select(model, alpha_mean, beta_mean)

brin_airt <- subset(airt_df,
                    species == "Brassica_insularis") %>%
                    dplyr::select(model, alpha_mean, beta_mean)

# Opuntia_imbricata
opim_airt <- subset(airt_df,
                    species == "Brassica_insularis") %>%
                    dplyr::select(model, alpha_mean, beta_mean)

# Cryptantha_flava_2
crfl      <- subset(prec_df,
                    species == "Cryptantha_flava_2") %>%
                dplyr::select(model, alpha_mean, beta_mean)

# Purshia_subintegra
pusu      <- subset(prec_df,
                    species == "Purshia_subintegra") %>%
                dplyr::select(model, alpha_mean, beta_mean)

# Helianthemum_juliae
heju_prec <- subset(prec_df, 
                    species == "Helianthemum_juliae") %>%
                dplyr::select(model, alpha_mean, beta_mean)

heju_airt <- subset(airt_df, 
                    species == "Helianthemum_juliae") %>%
                dplyr::select(model, alpha_mean, beta_mean)

# Daphne_rodriguezii
daro      <- subset(prec_df,
                    species == "Daphne_rodriguezii") %>%
                dplyr::select(model, alpha_mean, beta_mean)
