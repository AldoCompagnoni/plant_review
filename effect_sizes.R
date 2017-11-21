setwd("C:/cloud/MEGA/Projects/sApropos/results_plant_review")
library(tidyverse)
options(stringsAsFactors = F)

# literature effect sizes
lit_eff   <- read.csv("C:/cloud/MEGA/Projects/sApropos/data_plant_review/lit_effect_sizes_pr.csv")

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

# to "pipe" Reduce + rbind
reduce_rbind_pipe <- function(x){
  Reduce(function(...) rbind(...),x)
}

# list of data frames
airt_df  <- Map(format_dfs, airt_f, airt_spp, "airt") %>% reduce_rbind_pipe
prec_df  <- Map(format_dfs, prec_f, prec_spp, "precip") %>% reduce_rbind_pipe


# effect sizes species by specis ----------------------------------------------------------

# best model indexes
best_id <- which( (c(1:14) %% 2) == 1)

# select betas of best moving window model
beta_get <- function(x, best_id, clim_var){
  
  x %>% 
    subset( (model == "expp" | model == "gaus") ) %>% 
    dplyr::select(species, beta_mean) %>%
    rename_( .dots = setNames("beta_mean",paste0(clim_var,"_beta")) ) %>%
    .[best_id,]
    
}

prec_betas <- beta_get(prec_df, best_id, "prec")
airt_betas <- beta_get(airt_df, best_id, "airt")

comp_tab   <- Reduce(function(...) left_join(...), list(prec_betas,airt_betas, lit_eff) ) %>%
                dplyr::select(species, lit_effect_prec, prec_beta,
                                       lit_effect_airt, airt_beta)

write.csv(comp_tab, "compare_beta_pr.csv", row.names=F)


# airt_betas <- airt_df %>% 
#                 subset( (model == "expp" | model == "gaus") ) %>% 
#                 dplyr::select(species, model, beta_mean) 
# 
# prec_betas <- prec_df %>% 
#                 subset( (model == "expp" | model == "gaus") ) %>% 
#                 dplyr::select(species, model, beta_mean)

# 
# # Astragalus_cremnophylax_var._cremnophylax 
# ascr      <- subset(prec_df,
#                     species == "Astragalus_cremnophylax_var._cremnophylax") %>%
#                     dplyr::select(model, alpha_mean, beta_mean)
# 
# # Brassica_insularis
# brin_prec <- subset(prec_df,
#                 species == "Brassica_insularis") %>%
#                 dplyr::select(model, alpha_mean, beta_mean)
# 
# brin_airt <- subset(airt_df,
#                     species == "Brassica_insularis") %>%
#                     dplyr::select(model, alpha_mean, beta_mean)
# 
# # Opuntia_imbricata
# opim_airt <- subset(airt_df,
#                     species == "Brassica_insularis") %>%
#                     dplyr::select(model, alpha_mean, beta_mean)
# 
# # Cryptantha_flava_2
# crfl      <- subset(prec_df,
#                     species == "Cryptantha_flava_2") %>%
#                 dplyr::select(model, alpha_mean, beta_mean)
# 
# # Purshia_subintegra
# pusu      <- subset(prec_df,
#                     species == "Purshia_subintegra") %>%
#                 dplyr::select(model, alpha_mean, beta_mean)
# 
# # Helianthemum_juliae
# heju_prec <- subset(prec_df, 
#                     species == "Helianthemum_juliae") %>%
#                 dplyr::select(model, alpha_mean, beta_mean)
# 
# heju_airt <- subset(airt_df, 
#                     species == "Helianthemum_juliae") %>%
#                 dplyr::select(model, alpha_mean, beta_mean)
# 
# # Daphne_rodriguezii
# daro      <- subset(prec_df,
#                     species == "Daphne_rodriguezii") %>%
#                 dplyr::select(model, alpha_mean, beta_mean)
