setwd("C:/cloud/MEGA/Projects/sApropos")
library(dplyr)
library(tidyr)
library(testthat)
load("COMPADRE_v.X.X.X.2.Rdata")

# plant review species
pr_spp <- c("Poa_secunda", "Hesperostipa_comata", 
            "Artemisia_tripartita", "Pseudoroegneria_spicata",
            "Helianthemum_juliae", "Purshia_subintegra",
            "Astragalus_cremnophylax_var._cremnophylax",
            "Opuntia_imbricata", 
            "Brassica_insularis", "Cryptantha_flava_2")
  
future_spp <- "Dracocephalum_austriacum"


tmp <- subset(compadre$metadata , 
              SpeciesAuthor == "Astragalus_cremnophylax_var._cremnophylax" &
              MatrixComposite == "Individual" )$MatrixPopulation


# species lat lon ------------------------------------------------------------------
grouped_data    <- data.frame(  SpeciesAuthor  = pr_spp,
                                MatrixPopulation = c( rep("Dubois sheep station",4),
                                                      NA, NA,   
                                                      "Grand Canyon National Park",
                                                      "Sevilleta National Wildlife Refuge",
                                                      NA,
                                                      "Redfleet State Park"),
                                lat      = c( rep(44.2,  4), NA, NA,   36.1007,   34.33480556,
                                              NA, 40.500  ),
                                lon      = c( rep(-112.1,4), NA, NA, -112.1013, -106.63138889,
                                              NA, -109.375),
                                end_year = c( rep(1957,  4),   NA, NA,    1994,          2013,
                                              NA,     2012)
                              ) %>%
                      subset( ( !is.na(MatrixPopulation) & 
                                 MatrixPopulation != "Dubois sheep station") )


# (nested) functions to fecth climate --------------------------------------------

# fetch daily climate
fetch_daily_clim <- function(yr, var, sp){
  
  fc_raw    <- fcTimeSeriesDaily(variable = var,
                                 latitude = sp$lat, longitude = sp$lon,
                                 firstYear = yr, lastYear = yr)
  
  # format data into a data frame
  fc_out    <- data.frame(species = as.character(sp$SpeciesAuthor),
                          population = sp$MatrixPopulation,
                          year = as.integer(yr), 
                          day = as.integer(fc_raw$days), 
                          clim_var = as.numeric(fc_raw$values[1,]),
                          stringsAsFactors = F)
  
  # change name of variable 
  fc_out    <- rename_(fc_out, .dots = setNames("clim_var", var) )
  
  return(fc_out)
  
}

# fetch climate across species
climate_spp <- function(sp_i, var, grouped_data, yr_back){
  
  sp      <- grouped_data[sp_i,]
  yr_r    <- seq(sp$end_year-yr_back, sp$end_year, by = 1)
  tmp     <- lapply(yr_r, fetch_daily_clim, var, sp)
  sp_clim <- Reduce(function(...) rbind(...), tmp) %>%
    as.data.frame(stringsAsFactors = F)
  return(sp_clim)
  
}

# download data separately (computer crashes otherwise) --------------------------------------

# air temperature data
spp_airt  <- lapply(1:3, climate_spp, "airt", grouped_data, 49)
airt_1    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_1, "C:/cloud/MEGA/Projects/sApropos/airt_fc_plant_rev.csv", row.names = F)

# precipitation data
spp_airt  <- lapply(1:3, climate_spp, "prate", grouped_data, 49)
airt_1    <- Reduce(function(...) rbind(...), spp_airt)
write.csv(airt_1, "C:/cloud/MEGA/Projects/sApropos/precip_fc_plant_rev.csv", row.names = F)



# Idaho data sets ----------------------------------------------------------------------------

# extract string characters
extr <- function(string,pattern){
  regmatches(string, regexpr(pattern,string))
}

# idaho data
dubois_df <- read.csv("data/idaho_matrices/climate/UCC_ghcn_USC00102707_2017_11_10_1510325677.csv",
                      stringsAsFactors = F, skip=15)


# introduce NAs --------------------------------------------------------------------

# change to character
dubois_df[] <- lapply(dubois_df[], function(x) x %<>% as.character(x) )

# substitute to NA
replace_char <- function(x, field){ x <- replace(x, x == field, NA) } 

# na flags
na_flags     <- sapply(dubois_df, unique) %>% 
  unique %>% 
  unlist %>% 
  extr("[A-Z]") %>%
  unique

# introduce NAs
dubois_df[] <- lapply(dubois_df[],  replace_char, "T")
dubois_df[] <- lapply(dubois_df[],  replace_char, "M")

# rest that flags have disappeared
na_flags    <- sapply(dubois_df, unique) %>% 
  unique %>% 
  unlist %>% 
  extr("[A-Z]") %>%
  unique
# there should be 0 flags 
expect_equal(length(na_flags), 0 )


# format dates --------------------------------------------------------------

# format Data in three separate columns
dubois_df   <- separate(dubois_df, Day, c("year","month","day"), sep = "-")

# make everything "numeric"
dubois_df[] <- lapply(dubois_df[], function(x) x %<>% as.numeric(x) )

# Make day vary from 1 to 365-366
dubois_l    <- split(dubois_df, as.factor(dubois_df$year) )
dubois_df_l <- lapply(dubois_l, function(x) {
  x$day <- c(1:nrow(x))
  return(x)
} )
dubois_df   <- Reduce(function(...) rbind(...), dubois_df_l)


# "sub in" dubois data into fetch climate, FOR DAILGLEISH SPP. --------------------

# calculate airt (mean daily temperature), PET, precip
dubois_df   <- dubois_df %>% 
  dplyr::select(-Snow.Depth,-Snow.Fall) %>%
  mutate( airt = (Min.Temperature + Max.Temperature) / 2 ) %>%
  dplyr::select(-Min.Temperature,-Max.Temperature) %>% 
  rename( ppt = Precipitation,
          pet = Ref.Evapotranspiration ) %>%
  mutate( population = "Dubois sheep station" ) %>%
  dplyr::select( -month ) 


# ciao ----------------------------------------------------------------------------

# introduce species names
intro_spp <- function(spp_n, x){
  
  x %>%
    mutate( species = spp_n ) %>%
    dplyr::select(species, population, year, day, everything() )
  
}

# intro spp names
dubois_o_l <- lapply( c("Artemisia_tripartita",
                        "Pseudoroegneria_spicata",
                        "Hesperostipa_comata",
                        "Poa_secunda"),
                      intro_spp, dubois_df )

# test that all data frames have same length
expect_equal( lapply(dubois_o_l, nrow) %>% 
                unlist %>%
                unique %>% 
                length, 1)

# full climate data for dubois
dubois_out  <- Reduce(function(...) rbind(...), dubois_o_l)


# write -----------------------------------------------------------------

# air temperature data
data.table::fwrite(dplyr::select(dubois_out,-airt,-pet), 
                   "precip_dubois.csv")

# air temperature data
data.table::fwrite(dplyr::select(dubois_out,-ppt,-pet), 
                   "airt_dubois.csv")


# put all climate data together -----------------------------------------

# plant review airt data
pr_new    <- read.csv("airt_fc_plant_rev.csv")
pr_dubois <- read.csv("airt_dubois.csv")
pr_airt   <- rbind(pr_dubois,  pr_new)

# plant review precip data
pr_new    <- read.csv("precip_fc_plant_rev.csv") %>% rename( ppt = prate )
pr_dubois <- read.csv("precip_dubois.csv")
pr_precip <- rbind(pr_dubois,  pr_new)


# TEST: Compare species with the most recent climate file
all_spp  <- data.table::fread("precip_fc_hays.csv") %>% .$species %>% unique %>% sort

# This should be equal to the list of our species
pr_spp_l <- pr_airt$species %>% unique %>% as.character %>% sort

expect_equal( setdiff(pr_spp_l, all_prec$species),
              pr_spp_l )


# store updated climatic data --------------------------------------------------------

# precipitation
prec_all  <- data.table::fread("precip_fc_hays.csv") 
prec_updt <- rbind(prec_all, pr_precip) %>% arrange(species, population, year, day)

# air temperature
airt_all  <- data.table::fread("airt_fc_hays.csv") 
airt_updt <- rbind(airt_all, pr_airt)   %>% arrange(species, population, year, day)

# store results 
data.table::fwrite(prec_updt, "precip_fc_11.10.17.csv")
data.table::fwrite(airt_updt, "airt_fc_11.10.17.csv")
