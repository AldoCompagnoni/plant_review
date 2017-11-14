# create file with response variables (lambda) for plant_review species
setwd("C:/cloud/MEGA/Projects/sApropos")
options( stringsAsFactors = F )

# read 
all_demog <- read.csv("all_demog_6tr.csv")

# plant review species
pr_spp     <- c("Poa_secunda", "Hesperostipa_comata", 
                "Artemisia_tripartita", "Pseudoroegneria_spicata",
                "Helianthemum_juliae", "Purshia_subintegra",
                "Astragalus_cremnophylax_var._cremnophylax",
                "Opuntia_imbricata", 
                "Brassica_insularis", "Cryptantha_flava_2")

# species whose demographic information need be updated
spp_update <-  setdiff( pr_spp,
                        unique(all_demog$SpeciesAuthor) )


# Calculate lambdas 

# lambdas for Astragalus -------------------------------------------------------------------
load("COMPADRE_v.X.X.X.2.Rdata")
id        <- which(compadre$metadata$SpeciesAuthor == "Astragalus_cremnophylax_var._cremnophylax" &
                   compadre$metadata$MatrixComposite == "Individual" )

# calculate lambdas from ii indexes
lam_calc  <- function(ii){ 
  mat <- compadre$mat[[ii]]["matA"][[1]]
  Re(eigen(mat)$values[1])
}

# lambdas and metadata for Astragalus_cremnophylax_var._cremnophylax
# Matrix end month assumes sampling at flowering. Phenology from link below:
# https://www.fws.gov/southwest/es/arizona/Documents/SpeciesDocs/Sentry/Sentry%20milk-vetch%20facts.pdf
ascrcr  <- compadre$metadata[id,] %>%
            mutate( lambda = lapply(id, lam_calc) %>% unlist ) %>%
            mutate( MatrixEndMonth = 5)
fact_id <- which( sapply(ascrcr, class) == "factor" )

# convert factor variables to character 
for(ii in 1:length(fact_id) ){
  
  ascrcr[,fact_id[ii]] <- as.character( ascrcr[,fact_id[ii]] )
  
}


# lambdas for Cryptantha -------------------------------------------------------------------

# read Cryptantha data
crfl_raw  <- read.csv("results/compadre_updated_files/Cryptantha_flava_2.csv") %>%
                rename( Lat = LatNS,
                        Lon = LonWE) %>%
                mutate( Lat = 40.50000000,  
                        Lon = -109.36666667 )
crfl_sel  <- dplyr::select(crfl_raw, MatrixTreatment, MatrixEndYear) %>% unique
crfl_des  <- dplyr::select(crfl_raw, EnteredBy:MatrixClassOrganized) %>% unique

# test thast crfl_sel and crfl_des hold the same information
expect_true( all.equal(crfl_sel, 
                       dplyr::select(crfl_des, MatrixTreatment, MatrixEndYear) )
            )

# cal lambda by selecting
crfl_lam <- function(ii, x, crfl_des){
  
  # produce lambdas
  lambda <- subset(x, MatrixTreatment == crfl_des[ii,"MatrixTreatment"] &
                      MatrixEndYear   == crfl_des[ii,"MatrixEndYear"] ) %>%
              dplyr::select(A1:A10) %>%
              eigen %>%
              .$values %>%
              .[1] %>%
              Re

  # combine with metadata (or "design")
  out <- cbind(crfl_des[ii,], lambda)
  
  return( out )
    
}

crfl_lam_l  <- lapply(1:nrow(crfl_des), crfl_lam, crfl_raw, crfl_des)
crfl_lam_r  <- Reduce(function(...) rbind(...), crfl_lam_l)
# select relevant variables
sel_crfl    <- intersect(names(compadre$metadata), names(crfl_lam_r) )
crfl        <- dplyr::select(crfl_lam_r, sel_crfl, lambda)


# Species formatted by Aldo ----------------------------------------------------------------------

# general function: just provide file name
lambdas_by_file <- function(file_name){
  
  # read Cryptantha data
  dfrm_raw  <- read.csv( paste0("results/compadre_updated_files/", file_name) )
  dfrm_des  <- dplyr::select(dfrm_raw, EnteredBy:MatrixEndYear) %>% unique
  
  # cal lambda by selecting
  calc_lam <- function(ii, x, dfrm_des){
    
    mat_nam <- names(x)[ grep( "matA", names(x) ) ]
    
    # produce lambdas ( should crash if selection does not fetch a square matrix)
    lambda  <- subset(x, MatrixTreatment == dfrm_des[ii,"MatrixTreatment"] &
                         MatrixEndYear   == dfrm_des[ii,"MatrixEndYear"] ) %>%
                dplyr::select(mat_nam) %>%
                eigen %>%
                .$values %>%
                .[1] %>%
                Re
    
    # combine with metadata (or "design")
    out <- cbind(dfrm_des[ii,], lambda)
    
    return( out )
    
  }
  
  # data frame including lambda values 
  dfrm_lam_l  <- lapply(1:nrow(dfrm_des), calc_lam, dfrm_raw, dfrm_des)
  dfrm_lam_r  <- Reduce(function(...) rbind(...), dfrm_lam_l) %>%
                    rename( Lat = LatNS,
                            Lon = LonWE)
  # select relevant variables
  sel_dfrm    <- intersect(names(compadre$metadata), names(dfrm_lam_r) )
  out         <- dplyr::select(dfrm_lam_r, sel_dfrm, lambda)
  
  return(out)

}

aldo_form_spp   <- pr_spp[c(1:4,8)]
lam_aldo_spp_l  <- lapply(paste0(aldo_form_spp,".csv"), lambdas_by_file)
lam_aldo_spp_df <- Reduce(function(...) rbind(...), lam_aldo_spp_l)


# put it all together -------------------------------------------------------------- 

# new species for plant review 
new_lambdas <- Reduce(function(...) bind_rows(...), list(lam_aldo_spp_df, ascrcr, crfl) )
new_lambdas <- new_lambdas %>% mutate( log_lambda = log(lambda) )

# add species already in all_demog object
spp_add    <- setdiff(pr_spp, unique(new_lambdas$SpeciesAuthor) )
common_var <- intersect(names(new_lambdas),
                        names(all_demog ) )
spp_add_df <- all_demog %>%
                subset( SpeciesAuthor %in% spp_add ) %>%
                dplyr::select( common_var )

# put it all together
pr_lambdas <- rbind(new_lambdas, spp_add_df )

write.csv(pr_lambdas, "data_plant_review/pr_lambdas.csv", row.names=F)
