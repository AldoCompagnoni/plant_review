# Formatting functions for climate and lambda data-------------------------------------

# format species 
format_species <- function(spp_name, lam, response = "lambda"){
  
  # fetch what you need from 'lam' object
  lam_sel <- lam %>%
                subset( SpeciesAuthor == spp_name ) %>%
                dplyr::select( c("MatrixEndYear", "MatrixEndMonth", "MatrixPopulation", response) ) %>%
                setNames( c("year","month","population", response) ) %>%
                mutate( population = as.factor(population) )

  # list of pop-specific lambdas
  lam_l   <- lam_sel %>%
                dplyr::select( -population ) %>%
                split( lam_sel$population )
  
  return(lam_l)
  
}

# separate climate variables by population 
clim_list <- function(spp_name, clim, lam_spp){#clim_var, 
 
  clim_spp  <- clim %>%
                  subset( species == spp_name) %>%
                  mutate( population = as.factor(population) ) 
           
  # climate variables list
  clim_l    <- clim_spp %>%
                  dplyr::select(-population) %>%
                  split( clim_spp$population )
  
  return(clim_l)
  
}

# detrend population-level climate; put it in "long" form
clim_detrend <- function(clim_x, clim_var = "precip", st_dev = FALSE ){ 

  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y" ) 

  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
                  as.character %>%
                  as.data.frame(stringsAsFactors=F) %>%
                  separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                  bind_cols(clim_x) %>%
                  dplyr::select(-year,-day) %>%
                  setNames( c("year", "month", "day", "species", "ppt") )
    
  # # if climate_var airt, then do means, otherwise, do sums! 
  if( clim_var == "airt"){
    clim_m  <- clim_d %>%
      group_by(year, month) %>%
      summarise( ppt = mean(ppt, na.rm=T) )  %>%
      spread( month, ppt ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame()
  } else{
    clim_m   <- clim_d %>%
      group_by(year, month) %>%
      summarise( ppt = sum(ppt, na.rm=T) )  %>%
      spread( month, ppt ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame() #%>%
  }

  # if st_dev == T, this overrides the above conditional statements
  if(st_dev == T){
    clim_m   <- clim_d %>%
      group_by(year, month) %>%
      summarise( ppt = sd(ppt, na.rm=T) )  %>%
      spread( month, ppt ) %>%
      setNames( c("year",month.abb) ) %>%            
      as.data.frame() #%>%
  }
  
  # throw error
  if( !any( clim_var %in% c("precip","pet","airt","gdd")) ) {
    stop( paste0(clim_var," is not a supported varible") ) 
  }
  
  # detrend climate - but NOT if you are using GDD 
  if( clim_var != "gdd" ){
    d_clim   <- apply(clim_m[,-1], 2, FUN = scale, center = T, scale = T) %>%
                  as.data.frame() %>%
                  bind_cols( clim_m[,"year",drop=F] ) %>%
                  dplyr::select( c("year", month.abb) )
  }else{
    d_clim   <- clim_m
  }
  
  # Make NaNs 0
  for(c_i in 1:ncol(d_clim) ){
    d_clim[,c_i] <- replace(d_clim[,c_i], is.nan(d_clim[,c_i]), 0)
  } 
  
  return(d_clim)
  
}

# climate in long form 
clim_long <- function(clim_detr, lambda_data, m_back){
  
  # fecth year range, observation month
  years     <- lambda_data$year %>% unique
  yr_range  <- range(lambda_data$year)
  obs_month <- lambda_data$month %>% unique
  
  # detrended climate in "long" form
  long_out  <- clim_detr %>%
                  subset(year < (yr_range[2]+1) & year > (yr_range[1] - 6) ) %>%
                  gather(month, precip, Jan:Dec) %>%
                  setNames(c("year", "month", "clim_value")) %>% 
                  mutate(month_num = factor(month, levels = month.abb) ) %>% 
                  mutate(month_num = as.numeric(month_num)) %>% 
                  arrange(year, month_num)
  
  # select temporal extent
  clim_back <- function(yrs, obs_month, dat){
    
    obs_id    <- which(dat$year      == yrs & 
                       dat$month_num == obs_month)
    step_back <- m_back-1
    back_rows <- obs_id:(obs_id - step_back)
    
    return(dat[back_rows,"clim_value"])
  }

  # climate data in matrix form 
  mat_form<- function(dat, years){
    do.call(rbind, dat) %>% 
      as.data.frame %>%
      tibble::add_column(year = years, .before=1)
  }
  
  # arrange monthly climate values in wide form
  clim_obs  <- lapply(years, clim_back, obs_month, long_out)
  x_precip  <- mat_form(clim_obs, years)
  return(x_precip)
  
}

# combine climate data frames (if any)
lambda_plus_clim <- function(lambdas_l, clim_mat_l, response = "lambda"){
  
  # lambda and climate n. of populations correspond?
  if( length(lambdas_l) != length(clim_mat_l) ) stop("lambda and climate lists have differing lengths")
  
  # add population name to data frames
  population_add <- function(x, pop_name){
    tibble::add_column(x, population = pop_name)
  }
  lambdas_l   <- Map(population_add, lambdas_l,  names(lambdas_l) )
  clim_mat_l  <- Map(population_add, clim_mat_l, names(clim_mat_l) )
  
  # merge 
  if( length(lambdas_l) > 1){ # if n. of populations exceeds 1
    
    lambdas   <- Reduce(function(...) rbind(...), lambdas_l)
    climates  <- Reduce(function(...) rbind(...), clim_mat_l)
    clim_lam  <- merge(lambdas, climates)
    
  } else {
    
    clim_lam  <- merge(lambdas_l[[1]], clim_mat_l[[1]]) 
    
  }
  
  if( response == "lambda"){
    # order, and erase cases in which lambda == 0 (e.g. Eryngium_alpinum, BOU, year 2009)
    clim_lam    <- arrange(clim_lam, year, population)  %>%
                      subset( lambda != 0 ) 
    # erase any row containing NAs (for Dalgleish et al. 2010 data)
    r_id        <- lapply(clim_lam, function(x) which(is.na(x)) ) %>% unlist
    if( length(r_id) > 0 ) clim_lam  <- clim_lam[-r_id,]
    lam_out     <- dplyr::select(clim_lam, year:log_lambda)
    clim_out    <- dplyr::select(clim_lam, -c(year:log_lambda) )
    out         <- list(lambdas = lam_out, climate = clim_out)
  }else{
    # order, and erase cases in which lambda == 0 (e.g. Eryngium_alpinum, BOU, year 2009)
    clim_lam    <- arrange(clim_lam, year, population)
    # erase any row containing NAs (for Dalgleish et al. 2010 data)
    r_id        <- lapply(clim_lam, function(x) which(is.na(x)) ) %>% unlist
    if( length(r_id) > 0 ) clim_lam  <- clim_lam[-r_id,]
    eval(parse(n=1, text=paste0("lam_out <- dplyr::select(clim_lam, year:",response,")")))
    eval(parse(n=1, text=paste0("clim_out<- dplyr::select(clim_lam, -c(year:",response,"))")))
    out         <- list(lambdas = lam_out, climate = clim_out)
  }
  
  return(out)
  
}


# observed cliamtic range
observed_clim_range <- function(clim_x, lambda_d, spp_name, clim_var){

  # format day one
  day_one   <- as.Date( paste0("1/1/", first(clim_x$year) ), 
                        format="%d/%m/%Y" )
  
  # climate data
  clim_d    <- as.Date(1:nrow(clim_x), day_one-1) %>%
                  as.character %>%
                  as.data.frame(stringsAsFactors=F) %>%
                  separate_(col=".",into=c("year1","month1","day1"),sep="-") %>%
                  bind_cols(clim_x) %>%
                  dplyr::select(-year,-day) %>%
                  setNames( c("year", "month", "day", "species", "ppt") )
  
  # monthly climates
  # clim_var == "airt" calculate MEAN monthly air temperature 
  if( clim_var == "airt" ){
    clim_m   <- clim_d %>%
                    group_by(year, month) %>%
                    summarise( ppt = mean(ppt, na.rm=T) ) %>%
                    ungroup %>%
                    mutate( month = as.numeric(month) ) %>%
                    mutate( year = as.numeric(year) ) 
  }else{
    clim_m   <- clim_d %>%
                    group_by(year, month) %>%
                    summarise( ppt = sum(ppt, na.rm=T) ) %>%
                    ungroup %>%
                    mutate( month = as.numeric(month) ) %>%
                    mutate( year = as.numeric(year) ) 
  }
  
  # range of years
  max_yr    <- max(lambda_d$year)
  min_yr    <- min(lambda_d$year)
  month_i   <- unique(lambda_d$month)
  
  yearly_climate <- function(yrs){
    
    year_clim <- clim_m %>% subset( year == yrs & month < month_i + 1 )
    
    if( unique(lambda_d$month) != 12 ){
     
     append     <- clim_m %>% 
                      subset( year == (yrs-1) & month > month_i )
     year_clim  <- rbind(year_clim, append)
     
    }
  
    # Calculate means for "airt"
    if( clim_var == "airt"){
      out <- data.frame( year = yrs, ppt = mean(year_clim$ppt) )
    }else{
      out <- data.frame( year = yrs, ppt = sum(year_clim$ppt) )
    }
    
    return(out)
    
  }
  
  # yearly climates
  all_yrs   <- (max_yr-48):max_yr
  yr_clim_l <- lapply(all_yrs, yearly_climate)
  full_clim <- Reduce(function(...) rbind(...), yr_clim_l)
  obs_clim  <- subset(full_clim, year >= min_yr & year <= max_yr )
  
  # climate in full data set
  full_rng  <- max(full_clim$ppt) - min(full_clim$ppt)
  full_mean <- mean(full_clim$ppt)
  full_med  <- median(full_clim$ppt)
  full_sd   <- sd(full_clim$ppt)
  full_dev  <- range(abs(full_mean - full_clim$ppt))
  full_dev_r<- range(abs(full_med - full_clim$ppt))
  
  # observed range of climate anomalies
  obs_dev   <- range( abs(full_mean - obs_clim$ppt) )
  obs_dev_r <- range( abs(full_med - obs_clim$ppt) )
  obs_range <- max(obs_clim$ppt) - min(obs_clim$ppt)
  extr_yr_n <- sum(full_clim$ppt > max(obs_clim$ppt)) + sum(full_clim$ppt < min(obs_clim$ppt)) 
  obs_mean  <- mean(obs_clim$ppt)
  
  # proportions and means
  prop_yrs  <- (48-extr_yr_n) / 48 
  prop_rang <- obs_range / full_rng
  prop_var  <- (obs_dev[2] - obs_dev[1]) / (full_dev[2] - full_dev[1])
  prop_var_r<- (obs_dev_r[2] - obs_dev_r[1]) / (full_dev_r[2] - full_dev_r[1])
  mean_dev  <- (full_mean - obs_mean) / full_sd
  
  return( data.frame( species   = spp_name, 
                      prop_rang = prop_rang, 
                      prop_yrs  = prop_yrs,
                      prop_var  = prop_var,
                      prop_var_r= prop_var_r,
                      mean_dev  = mean_dev,
                      mean_clim = full_mean )
         )
  
}
  