#bjtwd<-"C:/Users/admin_bjt162/Dropbox/A.Current/Ongoing_Collab_Research/sApropos project/"
rm(list=ls())
setwd("C:/cloud/MEGA/Projects/sApropos/")
source("C:/CODE/plant_review/format_data.R")
library(tidyverse)
library(dismo)
library(mgcv)
library(testthat)
library(rstan)
library(loo)

# set rstan options
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )

# arguments from command line
args <- commandArgs(TRUE)

# climate predictor, response, months back, max. number of knots
clim_var <- args[4]
response <- args[5]
m_back   <- 24
expp_beta<- 20

# read data -----------------------------------------------------------------------------------------
lam       <- read.csv( paste0(args[2], 'pr_lambdas.csv'), stringsAsFactors = F)
clim      <- read.csv( paste0(args[2], clim_var,"_fc_11.10.17.csv"),  stringsAsFactors = F)
spp       <- lambdas$SpeciesAuthor %>% unique

# format data --------------------------------------------------------------------------------------

# set up model "family" based on response
if( response == "log_lambda") family = "normal"

# set species (I pick Sphaeraclea_coccinea)
ii            <- as.numeric(args[1])
spp_name      <- spp[ii]

# lambda data
spp_lambdas   <- format_species(spp_name, lambdas, response)

# climate data
clim_separate <- clim_list(spp_name, clim, spp_lambdas)
clim_detrnded <- lapply(clim_separate, clim_detrend, clim_var, st_dev)
clim_mats     <- Map(clim_long, clim_detrnded, spp_lambdas, m_back)

# model data
mod_data          <- lambda_plus_clim(spp_lambdas, clim_mats, response)
mod_data$climate  <- mod_data$climate #/ diff(range(mod_data$climate))

# throw error if not enough data
if( nrow(mod_data$lambdas) < 6 ) stop( paste0("not enough temporal replication for '", 
                                              spp_name, "' and response variable '" , response, "'") )

# Fit models ----------------------------------------------------------------------------------------

# organize data into list to pass to stan
dat_stan <- list(
  n_time     = nrow(mod_data$climate),
  n_lag      = ncol(mod_data$climate),
  y          = mod_data$lambdas[,response],
  clim       = mod_data$climate,
  clim_means = rowMeans(mod_data$climate),
  expp_beta  = expp_beta
)

# data for the 12 month and 13-24 months (previous year) models
dat_stan_yrt      <- dat_stan 
dat_stan_yrt1     <- dat_stan

dat_stan_yrt$clim        <- mod_data$climate[,1:12]
dat_stan_yrt$clim_means  <- rowMeans(mod_data$climate[,1:12])

dat_stan_yrt1$clim       <- mod_data$climate[,13:24]
dat_stan_yrt1$clim_means <- rowMeans(mod_data$climate[,13:24])


# simulation parameters
sim_pars <- list(
  warmup = 1000, 
  iter = 4000, 
  thin = 2, 
  chains = 4
)

# NULL model (model of the mean)
fit_ctrl1 <- stan(
  file = paste0("stan/",family,"_null.stan"),
  data = dat_stan,
  pars = c('alpha', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# average climate model
fit_ctrl2 <- stan(
  file = paste0("stan/",family,"_ctrl2.stan"),
  data = dat_stan,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# average climate model
fit_yrt <- stan(
  file = paste0("stan/",family,"_ctrl2.stan"),
  data = dat_stan_yrt,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# average climate model
fit_yrt1 <- stan(
  file = paste0("stan/",family,"_ctrl2.stan"),
  data = dat_stan_yrt1,
  pars = c('alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# gaussian moving window
fit_gaus <- stan(
  file = paste0("stan/",family,"_gaus.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# exponential power moving window
fit_expp <- stan(
  file = paste0("stan/",family,"_expp.stan"),
  data = dat_stan,
  pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'log_lik'),
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains#,
  #control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth = 20)
)

# parameter values and diagnostics ----------------------------------------------------------------

# list of model fits
mod_fit   <- list( ctrl1 = fit_ctrl1, ctrl2 = fit_ctrl2,
                   yr_t  = fit_yrt,   yr_t1 = fit_yrt1,
                   gaus  = fit_gaus,  expp  = fit_expp )

# parameter values
pars      <- c('sens_mu', 'sens_sd', 'alpha', 'beta', "y_sd")

# get central tendencies
pars_diag_extract <- function(x){
  
  # central tendencies
  tmp         <- rstan::extract(x)
  par_means   <- sapply(tmp, function(x) mean(x)) %>%
                    setNames( paste0(names(tmp),"_mean") )
  par_medians <- sapply(tmp, function(x) median(x)) %>%
                    setNames( paste0(names(tmp),"_median") )
  central_tend<- c(par_means, par_medians)
  
  # diagnostics
  diverg      <- do.call(rbind, args = get_sampler_params(x, inc_warmup = F))[,5]
  n_diverg    <- length(which(diverg == 1))
  df_summ     <- as.data.frame(summary(x)$summary)
  rhat_high   <- length(which(df_summ$Rhat > 1.1))
  n_eff       <- df_summ$n_eff / length(diverg)
  n_eff_low   <- length(which(n_eff < 0.1))
  mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  diagnostics <- c(n_diverg = n_diverg, rhat_high = rhat_high,
                   n_eff_low = n_eff_low, mcse_high = mcse_high)
  out         <- c( central_tend, diagnostics ) %>% t %>% as.data.frame
  
  rm(tmp) ; return(out)
  
}

# store posteriors
posterior_extract <- function(model_fit, model_name){
  
  # central tendencies
  tmp         <- rstan::extract(model_fit)
  post_df     <- do.call(cbind, tmp) %>% as.data.frame
  ll_id       <- grep("V", colnames(post_df) )
  new_names   <- paste0("log_lik_", 1:length(ll_id) )
  names(post_df)[ll_id] <- new_names # no way to do this in dplyr
  post_df     <- tibble::add_column(post_df,
                                    model = model_name, .before=1)
  
  rm(tmp) ; return(post_df)
  
}

# calculate central tendencies
pars_diag_l   <- lapply(mod_fit, pars_diag_extract)
mod_pars_diag <- Reduce(function(...) bind_rows(...), pars_diag_l) %>%
                    tibble::add_column(model = names(mod_fit), .before = 1)

# store posteriors
posts_l       <- Map(posterior_extract, mod_fit, names(mod_fit) )
posteriors    <- Reduce(function(...) bind_rows(...), posts_l)


# WAIC model comparison --------------------------------------------------------------------

# wAIC model selection using loo approximation (from library 'loo')
log_liks  <- lapply(mod_fit, extract_log_lik)
names(log_liks)
# leave-one-out estimates
loo_l      <- lapply(log_liks, loo) %>%
                setNames( c("loo_ctrl1",  "loo_ctrl2", "loo_yr_t",  "loo_yr_t1", "loo_gaus", "loo_expp") )
loo_df     <- loo::compare(loo_l$loo_ctrl1, loo_l$loo_ctrl2, 
                           loo_l$loo_yr_t, loo_l$loo_yr_t1,
                           loo_l$loo_gaus, loo_l$loo_expp) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("loo_","",names(loo_l) ), .before = 1)

# WAIC estimates
waic_l    <- lapply(log_liks, waic) %>%
                setNames(c("waic_ctrl1", "waic_ctrl2", 
                           "waic_yr_t",  "waic_yr_t1", 
                           "waic_gaus",  "waic_expp") )
waic_df   <- loo::compare(waic_l$waic_gaus, waic_l$waic_expp, 
                          waic_l$waic_yr_t, waic_l$waic_yr_t1, 
                          waic_l$waic_ctrl1, waic_l$waic_ctrl2) %>%
                as.data.frame %>%
                tibble::add_column(model = gsub("waic_","",names(waic_l) ), .before = 1)



# leave-one-out crossvalidation ------------------------------------------------------------------------

# crossvalidation function
CrossVal <- function(i, mod_data, response) {       # i is index for row to leave out
  
  # identify years
  uniq_yr           <- mod_data$lambdas$year %>% unique 
  test_i            <- which(mod_data$lambdas$year == uniq_yr[i])
  
  # put all in matrix form 
  x_clim            <- mod_data$climate
  x_clim_means      <- rowMeans(mod_data$climate)   # climate averages over entire window (for control model #2)
  
  # response variable
  y_train           <- mod_data$lambdas[-test_i, response]
  y_test            <- mod_data$lambdas[test_i, response]
  
  # climate variable
  clim_train        <- x_clim[-test_i,]
  clim_test         <- x_clim[test_i,] 
  
  # climate averages over full 24-month window (for control model #2)
  clim_means_train  <- x_clim_means[-test_i]
  clim_means_test   <- x_clim_means[test_i] 
  
  # data into a list for stan ----------------------------------------------------------------------------
  dat_stan_crossval <- list(
    n_train = length(y_train),  # number of data points in train set (length of response var)
    n_test  = length(y_test),   # number of data points in test set
    n_lag   = ncol(clim_train), # maximum lag
    y_train = array(y_train),
    y_test  = array(y_test),
    clim_train       = array(clim_train),
    clim_test        = array(clim_test),
    clim_means_train = array(clim_means_train), # climate averages over full 24-month window (for control model #2)
    clim_means_test  = array(clim_means_test),   # climate averages over full 24-month window (for control model #2)
    expp_beta        = expp_beta       # beta paramater for exponential power distribution
  )
  
  # data for "year t' and "t-1" ----------------------------------------------------------------------------
  
  # data for the 12 month
  dat_stan_yr_t        <- dat_stan_crossval 
  
  # climate averages for 12 months
  x_clim_means_t      <- rowMeans(mod_data$climate[,1:12])
  
  # climate averages over previous year
  dat_stan_yr_t$clim_means_train  <- array(x_clim_means_t[-test_i])
  dat_stan_yr_t$clim_means_test   <- array(x_clim_means_t[test_i])
  
  # data for the 13-24 months (previous year) models
  dat_stan_yr_t1       <- dat_stan_crossval 
  
  # climate averages for 12 months
  x_clim_means_t1     <- rowMeans(mod_data$climate[,13:24])
  
  # climate averages over previous year
  dat_stan_yr_t1$clim_means_train  <- array(x_clim_means_t1[-test_i])
  dat_stan_yr_t1$clim_means_test   <- array(x_clim_means_t1[test_i])
  
  
  # fit control 1 (intercept only)
  fit_ctrl1_crossval <- stan(
    file = paste0("stan/",family,"_null_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit control 2 (full window climate average)
  fit_ctrl2_crossval <- stan(
    file = paste0("stan/",family,"_ctrl2_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit current year mean climate
  fit_yr_t_crossval <- stan(
    file = paste0("stan/",family,"_ctrl2_crossval.stan"),
    data = dat_stan_yr_t,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit previous year mean climate
  fit_yr_t1_crossval <- stan(
    file = paste0("stan/",family,"_ctrl2_crossval.stan"),
    data = dat_stan_yr_t1,
    pars = c('alpha', 'beta', 'y_sd', 'pred_y', 'log_lik', 'log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains
  )
  
  # fit moving window, gaussian
  fit_gaus_crossval <- stan( 
    file = paste0("stan/",family,"_gaus_crossval.stan"), 
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.005, max_treedepth = 12)
  )
  
  # fit moving window, exponential power
  fit_expp_crossval <- stan( 
    file = paste0("stan/",family,"_expp_crossval.stan"),
    data = dat_stan_crossval,
    pars = c('sens_mu', 'sens_sd', 'alpha', 'beta', 'y_sd', 'pred_y', 'log_lik','log_lik_test'),
    warmup = sim_pars$warmup,
    iter = sim_pars$iter,
    thin = sim_pars$thin,
    chains = sim_pars$chains#,
    #control = list(adapt_delta = 0.999, stepsize = 0.1, max_treedepth = 12)
  )
  
  # posterior mean prediction for the out-of-sample value
  crossval_mods <- list( ctrl1 = fit_ctrl1_crossval, 
                         ctrl2 = fit_ctrl2_crossval,
                         yr_t  = fit_yr_t_crossval,
                         yr_t1 = fit_yr_t1_crossval,
                         gaus  = fit_gaus_crossval, 
                         expp  = fit_expp_crossval )
  
  # predictions
  mod_preds <- lapply(crossval_mods, function(x) rstan::extract(x, 'pred_y')$pred_y %>% apply(2,mean) )
  
  # Expected Log Predictive Density
  mod_elpds <- lapply(crossval_mods, function(x) rstan::extract(x, 'log_lik_test')$log_lik_test %>% 
                                                 apply(2,mean) )
  
  # diagnostics 
  diagnostics <- function(fit_obj, name_mod){
    
    diverg      <- do.call(rbind, args = get_sampler_params(fit_obj, inc_warmup = F))[,5]
    n_diverg    <- length(which(diverg == 1))
    df_summ     <- as.data.frame(summary(fit_obj)$summary)
    rhat_high   <- length(which(df_summ$Rhat > 1.1))
    n_eff       <- df_summ$n_eff / length(diverg)
    n_eff_low   <- length(which(n_eff < 0.1))
    mcse_high   <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
    out         <- data.frame(n_diverg, rhat_high, 
                              n_eff_low, mcse_high) 
    out         <- setNames(out, paste0(names(out),"_",name_mod) )
    return(out)
    
  }
  
  # store diagnostics
  gaus_expp   <- crossval_mods[c("gaus","expp")]
  diagnost_l  <- Map(diagnostics, gaus_expp, names(gaus_expp))
  diagnost_df <- do.call(cbind, diagnost_l) %>%
                    bind_cols( unique( dplyr::select(mod_data$lambdas[test_i,],year) ) )
  
  # function
  pred_elpd_df<- mod_data$lambdas[test_i,] %>%
                    mutate( # predictions
                            gaus_pred  = mod_preds$gaus,
                            expp_pred  = mod_preds$expp,
                            ctrl1_pred = mod_preds$ctrl1,
                            ctrl2_pred = mod_preds$ctrl2,
                            yr_t_pred  = mod_preds$yr_t,
                            yr_t1_pred = mod_preds$yr_t1,
                            # Expected Log Predictive Density
                            gaus_elpd  = mod_elpds$gaus,
                            expp_elpd  = mod_elpds$expp,
                            ctrl1_elpd = mod_elpds$ctrl1,
                            ctrl2_elpd = mod_elpds$ctrl2,
                            yr_t_elpd  = mod_elpds$yr_t,
                            yr_t1_elpd = mod_elpds$yr_t1 )
                  
  # df to return
  out         <- left_join(pred_elpd_df, diagnost_df)
  
  # remove stanfit objects (garbage collection)
  rm(fit_gaus_crossval )
  rm(fit_expp_crossval )
  rm(fit_ctrl1_crossval)
  rm(fit_ctrl2_crossval)
  rm(fit_yr_t_crossval )
  rm(fit_yr_t1_crossval)
  
  return(out)
  
}

# spp-specific cross validation
year_inds   <- seq_along(unique(mod_data$lambdas$year))
cxval_res   <- lapply( year_inds, CrossVal, mod_data, response)
cxval_pred  <- do.call(rbind, cxval_res) 

# measures of fit -------------------------------------------------------------------------- 

# calculate either mse or deviance
pred_perform <- function(x, mod_data, response, type){
  
  if( type == "mse"){
    res   <- (x - mod_data$lambdas[,response])^2 %>% mean
  }
  if(type == "deviance"){
    res   <-calc.deviance(x, mod_data$lambdas[,response], 
                          weights = rep(1, length(x) ),  
                          family="gaussian", calc.mean = TRUE)
  }
  
  return(res)
  
}

# format results into a data frame
perform_format <- function(x, var){
  
  x %>%
    unlist %>%
    t %>% 
    t %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "model") %>%
    mutate( model = gsub("mod_preds.", "", model))  %>%
    setNames( c("model", var) )
  
}

# mean squared error
mse <- cxval_pred %>% 
          dplyr::select(gaus_pred:yr_t1_pred) %>%
          lapply(pred_perform, mod_data, response, "mse") %>%
          perform_format("mse") %>%
          mutate( model = gsub("_pred","",model) )

# # deviance 
# devi <- cxval_pred %>% 
#           dplyr::select(gaus_pred:ctrl2_pred) %>%
#           lapply(pred_perform, mod_data, "deviance") %>%
#           perform_format("deviance") %>%
#           mutate( model = gsub("_pred","",model) )

# Expected Log Predictive Density
elpd <- cxval_pred %>% 
          dplyr::select(gaus_elpd:yr_t1_elpd) %>%
          apply(2, sum) %>% 
          as.matrix %>% 
          as.data.frame %>%
          tibble::add_column(model = rownames(.), .before=1) %>%
          mutate( model = gsub("_elpd","",model) ) %>%
          setNames( c("model", "elpd") )

# measures of fit
# mof  <- merge(mse, devi)
mof  <- merge(mse, elpd)

# store results ---------------------------------------------------------------------------
mod_summs <- Reduce(function(...) merge(...), 
                    list(mod_pars_diag, loo_df, waic_df, mof) ) %>%
                    arrange( mse )

write.csv(mod_summs,  paste0(args[3],"_mod_summaries_",      args[4],"_",spp_name,".csv"), row.names = F)
write.csv(posteriors, paste0(args[3],"_posterior_",          args[4],"_",spp_name,".csv"), row.names = F)
write.csv(cxval_pred, paste0(args[3],"_crossval_pred_diag_", args[4],"_",spp_name,".csv"), row.names = F)
