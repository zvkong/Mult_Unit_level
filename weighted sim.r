# Load required library
source('./code/wt share ran.r')
library(dplyr); library(ggplot2);library(readr); library(coda)
library(sampling);library(tidyr)
library(HDInterval)
library(survey);library(ltm); library(mase); library(srvyr)
library(dplyr)

pums21 <- read_csv('data/IL21.csv')

pums1 <- pums21 %>% dplyr::select(PUMA, PWGTP, SEX, 
                    RAC1P, POVPIP, PINCP, SCHL) %>% 
  remove_missing() %>%
  mutate(PWGTP=as.numeric(PWGTP), POVPIP=as.numeric(POVPIP), SEX=factor(SEX), 
        RACE=factor(RAC1P), LOGINCO=log(as.numeric(PINCP)), 
        INCOME = as.numeric(PINCP), degree = as.integer(SCHL),
        income = sqrt(PINCP)) %>% 
  mutate(POV=case_when(
    POVPIP <= 100 ~ 1,
    POVPIP > 100 ~ 0, 
  ),
  BACH = case_when(SCHL >= 21 ~ 1,
    SCHL < 21 ~ 0,
  ) %>% factor())


## remove INCO is zero (Or should add a small value)
pums1 <- pums1[pums1$LOGINCO != -Inf,]

pums <- pums1 %>% dplyr::select(PUMA, PWGTP, SEX, RACE, POV, LOGINCO, INCOME, BACH) %>% 
  remove_missing() %>%
  mutate(POV = as.numeric(POV), EDU = as.numeric(BACH)-1)
# pums$INCO <- (pums$LOGINCO - min(pums$LOGINCO))/(max(pums$LOGINCO) - min(pums$LOGINCO))
# pums$INCO <- pums$LOGINCO
# abline(0,1)
hist(pums$LOGINCO)
bc_result <- boxcox(INCOME ~ 1, data = pums,
       lambda = seq(-2, 2, by = 0.01))
lambda <- bc_result$x[which.max(bc_result$y)]
bc_income <- ((pums$INCOME^lambda - 1)/lambda)
pums$INCO <- (bc_income - min(bc_income))/(max(bc_income) - min(bc_income))

## Compute the true mean INCO and POV rate for each area
truth <- pums %>% group_by(PUMA) %>% 
          summarise(INCO = mean(INCO),
                    POV = mean(POV))


## Choose bachlor degree and gender as the covariates
pcells <- pums %>% group_by(PUMA, SEX, BACH) %>% summarise(popsize=n())
pcells_region <- pums %>% group_by(PUMA) %>% summarise(popsize=n())
predX <- model.matrix(~ SEX + BACH, data=pcells)
predPsi <- model.matrix(~as.factor(pcells$PUMA)-1)


n_sim <- 5
nsim <- 1000
nburn <- 2000
nthin <- 1

## start 100 simulation data sets
ubios_pre <- array(NA, dim = c(nrow(pcells_region), n_sim))
mbios_pre <- array(NA, dim = c(nrow(pcells_region), n_sim))
ugaus_pre <- array(NA, dim = c(nrow(pcells_region), n_sim))
mgaus_pre <- array(NA, dim = c(nrow(pcells_region), n_sim))
dgaus_pre <- array(NA, dim = c(nrow(pcells_region), n_sim))
dbio_pre <- array(NA, dim = c(nrow(pcells_region), n_sim))
ubios_qual <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
mbios_qual <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
ugaus_qual <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
mgaus_qual <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
dgaus_qual <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
dbio_qual <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
mgaus_pre_2 <- array(NA, dim = c(nrow(pcells_region), n_sim))
mbios_pre_2 <- array(NA, dim = c(nrow(pcells_region), n_sim))
mgaus_qual_2 <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
mbios_qual_2 <- array(NA, dim = c(nrow(pcells_region), 2, n_sim))
cor_set <- c()


for (k in 1:n_sim){
    set.seed(5+k)
    ss <- 2000

    prob <- inclusionprobabilities(pums$PWGTP*(1 + (pums$POV == '1')), ss)
    # prob <- inclusionprobabilities(pums$PWGTP, ss)
    # prob <- prob*ss/sum(prob)
    Ind <- UPsystematic(prob)
    samp <- pums[as.logical(Ind),]
    samp$P <- prob[as.logical(Ind)]
    sum(samp$P)
    samp$W <- 1/ samp$P
    # samp$W <- samp$PWGTP
    samp$scaledWGT <- (samp$W)*nrow(samp)/sum(samp$W)


compare_df <- samp %>% group_by(PUMA) %>% 
  summarise(unweighted_means = mean(INCO),
            unweighted_means_POV = mean(POV),
               unweighted_means_EDU = mean(EDU),
            weighted_means = weighted.mean(INCO, W),
            weighted_means_POV = weighted.mean(POV, W),
            weighted_means_EDU = weighted.mean(EDU, W))
    cor_set[k] <- cor(samp$POV, samp$INCO)

    ## Fit model
    modwgt <- samp$scaledWGT
    modX <- model.matrix(~ SEX + BACH, data=samp)
    modPsi <- model.matrix(~as.factor(samp$PUMA)-1)
    modY <- samp$INCO
    modZ <- samp$POV

    unis_wage <- unis_gaus(X = modX, Y = modY, S = modPsi, 
                            sig2b=1000, wgt = modwgt, n = NULL, 
                            predX = predX, predS = predPsi, 
                            nburn = nburn, nsim = nsim, nthin = 1, 
                            a = 0.5, b = 0.5, a_eps = 0.1, b_eps = 0.1)
    ### Unit Binomial Model
    unis_pov <- unis_bios(X = modX, Y = modZ, S = modPsi, 
                            sig2b=1000, wgt = modwgt, n = NULL, 
                            predX = predX, predS = predPsi, 
                            nburn = nburn, nsim = nsim, nthin = 1, 
                            a = 0.1, b = 0.1)

    ### Multi Model wt random
    mult_wage <- MTSM_gr(X_1 = modX, X_2 = modX, Z_1 = modY, Z_2 = modZ, S = modPsi,
                        sig2b = 1000, wgt = modwgt, n = NULL, 
                        predX = predX, predS = predPsi, 
                        nburn = nburn, nsim = nsim, nthin = 1, 
                        sig2t = 10, sig2e = 10, tau_1_init = 1, tau_2_init = -1,
                        a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1)

    ### Multi Model wt random Binomial
    mult_pov <- MTSM_br(X_1 = modX, X_2 = modX, Z_1 = modY, Z_2 = modZ, S = modPsi,
                        sig2b = 1000, wgt = modwgt, n = NULL, 
                        predX = predX, predS = predPsi, 
                        nburn = nburn, nsim = nsim, nthin = 1, 
                        sig2t = 10, sig2e = 10, tau_1_init = 1, tau_2_init = -1,
                        a_eps = 0.1, b_eps = 0.1)
  ## Estimation for Gaussian Response
    results_ug <- gaus_post(unis_wage$Preds, 
                            unis_wage$sig2.chain, 
                            truth$INCO,
                            region = pcells$PUMA,
                            popsize = pcells$popsize)
    results_mg <- gaus_post(mult_wage$preds_gaus.chain, 
                            mult_wage$sig2.chain,
                            truth$INCO,
                            region = pcells$PUMA,
                            popsize = pcells$popsize)
    results_mg_2 <- gaus_post(mult_pov$preds_gaus.chain, 
                            mult_pov$sig2.chain,
                            truth$INCO,
                            region = pcells$PUMA,
                            popsize = pcells$popsize)

    ## Binomial Response
    results_ub <- bios_post(unis_pov$Preds, truth$POV,
                        region = pcells$PUMA,
                        popsize = pcells$popsize)
    results_mb <- bios_post(mult_wage$preds_bios.chain, 
                        truth$POV,
                        region = pcells$PUMA,
                        popsize = pcells$popsize)
    results_mb_2 <- bios_post(mult_pov$preds_bios.chain,
                        truth$POV,
                        region = pcells$PUMA,
                        popsize = pcells$popsize)

## Save the prediction and quantile 0.025, 0.975
  ubios_pre[,k] <- results_ub$est
  mbios_pre[,k] <- results_mb$est
  dbio_pre[,k] <- compare_df$weighted_means_POV
  ugaus_pre[,k] <- results_ug$est
  mgaus_pre[,k] <- results_mg$est
    mgaus_pre_2[,k] <- results_mg_2$est
    mbios_pre_2[,k] <- results_mb_2$est
  dgaus_pre[,k] <- compare_df$weighted_means

  ubios_qual[,1,k] <- results_ub$lb; ubios_qual[,2,k] <- results_ub$ub
  mbios_qual[,1,k] <- results_mb$lb; mbios_qual[,2,k] <- results_mb$ub
  ugaus_qual[,1,k] <- results_ug$lb; ugaus_qual[,2,k] <- results_ug$ub
  mgaus_qual[,1,k] <- results_mg$lb; mgaus_qual[,2,k] <- results_mg$ub
  mgaus_qual_2[,1,k] <- results_mg_2$lb; mgaus_qual_2[,2,k] <- results_mg_2$ub
  mbios_qual_2[,1,k] <- results_mb$lb; mbios_qual_2[,2,k] <- results_mb$ub

  cat("Finished", k, "simulation dataset\n")
}

cr_ub <- numeric(n_sim); cr_mb <- numeric(n_sim)
cr_ug <- numeric(n_sim); cr_mg <- numeric(n_sim)
cr_dg <- numeric(n_sim); cr_db <- numeric(n_sim)
is_ub <- numeric(n_sim); is_mb <- numeric(n_sim)
is_ug <- numeric(n_sim); is_mg <- numeric(n_sim)
is_dg <- numeric(n_sim); is_db <- numeric(n_sim)
is_mb_2 <- numeric(n_sim); is_mg_2 <- numeric(n_sim)
mse_ub <- numeric(n_sim); mse_mb <- numeric(n_sim)
mse_ug <- numeric(n_sim); mse_mg <- numeric(n_sim)
mse_dg <- numeric(n_sim); mse_db <- numeric(n_sim)
mse_mg_2 <- numeric(n_sim); mse_mb_2 <- numeric(n_sim)
for (i in 1: n_sim){
  cr_ub[i] <- mean((ubios_qual[,1,i] < compare_df$weighted_mean_POV) & 
  (compare_df$weighted_mean_POV < ubios_qual[,2,i]))
  cr_mb[i] <- mean((mbios_qual[,1,i] < compare_df$weighted_mean_POV) & 
  (compare_df$weighted_mean_POV < mbios_qual[,2,i]))
  cr_ug[i] <- mean((ugaus_qual[,1,i] < compare_df$weighted_mean) 
              & (compare_df$weighted_mean < ugaus_qual[,2,i]))
  cr_mg[i] <- mean((mgaus_qual[,1,i] < compare_df$weighted_mean) 
              & (compare_df$weighted_mean < mgaus_qual[,2,i]))
#   cr_dg[i] <- mean((dgaus_qual[,1,i] < pop_means$weighted_mean) 
#               & (pop_means$weighted_mean < dgaus_qual[,2,i]))
#   cr_db[i] <- mean((dbio_qual[,1,i] < pop_means$weighted_mean_POV) 
#               & (pop_means$weighted_mean_POV < dbio_qual[,2,i]))
  is_ub[i]  <- interval_score(ubios_qual[,1,i], 
                            ubios_qual[,2,i], 
                            truth$POV)
  is_mb[i]  <- interval_score(mbios_qual[,1,i], 
                            mbios_qual[,2,i], 
                            truth$POV)
    is_mb_2[i]  <- interval_score(mbios_qual_2[,1,i], 
                            mbios_qual_2[,2,i], 
                            truth$POV)
  is_ug[i]  <- interval_score(ugaus_qual[,1,i], 
                              ugaus_qual[,2,i], 
                              truth$INCO)
  is_mg[i]  <- interval_score(mgaus_qual[,1,i], 
                              mgaus_qual[,2,i], 
                              truth$INCO)
    is_mg_2[i]  <- interval_score(mgaus_qual_2[,1,i], 
                              mgaus_qual_2[,2,i], 
                              truth$INCO)    
#   is_dg[i]  <- interval_score(dgaus_qual[,1,i], 
#                               dgaus_qual[,2,i], 
#                               pop_means$weighted_mean)
#   is_db[i]  <- interval_score(dbio_qual[,1,i], 
#                                 dbio_qual[,2,i], 
#                                 pop_means$weighted_mean_POV)
  mse_mg[i] <- mean((mgaus_pre[,i] - truth$INCO)^2,na.rm = TRUE)
  mse_ug[i] <- mean((ugaus_pre[,i] - truth$INCO)^2,na.rm = TRUE)
  mse_ub[i] <- mean((ubios_pre[,i] - truth$POV)^2,
                    na.rm = TRUE)
  mse_mb[i] <- mean((mbios_pre[,i] - truth$POV)^2,
                    na.rm = TRUE)
  mse_dg[i] <- mean((dgaus_pre[,i] - truth$INCO)^2,na.rm = TRUE)
  mse_db[i] <- mean((dbio_pre[,i] - truth$POV)^2,
                    na.rm = TRUE)
    mse_mg_2[i] <- mean((mgaus_pre_2[,i] - truth$INCO)^2,na.rm = TRUE)
    mse_mb_2[i] <- mean((mbios_pre_2[,i] - truth$POV)^2,
                        na.rm = TRUE)
}

mse_mg/mse_ug
mse_ug/mse_dg
mean(mse_mb/mse_ub)
mse_ub/mse_db
mse_mg_2/mse_ug
mean(mse_mb_2/mse_ub)
is_mg/is_ug
is_mg_2/is_ug
is_mb/is_ub
is_mb_2/is_ub
cor_set
# save.image("./data/simpop100_SET2.rdata")
