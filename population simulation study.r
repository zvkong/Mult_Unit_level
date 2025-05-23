# Load required library
source('./model.r')
library(dplyr); library(ggplot2);library(readr); library(coda)
library(sampling);library(tidyr)
library(HDInterval)
library(survey);library(ltm); library(mase); library(srvyr)
library(dplyr)

set.seed(2025)

# 1. Define parameters
J      <- 100             # Number of regions
mu     <- log(1000)       # Log-normal mean on log scale
sigma  <- 0.4             # Log-normal standard deviation

# Covariate probabilities
race_levels  <- c("White", "Hispanic", "Black", "Asian")
race_p       <- c(0.70, 0.15, 0.10, 0.05)

gender_levels <- c("Male", "Female")
gender_p      <- c(0.51, 0.49)

edu_levels   <- c("PhD", "Master", "Bachelor", "HighSchool", "Others")
edu_p        <- c(0.05, 0.10, 0.20, 0.55, 0.10)

# Income model coefficients
beta0      <- 50000    # Intercept
race_eff   <- c(White = 0, Hispanic = -5000, Black = -7000, Asian = 3000)
gender_eff <- c(Male = 0, Female = -2000)
edu_eff    <- c(PhD = 20000, Master = 10000, Bachelor = 5000, HighSchool = 0, Others = -2000)

sigma_u    <- 5000     # Standard deviation of region random effect
sigma_eps  <- 10000    # Standard deviation of individual error

# Total sample size for biased sampling
n <- 1000

# 2. Generate region population sizes N_j from a log-normal distribution
Nj <- round(exp(rnorm(J, mean = mu, sd = sigma)))

# 3. Build the full population data frame
pop_list <- lapply(seq_len(J), function(j) {
  n_j <- Nj[j]
  data.frame(
    region = j,
    race   = sample(race_levels,  size = n_j, replace = TRUE, prob = race_p),
    gender = sample(gender_levels, size = n_j, replace = TRUE, prob = gender_p),
    edu    = sample(edu_levels,    size = n_j, replace = TRUE, prob = edu_p)
  )
})
pop <- bind_rows(pop_list)

# 4. Add region random effects and generate income
u_vals <- rnorm(J, mean = 0, sd = sigma_u)
pop$u <- u_vals[pop$region]
pop$race_eff   <- race_eff[pop$race]
pop$gender_eff <- gender_eff[pop$gender]
pop$edu_eff    <- edu_eff[pop$edu]
pop$eps        <- rnorm(nrow(pop), mean = 0, sd = sigma_eps)

pop$income <- beta0 +
              pop$race_eff +
              pop$gender_eff +
              pop$edu_eff +
              pop$u +
              pop$eps

# 4.5 Add poverty status indicator: income below the 20th percentile is 1, else 0
poverty_threshold <- quantile(pop$income, 0.2)
pop$poverty <- ifelse(pop$income < poverty_threshold, 1, 0)
pop <- pop %>%
  mutate(
    poverty = ifelse(
      income < poverty_threshold,
      rbinom(n(), 1, prob = 0.85),  
      rbinom(n(), 1, prob = 0.05)   
    )
  )

# 5. Define sampling scores for biased (unequal-probability) sampling
pop$s <- with(pop,
  exp(
    0.2 +
    # 2.0 * (race == "White")    +  # Oversample White
    # 2.0 * (edu  == "Bachelor") +  # Oversample Bachelor
    # 0.5 * (gender == "Female") +  # Slightly oversample Female
    # 3 * (poverty == 1)  #     -  # Oversample Poor
    0.00008 * income                # Higher-income get higher score
  )
)
pop$income <- (pop$income - min(pop$income))/
            (max(pop$income) - min(pop$income))
# 6. Convert scores into inclusion probabilities π_i
# Total sample size for biased sampling
n <- 1000
pop$pi <- pop$s / sum(pop$s) * n
pop$pi <- pmin(pop$pi, 0.999)  # Truncate to avoid π_i ≥ 1

# 7. Draw an unequal-probability, no-replacement sample
set.seed(2025)
sample_idx <- sample(
  seq_len(nrow(pop)),
  size    = n,
  replace = FALSE,
  prob    = pop$pi
)
samp <- pop[sample_idx, ]
samp$w <- 1 / samp$pi  # Design weight = 1/π_i
samp$scaleWGT <- samp$w * nrow(samp) / sum(samp$w)  # Scaled weight
cor(samp$poverty, samp$income)
# 8. Compute true, unweighted, and HT-weighted means by region
pop_means <- pop %>%
  group_by(region) %>%
  summarise(true_mean = mean(income),
            true_mean_poverty = mean(poverty))

unweighted_means <- samp %>%
  group_by(region) %>%
  summarise(unweighted_mean = mean(income),
            unweighted_mean_poverty = mean(poverty))

weighted_means <- samp %>%
  group_by(region) %>%
  summarise(weighted_mean = sum(w * income) / sum(w),
            weighted_mean_poverty = sum(w * poverty) / sum(w))

# 9. Combine results and calculate biases
compare_df <- pop_means %>%
  left_join(unweighted_means, by = "region") %>%
  left_join(weighted_means,   by = "region") %>%
  mutate(
    bias_unweighted = unweighted_mean - true_mean,
    bias_weighted   = weighted_mean   - true_mean
  )

# 10. Display the comparison
print(compare_df)

# 11. Summarize biases across regions
summary(compare_df$bias_unweighted)
summary(compare_df$bias_weighted)

plot(compare_df$true_mean, compare_df$unweighted_mean, 
     xlab = "True Mean", ylab = "Unweighted Mean",
     main = "Unweighted vs True Means")
abline(0,1)

plot(compare_df$true_mean, compare_df$weighted_mean, 
     xlab = "True Mean", ylab = "HT-Weighted Mean",
     main = "HT-Weighted vs True Means")
abline(0,1)

plot(compare_df$true_mean_poverty, compare_df$unweighted_mean_poverty, 
     xlab = "True Mean Poverty", ylab = "Unweighted Mean Poverty",
     main = "Unweighted vs True Poverty Means")
abline(0,1)
plot(compare_df$true_mean_poverty, compare_df$weighted_mean_poverty, 
     xlab = "True Mean Poverty", ylab = "HT-Weighted Mean Poverty",
     main = "HT-Weighted vs True Poverty Means")
abline(0,1)

pcells <- pop %>% group_by(region, gender, race) %>% summarise(popsize=n())
pcells_region <- pop %>% group_by(region) %>% summarise(total_count = n())
predX <- model.matrix(~ gender + race - 1, data=pcells)
predPsi <- model.matrix(~as.factor(pcells$region)-1)

n_sim <- 10
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
cor_set <- c()


for (k in 1:n_sim){
  ### Subsample population 
    set.seed(k)
    sample_idx <- sample(
    seq_len(nrow(pop)),
    size    = n,
    replace = FALSE,
    prob    = pop$pi
    )
    samp <- pop[sample_idx, ]
    samp$w <- 1 / samp$pi  # Design weight = 1/π_i
    samp$scaleWGT <- samp$w * nrow(samp) / sum(samp$w)  # Scaled weight

    unweighted_means <- samp %>%
    group_by(region) %>%
    summarise(unweighted_mean = mean(income),
                unweighted_mean_poverty = mean(poverty))

    weighted_means <- samp %>%
    group_by(region) %>%
    summarise(weighted_mean = sum(w * income) / sum(w),
                weighted_mean_poverty = sum(w * poverty) / sum(w))

    # 9. Combine results and calculate biases
    compare_df <- pop_means %>%
    left_join(unweighted_means, by = "region") %>%
    left_join(weighted_means,   by = "region") %>%
    mutate(
        bias_unweighted = unweighted_mean - true_mean,
        bias_weighted   = weighted_mean   - true_mean
    )
  cor_set[k] <- cor(samp$poverty, samp$income)

  ## Fit model
    modwgt <- samp$scaleWGT
    modX <- model.matrix(~ gender + race - 1, data=samp)
    modPsi <- model.matrix(~as.factor(samp$region)-1)
    modPsi <- model.matrix(~ factor(region, 
                        levels=levels(as.factor(pcells$region)))-1,
                        data = samp)

    modY <- samp$income
    modZ <- samp$poverty
  
    ### Unit Gaussian Model
    unis_wage <- unis_gaus(X = modX, Y = modY, S = modPsi, 
                            sig2b=100, wgt = modwgt, n = NULL, predX = predX, predS = predPsi, 
                            nburn = nburn, nsim = nsim, nthin = 1, 
                            a = 0.5, b = 0.5, a_eps = 0.1, b_eps = 0.1)

    ### Unit Binomial Model
    unis_pov <- unis_bios(X = modX, Y = modZ, S = modPsi, 
                            sig2b=100, wgt = modwgt, n = NULL, predX = predX, predS = predPsi, 
                            nburn = nburn, nsim = nsim, nthin = 1, 
                            a = 0.1, b = 0.1)

    ### Multi Model wt binomial random
    mult_wage <- MTSM_br(X_1 = modX, X_2 = modX, Z_1 = modY, Z_2 = modZ, S = modPsi, 
                        sig2b = 100, wgt = modwgt, n = NULL, predX = predX, predS = predPsi, 
                        nburn = nburn, nsim = nsim, nthin = 1, 
                        sig2t = 1, sig2e = 1, tau_1_init = 1, tau_2_init = 1,  
                        a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1)

  ## Estimation for Gaussian Response
    results_ug <- gaus_post(unis_wage$Preds, 
                            unis_wage$sig2.chain, 
                            compare_df$true_mean,
                            region = pcells$region,
                            popsize = pcells$popsize)
    results_mg <- gaus_post(mult_wage$preds_gaus.chain, 
                            mult_wage$sig2.chain,
                            compare_df$true_mean,
                            region = pcells$region,
                            popsize = pcells$popsize)

    ## Binomial Response
    results_ub <- bios_post(unis_pov$Preds, compare_df$true_mean_poverty,
                            region = pcells$region,
                            popsize = pcells$popsize) 
    results_mb <- bios_post(mult_wage$preds_bios.chain, 
                            compare_df$true_mean_poverty,
                            region = pcells$region,
                            popsize = pcells$popsize)

## Save the prediction and quantile 0.025, 0.975
  ubios_pre[,k] <- results_ub$est
  mbios_pre[,k] <- results_mb$est
  dbio_pre[,k] <- compare_df$weighted_mean_poverty
  ugaus_pre[,k] <- results_ug$est
  mgaus_pre[,k] <- results_mg$est
  dgaus_pre[,k] <- compare_df$weighted_mean
  ubios_qual[,1,k] <- results_ub$lb; ubios_qual[,2,k] <- results_ub$ub
  mbios_qual[,1,k] <- results_mb$lb; mbios_qual[,2,k] <- results_mb$ub
  ugaus_qual[,1,k] <- results_ug$lb; ugaus_qual[,2,k] <- results_ug$ub
  mgaus_qual[,1,k] <- results_mg$lb; mgaus_qual[,2,k] <- results_mg$ub

  cat("Finished", k, "simulation dataset\n")
}

cr_ub <- numeric(n_sim); cr_mb <- numeric(n_sim)
cr_ug <- numeric(n_sim); cr_mg <- numeric(n_sim)
cr_dg <- numeric(n_sim); cr_db <- numeric(n_sim)
is_ub <- numeric(n_sim); is_mb <- numeric(n_sim)
is_ug <- numeric(n_sim); is_mg <- numeric(n_sim)
is_dg <- numeric(n_sim); is_db <- numeric(n_sim)
mse_ub <- numeric(n_sim); mse_mb <- numeric(n_sim)
mse_ug <- numeric(n_sim); mse_mg <- numeric(n_sim)
mse_dg <- numeric(n_sim); mse_db <- numeric(n_sim)
for (i in 1: n_sim){
  cr_ub[i] <- mean((ubios_qual[,1,i] < compare_df$true_mean_poverty) & 
  (compare_df$true_mean_poverty < ubios_qual[,2,i]))
  cr_mb[i] <- mean((mbios_qual[,1,i] < compare_df$true_mean_poverty) & 
  (compare_df$true_mean_poverty < mbios_qual[,2,i]))
  cr_ug[i] <- mean((ugaus_qual[,1,i] < compare_df$true_mean) 
              & (compare_df$true_mean < ugaus_qual[,2,i]))
  cr_mg[i] <- mean((mgaus_qual[,1,i] < compare_df$true_mean) 
              & (compare_df$true_mean < mgaus_qual[,2,i]))
#   cr_dg[i] <- mean((dgaus_qual[,1,i] < pop_means$true_mean) 
#               & (pop_means$true_mean < dgaus_qual[,2,i]))
#   cr_db[i] <- mean((dbio_qual[,1,i] < pop_means$true_mean_poverty) 
#               & (pop_means$true_mean_poverty < dbio_qual[,2,i]))
  is_ub[i]  <- interval_score(ubios_qual[,1,i], 
                            ubios_qual[,2,i], 
                            compare_df$true_mean_poverty)
  is_mb[i]  <- interval_score(mbios_qual[,1,i], 
                            mbios_qual[,2,i], 
                            compare_df$true_mean_poverty)
  is_ug[i]  <- interval_score(ugaus_qual[,1,i], 
                              ugaus_qual[,2,i], 
                              compare_df$true_mean)
  is_mg[i]  <- interval_score(mgaus_qual[,1,i], 
                              mgaus_qual[,2,i], 
                              compare_df$true_mean)
#   is_dg[i]  <- interval_score(dgaus_qual[,1,i], 
#                               dgaus_qual[,2,i], 
#                               pop_means$true_mean)
#   is_db[i]  <- interval_score(dbio_qual[,1,i], 
#                                 dbio_qual[,2,i], 
#                                 pop_means$true_mean_poverty)
  mse_mg[i] <- mean((mgaus_pre[,i] - compare_df$true_mean)^2,na.rm = TRUE)
  mse_ug[i] <- mean((ugaus_pre[,i] - compare_df$true_mean)^2,na.rm = TRUE)
  mse_ub[i] <- mean((ubios_pre[,i] - compare_df$true_mean_poverty)^2,
                    na.rm = TRUE)
  mse_mb[i] <- mean((mbios_pre[,i] - compare_df$true_mean_poverty)^2,
                    na.rm = TRUE)
  mse_dg[i] <- mean((dgaus_pre[,i] - compare_df$true_mean)^2,na.rm = TRUE)
  mse_db[i] <- mean((dbio_pre[,i] - compare_df$true_mean_poverty)^2,
                    na.rm = TRUE)
}

mse_mg/mse_ug
mse_ug/mse_dg
mse_mb/mse_ub
mse_ub/mse_db
cor_set
 # save.image("./data/simpop100_SET2.rdata")


boxplot(mse_mg/mse_ug, mse_mb/mse_ub, 
        names = c("Gaussian", "Binomial"), 
        main = "MSE Ratio: Gaussian vs Binomial")
plot(mult_wage$Beta_2.chain[1,], type = 'l')
plot(mult_wage$Beta_2.chain[2,], type = 'l')
