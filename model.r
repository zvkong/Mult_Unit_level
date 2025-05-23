library(Matrix)
library(MASS)
library(BayesLogit)
## Build functions

## Sum of subvectors
sms <- function(v, n){
# Initialize an empty vector to store the sums of the subvectors
  subvector_sums <- numeric(length(n))
  
  # Starting index for the first subvector
  start_index <- 1
  
  # Loop through each size in n to create subvectors and compute their sums
  for (i in seq_along(n)) {
    end_index <- start_index + n[i] - 1  # Calculate the ending index
    subvector_sums[i] <- sum(v[start_index:end_index])  # Compute sum and store
    start_index <- end_index + 1  # Update the start index for the next subvector
  }
  
  return(subvector_sums)
}
## Create design matrix
cdm <- function(s, n) {
  # s is the number of different Ji vectors
  # n is the vector of repeat counts for each Ji
  
  # Check if the length of n matches s
  if (length(n) != s) {
    stop("Length of n must be equal to s.")
  }
  
  # Initialize an empty matrix to store the final results
  total_rows <- sum(n)  # Total number of rows in the final matrix
  design_matrix <- matrix(0, nrow = total_rows, ncol = s)  # Start with all zeros
  
  # Create each Ji vector and repeat it ni times
  start_index <- 1
  for (i in 1:s) {
    end_index <- start_index + n[i] - 1
    Ji <- rep(0, s)
    Ji[i] <- 1  # Set the i-th position to 1
    design_matrix[start_index:end_index, ] <- matrix(rep(Ji, n[i]), nrow = n[i], byrow = TRUE)
    start_index <- end_index + 1
  }
  
  return(design_matrix)
}
## MSE
MSE <- function(x, y){
  return(mean((x - y)^2))
}
## Inverse logit function
inv.logit <- function(x){
  return(exp(x)/(1 + exp(x)))
}
## Draw from MVN with sparse variance matrix function
rmvn <- function(n, mu, Sigma) {
  p <- dim(Sigma)[1]
  L <- Matrix::chol(Sigma)
  z <- Matrix(rnorm(n * p), p, n)
  #z <- rnorm(p)
  result <- (t(mu +  L %*% z))
  return(result)
}

## Quantile function
QLU <- function(x){
  q <- apply(x, 1, quantile, probs = c(0.025, 0.975))
  return(q)
}

## Coverage rate
cr <- function(quantiles, true){
    rate <- c()
    for (i in 1:dim(quantiles)[1]){
    rate[i] <- mean(quantiles[i,1,] <= true[i] & quantiles[i,2,] >= true[i])
    }   
    return(rate)
}

# # Interval score function
interval_score <- function(lower, upper, truth, alpha = 0.05){
  width <- upper - lower
  penalty_low  <- (2/alpha) * pmax(0, lower - truth)
  penalty_high <- (2/alpha) * pmax(0, truth - upper)
  width + penalty_low + penalty_high
}

## Poststratification
gaus_post_mean <- function(preds,
                       true_mean,
                       region,
                       popsize) {  
  regions <- unique(region)
  R       <- ncol(preds)
  C       <- length(regions)
  
  idx_by_reg <- split(seq_along(region), region)
  N_by_reg   <- sapply(idx_by_reg, function(id) sum(popsize[id]))
  
  post <- matrix(NA, nrow = C, ncol = R, dimnames = list(regions, NULL))
  
  for (r in seq_len(R)) {
    mu_j <- preds[, r]
    for (c in seq_along(regions)) {
      ids         <- idx_by_reg[[ c ]]
      post[c, r]  <- sum(popsize[ids] * mu_j[ids]) / N_by_reg[c]
    }
  }
  
  est <- rowMeans(post)  
  lb  <- apply(post, 1, quantile, probs = 0.025)
  ub  <- apply(post, 1, quantile, probs = 0.975)
  cr  <- mean(lb <= true_mean & true_mean <= ub)
  
  list(est = est, lb = lb, ub = ub, cr = cr)
}
gaus_post <- function(preds, sig2chain, true_mean, region, popsize){
  df <- array(NA, dim=c(length(unique(region)), nsim))
  temp <- data.frame(region=region, mu=popsize * preds)
  temp2 <- data.frame(region=region, 
                  popsize = popsize, temp[,-1]) %>% 
                  group_by(region) %>% summarise_all(sum)
  sig2 <- sig2chain
  for(i in 1:length(unique(region))){
  temp3 <- rnorm(nsim, t(temp2[i,-(1:2)]/temp2$popsize[i]), 
          sqrt(sig2/temp2$popsize[i]))
  df[i,] <- temp3
  }
  est <- rowMeans(df)
  lb <- apply(df, 1, quantile, probs = c(0.025))
  ub <- apply(df, 1, quantile, probs = c(0.975))
  cr <- mean(lb <= true_mean & ub >= true_mean)
  list(est = est, lb = lb, ub = ub, cr = cr)
}

bios_post <- function(preds, true_mean, region, popsize){
  df <- array(NA, dim=c(length(unique(region)), dim(preds)[2]))
  for(j in 1:dim(preds)[2]){
      temp <- data.frame(region=region,
                         N=popsize, 
                         P=I(preds[,j]))
      temp2 <- (apply(cbind(temp$N, temp$P), 1, 
                function(x) rbinom(1, x[1], x[-1])))
      temp2 <- data.frame(region=region, 
                          popsize = popsize, 
                          temp2) %>% 
              group_by(region) %>% summarize_all(sum)
      temp3 <- as.matrix(temp2[,3]/temp2[,2], ncol = 1)
      df[,j] <- temp3
  }

  est <- rowMeans(df)
  lb <- apply(df, 1, quantile, probs = c(0.025))
  ub <- apply(df, 1, quantile, probs = c(0.975))
  cr <- mean(lb <= true_mean & ub >= true_mean)
  list(est = est, lb = lb, ub = ub, cr = cr)
}
## let w represent for weight and n represent for populations
## UNIS model for Binomial response
unis_bios <- function(X, Y, S, sig2b=1000, wgt = NULL, n = NULL, 
                    predX, predS, nburn = 1000, nsim = 5000, nthin = 1, 
                    a = 0.1, b = 0.1) {
    ## X is the covariates matrix for binomial response Z
    ## S is the basis function for spatial information
    ## wgt is the weights
    ## sig2b is the prior variance for beta 
    ## n is the total number parameter for binomial response, should be 1 if this is binomial
    ## nburn is burning in 
    ## nsim is the simulation number
    N <- length(Y)
    p <- dim(X)[2]
    r <- ncol(S)
    npred <- dim(predX)[1]
    if(is.null(wgt)) wgt <- rep(1, N)
    if(is.null(n)) n <- rep(1, N)

    k <- wgt * (Y - n/2)
    Ip <- Diagonal(p) ## identity matrix of p dimention
    #IN <- Diagonal(N) ## identity matrix of N dimention
    Ir <- Diagonal(r) ## identity matrix of r dimention

    # Initialize parameters
    Mu <- rep(1, N)  # Initialize theta with zeros
    Beta <- rep(1, p)  # Initialize beta as a column vector
    Sigma2_u <- 1 # Initialize variance for u
    U <- rep(1, r) # Initialize the spatial random effect u

    # Storage for parameters
    Beta.chain <- array(0, dim = c(p, nsim / nthin))
    Sigma2_u.chain <- rep(0, nsim / nthin)
    U.chain <- array(0, dim = c(r, nsim / nthin))
    Mu.chain <- array(0, dim = c(N, nsim / nthin))
    preds.chain <- array(0, dim = c(npred, nsim / nthin))
    logit_bios.chain <- array(0, dim = c(npred, nsim / nthin))
    print(paste0("Starting ",nsim / nthin, " iterations."))
    pb <- txtProgressBar(min=0, max=(nsim + nburn) / nthin, style=3)

    for (index in 1:(nsim + nburn)) {
        if (index %% 10000 == 0) cat(index, "\n")
        # Update latent variable w
        omega <- rpg.gamma(N, wgt, Mu)
        Omega <- Diagonal(x=omega)
        SO <- t(S) %*% Omega
        OM <- SO %*% S 
        gamma <- k/omega
        
        # Update Beta
        var.Beta <- solve(t(X) %*% Omega %*% X +  Ip/sig2b)
        mean.Beta <- var.Beta %*% (t(X) %*% Omega %*% (gamma - S %*% U))
        Beta <- as.vector(mvrnorm(1, mu = mean.Beta, Sigma = var.Beta))
        
        # Update Sigma2_u
        a_prime <- a + r / 2
        b_prime <- b + 0.5 * t(U) %*% U
        Sigma2_u <- 1 / rgamma(1, shape = a_prime, rate = b_prime) 
        
        # Update U
        var.U <- solve(OM  + Ir/Sigma2_u)
        mean.U <- var.U %*% SO %*%(gamma - X %*% Beta)
        U <- as.vector(rmvn(1, as.vector(mean.U), var.U))
        # U <- U - mean(U) ## sum to zero of random effect


        # Update mu
        Mu <- as.vector(X %*% Beta + S %*% U)

        # Pred
        logit_bios <- predX %*% Beta + predS %*% U
        preds <- plogis(logit_bios)

        setTxtProgressBar(pb, index)

        if (index > nburn && (index - nburn) %% nthin == 0) {
            Beta.chain[,(index - nburn) / nthin] <- Beta
            U.chain[,(index - nburn) / nthin] <- U
            Sigma2_u.chain[(index - nburn) / nthin] <- Sigma2_u
            Mu.chain[,(index-nburn)/nthin] <- Mu
            preds.chain[,(index-nburn)/nthin] <- preds
            logit_bios.chain[,(index-nburn)/nthin] <- predX %*% Beta + predS %*% U
        }
    }

    list(Beta.chain = Beta.chain, U.chain = U.chain, 
            Sigma2_u.chain = Sigma2_u.chain,
            Mu.chain = Mu.chain, Preds = preds.chain,
            logit_bios.chain = logit_bios.chain)
}

## UNIS model for Gaussian response
unis_gaus <- function(X, Y, S,  sig2b=1000, wgt = NULL, n = NULL, 
                    predX, predS, nburn = 1000, nsim = 5000, nthin = 1, 
                    a = 0.1, b = 0.1, a_eps = 0.1, b_eps = 0.1) {
    ## X is the covariates matrix for binomial response Z
    ## S is the basis function for spatial information
    ## wgt is the weights
    ## sig2b is the prior variance for beta 
    ## n is the total number parameter for binomial response, should be 1 if this is binomial
    ## nburn is burning in 
    ## nsim is the simulation number
    N <- length(Y)
    p <- dim(X)[2]
    r <- ncol(S)
    npred <- dim(predX)[1]
    if(is.null(wgt)) wgt <- rep(1, N)
    if(is.null(n)) n <- rep(1, N)
    w <- sum(wgt)
    Wgt <- Diagonal(x=wgt)
    Ip <- Diagonal(p) ## identity matrix of p dimention
    Ir <- Diagonal(r) ## identity matrix of r dimention


    # Initial values
    Beta <- rep(1, p)
    Sigma2_u <- 1
    U <- rep(0, r)
    Mu <- rep(0, N)
    sig2 <- 1
    
    # Chain containers
    Beta.chain <- array(0, dim = c(p, nsim / nthin))
    Sigma2_u.chain <- rep(0, nsim / nthin)
    U.chain <- array(0, dim = c(r, nsim / nthin))
    Mu.chain <- array(0, dim = c(N, nsim / nthin))
    sig2.chain <- rep(0, nsim / nthin)
    preds.chain <- array(0, dim = c(npred, nsim / nthin))
    print(paste0("Starting ",nsim / nthin, " iterations."))
    pb <- txtProgressBar(min=0, max=(nsim + nburn) / nthin, style=3)

    for (index in 1:(nsim + nburn)) {
        if (index %% 10000 == 0) cat(index, "\n")

        # Update sig2
        a_star <- as.numeric(a_eps + w/2)
        b_star <- as.numeric(b_eps + 0.5 * t(Y-Mu)%*%Wgt%*%(Y-Mu))
        sig2 <- 1/rgamma(1, shape = a_star, rate = b_star)
        d <- Wgt/sig2
        SD <- t(S) %*% d
        M <- SD %*% S
        XD <- t(X) %*% d %*% X
        
        # Update Beta
        var.Beta <- solve(XD + Ip/sig2b)
        mean.Beta <- var.Beta %*% t(X) %*% d %*% (Y - S %*% U)
        Beta <- as.vector(mvrnorm(1, mu = mean.Beta, Sigma = var.Beta))
        
        # Update Sigma2_u
        a_prime <- a + r / 2
        b_prime <- b + 0.5 * t(U) %*% U
        Sigma2_u <- 1 / rgamma(1, shape = a_prime, rate = b_prime)
        
        # Update U
        var.U <- solve(M + Ir/Sigma2_u)
        mean.U <- var.U %*% (SD %*% (Y - X %*% Beta))
        U <- as.vector(rmvn(1, as.vector(mean.U), var.U))
        # U <- U - mean(U) ## sum to zero 

        # Update mu
        Mu <- X %*% Beta + S %*% U
        
        # Pred
        preds <- predX %*% Beta + predS %*% U

        setTxtProgressBar(pb, index)
        
        if (index > nburn && (index - nburn) %% nthin == 0) {
        Beta.chain[, (index - nburn) / nthin] <- Beta
        U.chain[, (index - nburn) / nthin] <- U
        Sigma2_u.chain[(index - nburn) / nthin] <- Sigma2_u
        Mu.chain[,(index-nburn)/nthin] <- Mu
        sig2.chain[(index - nburn) / nthin] <- sig2
        preds.chain[,(index-nburn)/nthin] <- preds
        }
    }

  list(Beta.chain = Beta.chain, U.chain = U.chain, 
        Sigma2_u.chain = Sigma2_u.chain, Mu.chain = Mu.chain, 
        sig2.chain = sig2.chain, Preds = preds.chain)
}

## Multitype spatial model
MTSM_br <- function(X_1, X_2, Z_1, Z_2, S, sig2b = 1000, wgt = NULL, n = NULL, 
                predX, predS, n_preds, nburn = 1000, nsim = 5000, nthin = 1, 
                sig2t = 10, sig2e = 10, tau_1_init = 1, tau_2_init = 1,  
                a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1, 
                alambda = 0.1, blambda = 0.1) {
    ## Z_1 is Gaussian response,
    ## Z_2 is binomial response with effective sample size m
    ## X is the covariates matrix for binomial response Z
    ## S is the basis function for spatial information
    ## wgt is the weights
    ## sig2b is the prior variance for beta 
    ## n is the total number parameter for binomial response, should be 1 if this is binomial
    ## nburn is burning in 
    ## nsim is the simulation number
    ## aeta, beta is hyperameter for sig2e
    ## alambda, blambda is hyperameter for sig2l
    N <- dim(X_1)[1] # Number of observations(or units)
    r <- ncol(S)
    npred <- dim(predX)[1]
    if(is.null(wgt)) wgt <- rep(1, N)
    if(is.null(n)) n <- rep(1, N)
    w <- sum(wgt)
    Wgt <- Diagonal(x=wgt)
    p_1 <- dim(X_1)[2] # length for Beta_1
    p_2 <- dim(X_2)[2] # length for Beta_2
    r <- dim(S)[2] # length for random effect
    Ip1 <- Diagonal(p_1) ## identity matrix of p dimention
    Ip2 <- Diagonal(p_2) ## identity matrix of p2 dimention
    Ir <- Diagonal(r) ## identity matrix of r dimention
    t_X_1 <- t(X_1)
    t_X_2 <- t(X_2)
    k <- wgt * (Z_2 - 1/2)

    # Initial values
    tau_1 <- tau_1_init
    tau_2 <- tau_2_init
    Beta_1 <- rep(1, p_1)
    Beta_2 <- rep(1, p_2)
    lambda <- rep(0, r)
    eta <- rep(1, r)
    Mu_1 <- rep(1, N)
    Mu_2 <- rep(1, N)
    sig2 <- 1
    sig2e <- sig2e
    sig2l <- 1


    # Chain containers  
    tau_1.chain <- array(0, dim = c(1, nsim / nthin))
    tau_2.chain <- array(0, dim = c(1, nsim / nthin))
    Beta_1.chain <- array(0, dim = c(p_1, nsim / nthin))
    Beta_2.chain <- array(0, dim = c(p_2, nsim / nthin))
    Sigma2_lambda.chain <- rep(0, nsim / nthin)
    Sigma2_eta.chain <- rep(0, nsim / nthin)
    sig2.chain <- rep(0, nsim / nthin)
    lambda.chain <- array(0, dim = c(r, nsim / nthin))
    eta.chain <- array(0, dim = c(r, nsim / nthin))
    Mu_1.chain <- array(0, dim = c(N, nsim / nthin))
    Mu_2.chain <- array(0, dim = c(N, nsim / nthin))
    preds_gaus.chain <- array(0, dim = c(npred, nsim / nthin))
    preds_bios.chain <- array(0, dim = c(npred, nsim / nthin))
    logit_bios.chain <- array(0, dim = c(npred, nsim / nthin))
    print(paste0("Starting ",nsim / nthin, " iterations."))
    pb <- txtProgressBar(min=0, max=(nsim + nburn) / nthin, style=3)
    
    for (index in 1:(nsim + nburn)) {
        if (index %% 10000 == 0) cat(index, "\n")
        # Update sig2
        a_star <- as.numeric(a_eps + w/2)
        b_star <- as.numeric(b_eps + 0.5 * t(Z_1-Mu_1)%*%Wgt%*%(Z_1-Mu_1))
        sig2 <- 1/rgamma(1, shape = a_star, rate = b_star)
        d <- Wgt/sig2
        SD <- t(S) %*% d
        M <- SD %*% S
        XD <- t_X_1 %*% d %*% X_1
     
        # Update latent variable w
        omega <- rpg.gamma(N, wgt, Mu_2)
        Omega <- Diagonal(x=omega)
        SO <- t(S) %*% Omega
        OM <- SO %*% S 
        gamma <- k/omega
        
        # Update sigeta
        a_eta <- aeta + r/2
        b_eta <- beta + 0.5 * t(eta) %*% eta
        sig2e <- 1 / rgamma(1, shape = a_eta, rate = b_eta)

        # Update eta
        var.eta <- solve(tau_1 * M * tau_1 +
                        Ir/sig2e +
                        tau_2 * OM * tau_2, sparse = TRUE)
        mean.eta <- var.eta %*%
        (tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1) +
            tau_2 * SO %*% (gamma - X_2 %*% Beta_2 - S %*% lambda))
        eta <- as.vector(rmvn(1, as.vector(mean.eta), var.eta))
        # eta <- eta - mean(eta) ## sum to zero

        # Update Beta_1
        var.Beta_1 <- solve(XD + Ip1/sig2b)
        mean.Beta_1 <- var.Beta_1 %*% t_X_1 %*% d %*% (Z_1 - tau_1 * S %*% eta)
        Beta_1 <- as.vector(mvrnorm(1, mean.Beta_1, var.Beta_1))

        # Update sig2l
        a_l <- alambda + r / 2
        b_l <- blambda + 0.5 * t(lambda) %*% lambda
        sig2l <- 1 / rgamma(1, shape = a_l, rate = b_l)

        # Update lambda
        var.lambda <- solve(OM + Ir/sig2l, sparse = TRUE)
        mean.lambda <- var.lambda %*% (SO %*%(gamma -
                                       X_2 %*% Beta_2 - tau_2 * S %*% eta))
        lambda <- as.vector(rmvn(1, as.vector(mean.lambda), var.lambda))
        # lambda <- lambda - mean(lambda) ## sum to zero

        # Update Beta_2
        var.Beta_2 <- solve(t_X_2 %*% Omega %*% X_2 + Ip2/sig2b)
        mean.Beta_2 <- var.Beta_2 %*% t_X_2 %*% Omega %*% (gamma - 
                                      tau_2 * S %*% eta - S %*% lambda)
        Beta_2 <- as.vector(mvrnorm(1, mean.Beta_2, var.Beta_2))
        
        # Update regression parameter tau_1
        # var_tau_1 <- as.numeric(solve(t(eta) %*% M %*% eta +
        #                                 1/sig2t))
        # mean_tau_1 <- as.numeric(var_tau_1 * t(eta) %*% SD %*%
        #                         (Z_1 - X_1 %*% Beta_1))
        # tau_1 <- as.numeric(rnorm(1, mean_tau_1, sqrt(var_tau_1)))
        
        # Update regression parameter tau_2
        var_tau_2 <- as.numeric(solve(t(eta) %*% OM %*% eta +
                                        1/sig2t))
        mean_tau_2 <- as.numeric(var_tau_2 * t(eta) %*% SO %*%
                                (gamma - X_2 %*% Beta_2 - S %*% lambda))
        tau_2 <- as.numeric(rnorm(1, mean_tau_2, sqrt(var_tau_2)))
        # tau_2 <- 1

        # Update Mu_1
        Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
        
        # Update Mu_2
        Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + S %*% lambda)
        
        # Pred_gaus
        preds_gaus <- predX %*% Beta_1 + tau_1*predS %*% eta

        # Pred_bios predicted probabilities
        logit_bios <- predX %*% Beta_2 + tau_2*predS %*% eta + predS %*% lambda
        preds_bios <- plogis(logit_bios)

        setTxtProgressBar(pb, index)
        if (index > nburn && (index - nburn) %% nthin == 0) {
        tau_1.chain[(index - nburn) / nthin] <- tau_1
        tau_2.chain[(index - nburn) / nthin] <- tau_2
        Beta_1.chain[, (index - nburn) / nthin] <- Beta_1
        Beta_2.chain[, (index - nburn) / nthin] <- Beta_2
        lambda.chain[, (index - nburn) / nthin] <- lambda
        eta.chain[, (index - nburn) / nthin] <- eta
        sig2.chain[(index - nburn) / nthin] <- sig2
        Sigma2_lambda.chain[(index - nburn) / nthin] <- sig2l
        Sigma2_eta.chain[(index - nburn) / nthin] <- sig2e
        Mu_1.chain[,(index-nburn)/nthin] <- Mu_1
        Mu_2.chain[,(index-nburn)/nthin] <- Mu_2
        preds_gaus.chain[,(index-nburn)/nthin] <- preds_gaus
        preds_bios.chain[,(index-nburn)/nthin] <- preds_bios
        logit_bios.chain[,(index-nburn)/nthin] <- logit_bios
        }
    }
    
    list(Beta_1.chain = Beta_1.chain, Beta_2.chain = Beta_2.chain, 
        lambda.chain = lambda.chain, eta.chain = eta.chain,
        Sigma2_lambda.chain = Sigma2_lambda.chain, Sigma2_eta.chain = Sigma2_eta.chain, sig2.chain = sig2.chain,
        Mu_1.chain = Mu_1.chain, Mu_2.chain = Mu_2.chain, 
        preds_gaus.chain = preds_gaus.chain, preds_bios.chain = preds_bios.chain,
        logit_bios.chain = logit_bios.chain,
        tau_1.chain = tau_1.chain, tau_2.chain = tau_2.chain)
}

## Multitype spatial model
MTSM_gr <- function(X_1, X_2, Z_1, Z_2, S, sig2b = 1000, wgt = NULL, n = NULL, 
                predX, predS, n_preds, nburn = 1000, nsim = 5000, nthin = 1, 
                sig2t = 10, sig2e = 10, tau_1_init = 1, tau_2_init = 1,  
                a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1, 
                a = 0.1, b = 0.1) {
    ## Z_1 is Gaussian response,
    ## Z_2 is binomial response with effective sample size m
    ## X is the covariates matrix for binomial response Z
    ## S is the basis function for spatial information
    ## wgt is the weights
    ## sig2b is the prior variance for beta 
    ## n is the total number parameter for binomial response, should be 1 if this is binomial
    ## nburn is burning in 
    ## nsim is the simulation number
    ## aeta, beta is hyperameter for sig2e
    ## alambda, blambda is hyperameter for sig2l
    N <- dim(X_1)[1] # Number of observations(or units)
    r <- ncol(S)
    npred <- dim(predX)[1]
    if(is.null(wgt)) wgt <- rep(1, N)
    if(is.null(n)) n <- rep(1, N)
    w <- sum(wgt)
    Wgt <- Diagonal(x=wgt)
    p_1 <- dim(X_1)[2] # length for Beta_1
    p_2 <- dim(X_2)[2] # length for Beta_2
    r <- dim(S)[2] # length for random effect
    Ip1 <- Diagonal(p_1) ## identity matrix of p dimention
    Ip2 <- Diagonal(p_2) ## identity matrix of p2 dimention
    Ir <- Diagonal(r) ## identity matrix of r dimention
    t_X_1 <- t(X_1)
    t_X_2 <- t(X_2)
    k <- wgt * (Z_2 - 1/2)

    # Initial values
    tau_1 <- tau_1_init
    tau_2 <- tau_2_init
    Beta_1 <- rep(1, p_1)
    Beta_2 <- rep(1, p_2)
    eta <- rep(1, r)
    Mu_1 <- rep(1, N)
    Mu_2 <- rep(1, N)
    sig2 <- 1
    sig2e <- sig2e
    Zeta <- rep(0, r)
    Sig2_Zeta <- 1


    # Chain containers  
    tau_1.chain <- array(0, dim = c(1, nsim / nthin))
    tau_2.chain <- array(0, dim = c(1, nsim / nthin))
    Beta_1.chain <- array(0, dim = c(p_1, nsim / nthin))
    Beta_2.chain <- array(0, dim = c(p_2, nsim / nthin))
    Sigma2_eta.chain <- rep(0, nsim / nthin)
    Zeta.chain <- array(0, dim = c(r, nsim / nthin))
    Sig2_Zeta.chain <- rep(0, nsim / nthin)
    sig2.chain <- rep(0, nsim / nthin)
    eta.chain <- array(0, dim = c(r, nsim / nthin))
    Mu_1.chain <- array(0, dim = c(N, nsim / nthin))
    Mu_2.chain <- array(0, dim = c(N, nsim / nthin))
    preds_gaus.chain <- array(0, dim = c(npred, nsim / nthin))
    preds_bios.chain <- array(0, dim = c(npred, nsim / nthin))
    logit_bios.chain <- array(0, dim = c(npred, nsim / nthin))
    print(paste0("Starting ",nsim / nthin, " iterations."))
    pb <- txtProgressBar(min=0, max=(nsim + nburn) / nthin, style=3)
    
    for (index in 1:(nsim + nburn)) {
        if (index %% 10000 == 0) cat(index, "\n")
        # Update sig2
        a_star <- as.numeric(a_eps + w/2)
        b_star <- as.numeric(b_eps + 0.5 * t(Z_1-Mu_1)%*%Wgt%*%(Z_1-Mu_1))
        sig2 <- 1/rgamma(1, shape = a_star, rate = b_star)
        d <- Wgt/sig2
        SD <- t(S) %*% d
        M <- SD %*% S
        XD <- t_X_1 %*% d %*% X_1
     
        # Update latent variable w
        omega <- rpg.gamma(N, wgt, Mu_2)
        Omega <- Diagonal(x=omega)
        SO <- t(S) %*% Omega
        OM <- SO %*% S 
        gamma <- k/omega

        # Update sigeta
        a_eta <- aeta + r/2
        b_eta <- beta + 0.5 * t(eta) %*% eta
        sig2e <- 1 / rgamma(1, shape = a_eta, rate = b_eta)

        # Update eta
        var.eta <- solve(tau_1 * M * tau_1 +
                        Ir/sig2e +
                        tau_2 * OM * tau_2, sparse = TRUE)
        mean.eta <- var.eta %*%
        (tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1 - S %*% Zeta) +
            tau_2 * SO %*% (gamma - X_2 %*% Beta_2))
        eta <- as.vector(rmvn(1, as.vector(mean.eta), var.eta))
        # eta <- eta - mean(eta) ## sum to zero

        # Update Beta_1
        var.Beta_1 <- solve(XD + Ip1/sig2b)
        mean.Beta_1 <- var.Beta_1 %*% t_X_1 %*% d %*% (Z_1 - tau_1 * S %*% eta - S %*% Zeta)
        Beta_1 <- as.vector(mvrnorm(1, mean.Beta_1, var.Beta_1))

        # Update Beta_2
        var.Beta_2 <- solve(t_X_2 %*% Omega %*% X_2 + Ip2/sig2b)
        mean.Beta_2 <- var.Beta_2 %*% t_X_2 %*% Omega %*% (gamma - 
                                                                  tau_2 * S %*% eta)
        Beta_2 <- as.vector(mvrnorm(1, mean.Beta_2, var.Beta_2))
        
        # Update regression parameter tau_1
        var_tau_1 <- as.numeric(solve(t(eta) %*% M %*% eta +
                                        1/sig2t))
        mean_tau_1 <- as.numeric(var_tau_1 * t(eta) %*% SD %*%
                                (Z_1 - X_1 %*% Beta_1 - S %*% Zeta))
        tau_1 <- as.numeric(rnorm(1, mean_tau_1, sqrt(var_tau_1)))
        # tau_1 <- 1

        # Update regression parameter tau_2
        # var_tau_2 <- as.numeric(solve(t(eta) %*% OM %*% eta +
        #                                 1/sig2t))
        # mean_tau_2 <- as.numeric(var_tau_2 * t(eta) %*% SO %*%
        #                                         (gamma - X_2 %*% Beta_2))
        # tau_2 <- as.numeric(rnorm(1, mean_tau_2, sqrt(var_tau_2)))
        tau_2 <- 1
        
        # Update Zeta
        var.Zeta <- solve(M + Ir/Sig2_Zeta)
        mean.Zeta <- var.Zeta %*% SD %*% (Z_1 - X_1 %*% Beta_1 - tau_1 * S %*% eta)
        Zeta <- as.vector(rmvn(1, as.vector(mean.Zeta), var.Zeta))
        # Zeta <- Zeta - mean(Zeta)

        # Update Sig2_Zeta
        a_zeta <- a + r/2
        b_zeta <- b + 0.5*t(Zeta) %*% Zeta
        Sig2_Zeta <- 1 / rgamma(1, shape = a_zeta, rate = b_zeta)

        # Update Mu_1
        Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + S %*% Zeta)
        
        # Update Mu_2
        Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta)
        
        # Pred_gaus
        preds_gaus <- predX %*% Beta_1 + tau_1 * predS %*% eta + predS %*% Zeta

        # Pred_bios predicted probabilities
        logit_bios <- predX %*% Beta_2 + tau_2 * predS %*% eta
        preds_bios <- plogis(logit_bios)

        setTxtProgressBar(pb, index)
        if (index > nburn && (index - nburn) %% nthin == 0) {
        tau_1.chain[(index - nburn) / nthin] <- tau_1
        tau_2.chain[(index - nburn) / nthin] <- tau_2
        Beta_1.chain[, (index - nburn) / nthin] <- Beta_1
        Beta_2.chain[, (index - nburn) / nthin] <- Beta_2
        eta.chain[, (index - nburn) / nthin] <- eta
        sig2.chain[(index - nburn) / nthin] <- sig2
        Sigma2_eta.chain[(index - nburn) / nthin] <- sig2e
        Zeta.chain[,(index-nburn)/nthin] <- Zeta
        Sig2_Zeta.chain[(index-nburn)/nthin] <- Sig2_Zeta
        Mu_1.chain[,(index-nburn)/nthin] <- Mu_1
        Mu_2.chain[,(index-nburn)/nthin] <- Mu_2
        preds_gaus.chain[,(index-nburn)/nthin] <- preds_gaus
        preds_bios.chain[,(index-nburn)/nthin] <- preds_bios
        logit_bios.chain[,(index-nburn)/nthin] <- logit_bios
        }
    }
    
    list(Beta_1.chain = Beta_1.chain, Beta_2.chain = Beta_2.chain, 
        eta.chain = eta.chain, Zeta.chain = Zeta.chain, Sig2_Zeta.chain = Sig2_Zeta.chain,
        Sigma2_eta.chain = Sigma2_eta.chain, sig2.chain = sig2.chain,
        Mu_1.chain = Mu_1.chain, Mu_2.chain = Mu_2.chain, 
        preds_gaus.chain = preds_gaus.chain, preds_bios.chain = preds_bios.chain,
        logit_bios.chain = logit_bios.chain,
        tau_1.chain = tau_1.chain, tau_2.chain = tau_2.chain)
}


# MTSM <- function(X_1, X_2, Z_1, Z_2, S, sig2b = 1000, wgt = NULL, n = NULL, 
#                 predX, predS, n_preds, nburn = 1000, nsim = 5000, nthin = 1, 
#                 sig2t = 10, sig2e = 10, tau_1_init = 1, tau_2_init = 1,  
#                 a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1, 
#                 alambda = 0.1, blambda = 0.1, a = 0.1, b = 0.1,
#                 l1 = 1,l2 = 1) {
#     ## Z_1 is Gaussian response,
#     ## Z_2 is binomial response with effective sample size m
#     ## X is the covariates matrix for binomial response Z
#     ## S is the basis function for spatial information
#     ## wgt is the weights
#     ## sig2b is the prior variance for beta 
#     ## n is the total number parameter for binomial response, should be 1 if this is binomial
#     ## nburn is burning in 
#     ## nsim is the simulation number
#     ## aeta, beta is hyperameter for sig2e
#     ## alambda, blambda is hyperameter for sig2l
#     N <- dim(X_1)[1] # Number of observations(or units)
#     r <- ncol(S)
#     npred <- dim(predX)[1]
#     if(is.null(wgt)) wgt <- rep(1, N)
#     if(is.null(n)) n <- rep(1, N)
#     w <- sum(wgt)
#     Wgt <- Diagonal(x=wgt)
#     p_1 <- dim(X_1)[2] # length for Beta_1
#     p_2 <- dim(X_2)[2] # length for Beta_2
#     r <- dim(S)[2] # length for random effect
#     Ip1 <- Diagonal(p_1) ## identity matrix of p dimention
#     Ip2 <- Diagonal(p_2) ## identity matrix of p2 dimention
#     Ir <- Diagonal(r) ## identity matrix of r dimention
#     t_X_1 <- t(X_1)
#     t_X_2 <- t(X_2)
#     k <- wgt * (Z_2 - 1/2)

#     # Initial values
#     tau_1 <- tau_1_init
#     tau_2 <- tau_2_init
#     Beta_1 <- rep(1, p_1)
#     Beta_2 <- rep(1, p_2)
#     lambda <- rep(0, r)
#     eta <- rep(0, r)
#     Mu_1 <- rep(1, N)
#     Mu_2 <- rep(1, N)
#     sig2 <- 1
#     sig2e <- sig2e
#     sig2l <- 1
#     Zeta <- rep(0, r)
#     Sig2_Zeta <- 1


#     # Chain containers  
#     tau_1.chain <- array(0, dim = c(1, nsim / nthin))
#     tau_2.chain <- array(0, dim = c(1, nsim / nthin))
#     Beta_1.chain <- array(0, dim = c(p_1, nsim / nthin))
#     Beta_2.chain <- array(0, dim = c(p_2, nsim / nthin))
#     Sigma2_lambda.chain <- rep(0, nsim / nthin)
#     Sigma2_eta.chain <- rep(0, nsim / nthin)
#     Zeta.chain <- array(0, dim = c(r, nsim / nthin))
#     Sig2_Zeta.chain <- rep(0, nsim / nthin)
#     sig2.chain <- rep(0, nsim / nthin)
#     lambda.chain <- array(0, dim = c(r, nsim / nthin))
#     eta.chain <- array(0, dim = c(r, nsim / nthin))
#     Mu_1.chain <- array(0, dim = c(N, nsim / nthin))
#     Mu_2.chain <- array(0, dim = c(N, nsim / nthin))
#     preds_gaus.chain <- array(0, dim = c(npred, nsim / nthin))
#     preds_bios.chain <- array(0, dim = c(npred, nsim / nthin))
#     print(paste0("Starting ",nsim / nthin, " iterations."))
#     pb <- txtProgressBar(min=0, max=(nsim + nburn) / nthin, style=3)
    
#     for (index in 1:(nsim + nburn)) {
#         if (index %% 10000 == 0) cat(index, "\n")
#         # Update sig2
#         a_star <- as.numeric(a_eps + w/2)
#         b_star <- as.numeric(b_eps + 0.5 * t(Z_1-Mu_1)%*%Wgt%*%(Z_1-Mu_1))
#         sig2 <- 1/rgamma(1, shape = a_star, rate = b_star)
#         d <- Wgt/sig2
#         SD <- t(S) %*% d
#         M <- SD %*% S
#         XD <- t_X_1 %*% d %*% X_1
     
#         # Update latent variable w
#         omega <- rpg.gamma(N, wgt, Mu_2)
#         Omega <- Diagonal(x=omega)
#         SO <- t(S) %*% Omega
#         OM <- SO %*% S 
#         gamma <- k/omega
                
#         # # Update sig2eta
#         # a_eta <- aeta + r / 2
#         # b_eta <- beta + 0.5 * t(eta) %*% eta
#         # sig2e <- 1 / rgamma(1, shape = a_eta, rate = b_eta)

#         # Update eta
#         var.eta <- solve(tau_1 * M * tau_1 +
#                         Ir/sig2e +
#                         tau_2 * OM * tau_2, sparse = TRUE)
#         mean.eta <- var.eta %*%
#         (tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1 - l1*S %*% Zeta) +
#             tau_2 * SO %*% (gamma - X_2 %*% Beta_2 - l2*S %*% lambda))
#         eta <- as.vector(rmvn(1, as.vector(mean.eta), var.eta))
#         eta <- eta - mean(eta) ## sum to zero

#         # Update Beta_1
#         var.Beta_1 <- solve(XD + Ip1/sig2b)
#         mean.Beta_1 <- var.Beta_1 %*% t_X_1 %*% d %*% (Z_1 - tau_1 * S %*% eta - l1*S %*% Zeta)
#         Beta_1 <- as.vector(mvrnorm(1, mean.Beta_1, var.Beta_1))
                
#         # Update Zeta
#         var.Zeta <- solve(M + Ir/Sig2_Zeta)
#         mean.Zeta <- var.Zeta %*% SD %*% (Z_1 - X_1 %*% Beta_1 - tau_1 * S %*% eta)
#         Zeta <- as.vector(mvrnorm(1, as.vector(mean.Zeta), var.Zeta))

#         # Update Sig2_Zeta
#         a_zeta <- a + r/2
#         b_zeta <- b + 0.5*t(Zeta)*Zeta
#         Sig2_Zeta <- 1 / rgamma(1, shape = a_zeta, rate = b_zeta)

#         # Update sig2l
#         a_l <- alambda + r / 2
#         b_l <- blambda + 0.5 * t(lambda) %*% lambda
#         sig2l <- 1 / rgamma(1, shape = a_l, rate = b_l)

#         # Update lambda
#         var.lambda <- solve(OM + Ir/sig2l, sparse = TRUE)
#         mean.lambda <- var.lambda %*% (SO %*%(gamma -
#                                                       X_2 %*% Beta_2 - tau_2 * S %*% eta))
#         lambda <- as.vector(rmvn(1, as.vector(mean.lambda), var.lambda))
#         lambda <- lambda - mean(lambda) ## sum to zero

#         # Update Beta_2
#         var.Beta_2 <- solve(t_X_2 %*% Omega %*% X_2 + Ip2/sig2b)
#         mean.Beta_2 <- var.Beta_2 %*% t_X_2 %*% Omega %*% (gamma - 
#                                                                   tau_2 * S %*% eta - l2*S %*% lambda)
#         Beta_2 <- as.vector(mvrnorm(1, mean.Beta_2, var.Beta_2))
        
#         # Update regression parameter tau_1
#         var_tau_1 <- as.numeric(solve(t(eta) %*% M %*% eta +
#                                         1/sig2t))
#         mean_tau_1 <- as.numeric(var_tau_1 * t(eta) %*% SD %*%
#                                 (Z_1 - X_1 %*% Beta_1 - l1*S %*% Zeta))
#         tau_1 <- as.numeric(rnorm(1, mean_tau_1, sqrt(var_tau_1)))
        
#         # Update regression parameter tau_2
#         var_tau_2 <- as.numeric(solve(t(eta) %*% OM %*% eta +
#                                         1/sig2t))
#         mean_tau_2 <- as.numeric(var_tau_2 * t(eta) %*% SO %*%
#                                                 (gamma - X_2 %*% Beta_2 -l2* S %*% lambda))
#         tau_2 <- as.numeric(rnorm(1, mean_tau_2, sqrt(var_tau_2)))
        
#         # Update Mu_1
#         Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + l1*S %*% Zeta)
        
#         # Update Mu_2
#         Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + l2*S %*% lambda)
        
#         # Pred_gaus
#         preds_gaus <- predX %*% Beta_1 + tau_1*predS %*% eta + l1*predS %*% Zeta

#         # Pred_bios predicted probabilities
#         preds_bios <- plogis(predX %*% Beta_2 + tau_2*predS %*% eta + l2*predS %*% lambda)

#         setTxtProgressBar(pb, index)
#         if (index > nburn && (index - nburn) %% nthin == 0) {
#         tau_1.chain[(index - nburn) / nthin] <- tau_1
#         tau_2.chain[(index - nburn) / nthin] <- tau_2
#         Beta_1.chain[, (index - nburn) / nthin] <- Beta_1
#         Beta_2.chain[, (index - nburn) / nthin] <- Beta_2
#         lambda.chain[, (index - nburn) / nthin] <- lambda
#         eta.chain[, (index - nburn) / nthin] <- eta
#         sig2.chain[(index - nburn) / nthin] <- sig2
#         Sigma2_lambda.chain[(index - nburn) / nthin] <- sig2l
#         Sigma2_eta.chain[(index - nburn) / nthin] <- sig2e
#         Zeta.chain[,(index-nburn)/nthin] <- Zeta
#         Sig2_Zeta.chain[(index-nburn)/nthin] <- Sig2_Zeta
#         Mu_1.chain[,(index-nburn)/nthin] <- Mu_1
#         Mu_2.chain[,(index-nburn)/nthin] <- Mu_2
#         preds_gaus.chain[,(index-nburn)/nthin] <- preds_gaus
#         preds_bios.chain[,(index-nburn)/nthin] <- preds_bios
#         }
#     }
    
#     list(Beta_1.chain = Beta_1.chain, Beta_2.chain = Beta_2.chain, 
#         lambda.chain = lambda.chain, eta.chain = eta.chain, Zeta.chain = Zeta.chain, Sig2_Zeta.chain = Sig2_Zeta.chain,
#         Sigma2_lambda.chain = Sigma2_lambda.chain, Sigma2_eta.chain = Sigma2_eta.chain, sig2.chain = sig2.chain,
#         Mu_1.chain = Mu_1.chain, Mu_2.chain = Mu_2.chain, 
#         preds_gaus.chain = preds_gaus.chain, preds_bios.chain = preds_bios.chain,
#         tau_1.chain = tau_1.chain, tau_2.chain = tau_2.chain)
# }
MTSM <- function(X_1, X_2, Z_1, Z_2, S, sig2b = 1000, wgt = NULL, n = NULL, predX, predS, nburn = 1000, nsim = 5000, nthin = 1,  tau_1_init = 1, tau_2_init = 1, tau_3_init = 1, sig2t = 10, sig2e = 1, a_eps = 0.1, b_eps = 0.1, aeta = 0.1, beta = 0.1, alambda = 0.1, blambda = 0.1) {
    ## Z_1 is Gaussian response,
    ## Z_2 is binomial response with effective sample size m
    ## X is the covariates matrix for binomial response Z
    ## S is the basis function for spatial information
    ## wgt is the weights
    ## sig2b is the prior variance for beta 
    ## n is the total number parameter for binomial response, should be 1 if this is binomial
    ## nburn is burning in 
    ## nsim is the simulation number
    ## aeta, beta is hyperameter for sig2e
    ## alambda, blambda is hyperameter for sig2l
    N <- dim(X_1)[1] # Number of observations(or units)
    r <- ncol(S)
    npred <- dim(predX)[1]
    if(is.null(wgt)) wgt <- rep(1, N)
    if(is.null(n)) n <- rep(1, N)
    w <- sum(wgt)
    Wgt <- Diagonal(x=wgt)
    p_1 <- dim(X_1)[2] # length for Beta_1
    p_2 <- dim(X_2)[2] # length for Beta_2
    r <- dim(S)[2] # length for random effect
    Ip1 <- Diagonal(p_1) ## identity matrix of p dimention
    Ip2 <- Diagonal(p_2) ## identity matrix of p2 dimention
    Ir <- Diagonal(r) ## identity matrix of r dimention
    t_X_1 <- t(X_1)
    t_X_2 <- t(X_2)
    k <- wgt * (Z_2 - n/2)

    # Initial values
    tau_1 <- tau_1_init
    tau_2 <- tau_2_init
    tau_3 <- tau_3_init
    Beta_1 <- rep(1, p_1)
    Beta_2 <- rep(1, p_2)
    lambda <- rep(0, r)
    Zeta <- rep(0, r)
    eta <- rep(1, r)
    Mu_1 <- rep(1, N)
    Mu_2 <- rep(1, N)
    sig2 <- 10
    sig2e <- sig2e
    sig2l <- 10
    sig2t <- sig2t
    sig2z <- 10


    # Chain containers  
    tau_1.chain <- array(0, dim = c(1, nsim / nthin))
    tau_2.chain <- array(0, dim = c(1, nsim / nthin))
    tau_3.chain <- array(0, dim = c(1, nsim / nthin))
    Beta_1.chain <- array(0, dim = c(p_1, nsim / nthin))
    Beta_2.chain <- array(0, dim = c(p_2, nsim / nthin))
    Sigma2_lambda.chain <- rep(0, nsim / nthin)
    Sigma2_Zeta.chain <- rep(0, nsim / nthin)
    Sigma2_eta.chain <- rep(0, nsim / nthin)
    sig2.chain <- rep(0, nsim / nthin)
    lambda.chain <- array(0, dim = c(r, nsim / nthin))
    Zeta.chain <- array(0, dim = c(r, nsim / nthin))
    eta.chain <- array(0, dim = c(r, nsim / nthin))
    Mu_1.chain <- array(0, dim = c(N, nsim / nthin))
    Mu_2.chain <- array(0, dim = c(N, nsim / nthin))
    preds_gaus.chain <- array(0, dim = c(npred, nsim / nthin))
    preds_bios.chain <- array(0, dim = c(npred, nsim / nthin))
    print(paste0("Starting ",nsim / nthin, " iterations."))
    pb <- txtProgressBar(min=0, max=(nsim + nburn) / nthin, style=3)
    
    for (index in 1:(nsim + nburn)) {
        if (index %% 10000 == 0) cat(index, "\n")
        # Update sig2
        a_star <- as.numeric(a_eps + w/2)
        b_star <- as.numeric(b_eps + 0.5 * t(Z_1-Mu_1)%*%Wgt%*%(Z_1-Mu_1))
        sig2 <- 1/rgamma(1, shape = a_star, rate = b_star)
        d <- Wgt/sig2
        SD <- t(S) %*% d
        M <- SD %*% S
        XD <- t_X_1 %*% d %*% X_1
     
        # Update latent variable w
        omega <- rpg.gamma(N, wgt, Mu_2)
        Omega <- Diagonal(x=omega)
        SO <- t(S) %*% Omega
        OM <- SO %*% S 
        gamma <- k/omega
                
        # Update sig2eta
        # a_eta <- aeta + r / 2
        # b_eta <- beta + 0.5 * t(eta) %*% eta
        # sig2e <- 1 / rgamma(1, shape = a_eta, rate = b_eta)

        # Update eta
        var.eta <- solve(tau_1 * M * tau_1 +
                        Ir/sig2e +
                        tau_2 * OM * tau_2, sparse = TRUE)
        mean.eta <- var.eta %*%
        (tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1 - S %*% Zeta) +
            tau_2 * SO %*% (gamma - X_2 %*% Beta_2 - S %*% lambda))
        eta <- as.vector(rmvn(1, as.vector(mean.eta), var.eta))
        eta <- eta - mean(eta) ## sum to zero
        # var.eta <- solve(tau_1 * M * tau_1 +
        #                 Ir/sig2e +
        #                 tau_2 * OM * tau_2, sparse = TRUE)
        # mean.eta <- var.eta %*%
        # (tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1) +
        #     tau_2 * SO %*% (gamma - X_2 %*% Beta_2 - S %*% lambda))
        # eta <- as.vector(rmvn(1, as.vector(mean.eta), var.eta))
        # eta <- eta - mean(eta) ## sum to zero

        # Update Beta_1
        var.Beta_1 <- solve(XD + Ip1/sig2b)
        mean.Beta_1 <- var.Beta_1 %*% t_X_1 %*% d %*% (Z_1 - tau_1 * S %*% eta - S %*% Zeta)
        Beta_1 <- as.vector(mvrnorm(1, mean.Beta_1, var.Beta_1))
        # var.Beta_1 <- solve(XD + Ip1/sig2b)
        # mean.Beta_1 <- var.Beta_1 %*% t_X_1 %*% d %*% (Z_1 - tau_1 * S %*% eta)
        # Beta_1 <- as.vector(mvrnorm(1, mean.Beta_1, var.Beta_1))
        
        # Update Beta_2
        var.Beta_2 <- solve(t_X_2 %*% Omega %*% X_2 + Ip2/sig2b)
        mean.Beta_2 <- var.Beta_2 %*% t_X_2 %*% Omega %*% (gamma - tau_2 * S %*% eta - S %*% lambda)
        Beta_2 <- as.vector(mvrnorm(1, mean.Beta_2, var.Beta_2))

        # Update sig2l
        a_l <- alambda + r / 2
        b_l <- blambda + 0.5 * t(lambda) %*% lambda
        sig2l <- 1 / rgamma(1, shape = a_l, rate = b_l)

        # Update lambda
        var.lambda <- solve(tau_3 * OM * tau_3  + Ir/sig2l, sparse = TRUE)
        mean.lambda <- var.lambda %*% (tau_3 * SO %*%(gamma -
                                                      X_2 %*% Beta_2 - tau_2 * S %*% eta))
        lambda <- as.vector(rmvn(1, as.vector(mean.lambda), var.lambda))
        lambda <- lambda - mean(lambda) ## sum to zero

        # Update Sig2_Zeta
        a_zeta <- alambda + r/2
        b_zeta <- blambda + 0.5*t(Zeta) %*% Zeta
        Sig2_Zeta <- 1 / rgamma(1, shape = a_zeta, rate = b_zeta)
        
        # Update Zeta
        var.Zeta <- solve(M + Ir/Sig2_Zeta)
        mean.Zeta <- var.Zeta %*% SD %*% (Z_1 - X_1 %*% Beta_1 - tau_1 * S %*% eta)
        Zeta <- as.vector(rmvn(1, as.vector(mean.Zeta), var.Zeta))
        Zeta <- Zeta - mean(Zeta)
        
        # Update regression parameter tau_1
        var_tau_1 <- as.numeric(solve(t(eta) %*% M %*% eta +
                                        1/sig2t))
        mean_tau_1 <- as.numeric(var_tau_1 * t(eta) %*% SD %*%
                                (Z_1 - X_1 %*% Beta_1 - S %*% Zeta))
        # var_tau_1 <- as.numeric(solve(t(eta) %*% M %*% eta +
        #                                 1/sig2t))
        # mean_tau_1 <- as.numeric(var_tau_1 * t(eta) %*% SD %*%
        #                         (Z_1 - X_1 %*% Beta_1))
        tau_1 <- as.numeric(rnorm(1, mean_tau_1, sqrt(var_tau_1)))
        tau_1 <- 1

        # Update regression parameter tau_2
        var_tau_2 <- as.numeric(solve(t(eta) %*% OM %*% eta +
                                        1/sig2t))
        mean_tau_2 <- as.numeric(var_tau_2 * t(eta) %*% SO %*%
                                                (gamma - X_2 %*% Beta_2 - S %*% lambda))
        tau_2 <- as.numeric(rnorm(1, mean_tau_2, sqrt(var_tau_2)))
        tau_2 <- 0

        # Update Mu_1
        Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + S %*% Zeta)
        # Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
        # Update Mu_2
        Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + S %*% lambda)
        
        # Pred_gaus
        preds_gaus <- predX %*% Beta_1 + tau_1 * predS %*% eta + predS %*% Zeta
        # preds_gaus <- predX %*% Beta_1 + tau_1 * predS %*% eta

        # Pred_bios predicted probabilities
        preds_bios <- plogis(predX %*% Beta_2 + tau_2 * predS %*% eta + predS %*% lambda)

        setTxtProgressBar(pb, index)
        if (index > nburn && (index - nburn) %% nthin == 0) {
        tau_1.chain[(index - nburn) / nthin] <- tau_1
        tau_2.chain[(index - nburn) / nthin] <- tau_2
        tau_3.chain[(index - nburn) / nthin] <- tau_3
        Beta_1.chain[, (index - nburn) / nthin] <- Beta_1
        Beta_2.chain[, (index - nburn) / nthin] <- Beta_2
        lambda.chain[, (index - nburn) / nthin] <- lambda
        Zeta.chain[,(index-nburn)/nthin] <- Zeta
        eta.chain[, (index - nburn) / nthin] <- eta
        sig2.chain[(index - nburn) / nthin] <- sig2
        Sigma2_lambda.chain[(index - nburn) / nthin] <- sig2l
        Sigma2_Zeta.chain[(index-nburn)/nthin] <- Sig2_Zeta
        Sigma2_eta.chain[(index - nburn) / nthin] <- sig2e
        Mu_1.chain[,(index-nburn)/nthin] <- Mu_1
        Mu_2.chain[,(index-nburn)/nthin] <- Mu_2
        preds_gaus.chain[,(index-nburn)/nthin] <- preds_gaus
        preds_bios.chain[,(index-nburn)/nthin] <- preds_bios
        }
    }
    
    list(Beta_1.chain = Beta_1.chain, Beta_2.chain = Beta_2.chain, 
        lambda.chain = lambda.chain, Zeta.chain = Zeta.chain, eta.chain = eta.chain, 
        Sigma2_lambda.chain = Sigma2_lambda.chain, Sigma2_eta.chain = Sigma2_eta.chain, 
        Sigma2_zeta.chain = Sigma2_Zeta.chain,sig2.chain = sig2.chain,
        Mu_1.chain = Mu_1.chain, Mu_2.chain = Mu_2.chain, 
        preds_gaus.chain = preds_gaus.chain, preds_bios.chain = preds_bios.chain,
        tau_1.chain = tau_1.chain, tau_2.chain = tau_2.chain, 
        tau_3.chain = tau_3.chain)
}

