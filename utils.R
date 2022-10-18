
library(Matrix)
library(psych)

# collapse parameters into vector
param2vect <- function(param, constraint)
{
  p <- length(param$psi_s[[1]])
  S <- length(param$psi_s)
  phi_vals <- as.vector(param$Phi[lower.tri(param$Phi, diag = TRUE)])
  k <- ncol(param$Phi)
  lambda_vals <- psi_vals <- j_s <- c()

  if(constraint=="block_lower1"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]][-(1:k),]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$psi_s[[s]])
      j_s[[s]] <- ncol(param$Lambda_s[[s]])}
  }
  
  if(constraint=="block_lower2"){
    for(s in 1:S){
      Lambda_s <- param$Lambda_s[[s]]
      lambda_vals <- c(lambda_vals, as.vector(Lambda_s[lower.tri(Lambda_s, diag = TRUE)]))
      psi_vals <- c(psi_vals, param$psi_s[[s]])
      j_s[[s]] <- ncol(param$Lambda_s[[s]])}
  }
  
  theta <- c(phi_vals, lambda_vals, psi_vals)
  
  return(theta)
}

# reconstitute parameters from vector
vect2param <- function(vect, param.struct, constraint, p, k, j_s)
{
  S <- length(j_s)
  nP <- k * p - k * ( k - 1) / 2
  if (constraint == "block_lower1")  { nL <- j_s * (p - k)  - j_s *  (j_s - 1) / 2}
  if (constraint == "block_lower2")  { nL <- j_s * p  - j_s *  (j_s - 1) / 2}
  phi_vals <- vect[1:nP]
  Phi <- matrix(0, nrow=p, ncol=k)
  Phi[lower.tri(Phi, diag = TRUE)] <- phi_vals
  Lambda_s <- param.struct$Lambda_s
  psi_s <- param.struct$psi_s
  
  for(s in 1:S){
    nL_s  <- if(s==1) 0 else sum(nL[1:(s-1)])
    ind <-  (nP + nL_s + 1):(nP + nL_s + nL[s])
    lambda_vals_s <-  vect[ind]
    Lambda_s[[s]][Lambda_s[[s]]!=0] <- lambda_vals_s
    ind_s <- (nP + sum(nL) + p * (s-1) + 1):(nP + sum(nL) + p * s)
    psi_vals_s <- vect[ind_s]
    psi_s[[s]] <-  psi_vals_s
  }
  
  return(list(Phi=Phi, Lambda_s=Lambda_s, psi_s=psi_s))
}

# E-Step
exp_values_msfam <- function(Phi, Lambda_s, 
                             Psi_s, Psi_s1, 
                             mu_s, n_s, X_s, 
                             incomplete_cols, 
                             incomplete_obs, 
                             getdet = FALSE,
                             to_calculate=c("Txis", "Txsxs", "Txsfs",
                                            "Txsls", "Tfsfs", "Tlsls",
                                            "Tfsls", "Txs",
                                            "Tfs", "Tls")
)
{
  k <- dim(Phi)[2]
  S <- length(Lambda_s)
  p <- dim(Phi)[1]

  # define objects
  j_s <- numeric(S)
  Sig_s <- vector(mode="list", length=S)
  ds_s <- vector(mode="list", length=S)
  Sig_s1 <- vector(mode="list", length=S)
  Txis <- vector(mode="list", length=S)  # this one is the raw data, imputed with conditional mean
  Txsxs <- vector(mode="list", length=S)
  Txsfs <- vector(mode="list", length=S)
  Txsls <- vector(mode="list", length=S)
  Tfsfs <- vector(mode="list", length=S)
  Tlsls <- vector(mode="list", length=S)
  Tfsls <- vector(mode="list", length=S)
  Txs <- vector(mode="list", length=S)
  Tfs <- vector(mode="list", length=S)
  Tls <- vector(mode="list", length=S)

  for (s in 1:S){

    num_complete_cases <- length(incomplete_obs[[s]][["complete"]])

    dat <- X_s[[s]]

    j_s[[s]] <- ncol(Lambda_s[[s]])

    Sig_s[[s]] <- Phi %*% t(Phi) +
      Lambda_s[[s]] %*% t(Lambda_s[[s]]) +
      Psi_s[[s]]

    ds_s[[s]] <- determinant(Sig_s[[s]], logarithm=FALSE)$modulus[[1]]

    Sig_s1[[s]] <- chol2inv(chol(Sig_s[[s]]))
    # ^-- chol2inv(chol()) is faster than solve()
    # Sig_s1[[s]] <- solve(Sig_s[[s]])

    # conditional expectation at individual level
    # will later take mean over all
    Txisxis <- Txisf <- Txisl <- Tff <- Tll <- Tfl <- Tf <- Tl <-
      vector(mode="list", length=n_s[s])

    for (i in incomplete_obs[[s]][["incomplete"]]) {
      # columns that are missing for person i
      miss.idx.i <- incomplete_cols[[s]][["missing"]][[i]]

      # columns that are observed for person i
      obs.idx.i <- incomplete_cols[[s]][["observed"]][[i]]

      # record order of missing and observed variable indices for person i
      correct_order <- order(c(miss.idx.i, obs.idx.i))

      # if any columns are not observed,
      # note which columns and
      # generate the "latent adjustment"
      # Exihs = mu_hs + Sigma_hvs Sigma_vvs^{-1} (x_ivs - mu_vs)
      # == conditional expectation of the missing value

      xivs <- dat[i, obs.idx.i]
      mu_hs <- mu_s[[s]][miss.idx.i]
      mu_vs <- mu_s[[s]][obs.idx.i]
      Sig_hhs <- Sig_s[[s]][miss.idx.i, miss.idx.i]
      Sig_hvs <- Sig_s[[s]][miss.idx.i, obs.idx.i]
      Sig_vvs_inv <- chol2inv(chol(Sig_s[[s]][obs.idx.i, obs.idx.i]))
      # ^-- chol2inv(chol()) is faster than solve()
      # Sig_vvs_inv <- solve(Sig_s[[s]][obs.idx.i, obs.idx.i])
      Sig_vhs <- Sig_s[[s]][obs.idx.i, miss.idx.i]

      # latent adjustment
      # == Exihs
      adj_i <- mu_hs +
        Sig_hvs %*%
        Sig_vvs_inv %*%
        (xivs - mu_vs)

      # Txisxis
      if ("Txsxs" %in% to_calculate) {
        Exihs <- adj_i
        Exihsxihs <- Sig_hhs -
          Sig_hvs %*% Sig_vvs_inv %*% Sig_vhs +
          adj_i %*% t(adj_i)
        Txisxis[[i]] <- matrix(NA, nrow=p, ncol=p)
        Txisxis[[i]][1:length(miss.idx.i), 1:length(miss.idx.i)] <- Exihsxihs
        Txisxis[[i]][(length(miss.idx.i) + 1):p, (length(miss.idx.i) + 1):p] <- xivs %*% t(xivs)
        Txisxis[[i]][(length(miss.idx.i) + 1):p, 1:length(miss.idx.i)] <- xivs %*% t(Exihs)
        Txisxis[[i]][upper.tri(Txisxis[[i]], diag=FALSE)] <- t(Txisxis[[i]])[upper.tri(Txisxis[[i]], diag=FALSE)]
        Txisxis[[i]] <- Txisxis[[i]][correct_order, correct_order]
      }

      # Txisf
      if ("Txsfs" %in% to_calculate) {
        Exihsf <- Phi[miss.idx.i, ] -
          Sig_hvs %*% Sig_vvs_inv %*% Phi[obs.idx.i, ] +
          adj_i %*% t(t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))

        Exivsf <- xivs %*% t(t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))

        Txisf[[i]] <- rbind(Exihsf, Exivsf)
        Txisf[[i]] <- Txisf[[i]][correct_order, ]
      }

      # Txisl
      if ("Txsls" %in% to_calculate) {
        Exihsl <- Lambda_s[[s]][miss.idx.i, ] -
          Sig_hvs %*% Sig_vvs_inv %*% Lambda_s[[s]][obs.idx.i, ] +
          adj_i %*% t(t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))

        Exivsl <- xivs %*% t(t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))

        Txisl[[i]] <- rbind(Exihsl, Exivsl)
        Txisl[[i]] <- Txisl[[i]][correct_order, ]
      }

      # Tff
      if ("Tfsfs" %in% to_calculate) {
        Tff[[i]] <- diag(k) -
          t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% Phi[obs.idx.i, ] +
          (t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs)) %*%
          t(t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))
      }

      # Tll
      if ("Tlsls" %in% to_calculate) {
        Tll[[i]] <- diag(j_s[[s]]) -
          t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% Lambda_s[[s]][obs.idx.i, ] +
          (t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs)) %*%
          t(t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))
      }

      # Tfl
      if ("Tfsls" %in% to_calculate) {
        Tfl[[i]] <- -t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% Lambda_s[[s]][obs.idx.i, ] +
          (t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs)) %*%
          t(t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs))
      }

      # Txis
      if (("Txs" %in% to_calculate) | ("Txis" %in% to_calculate)) {
        Exihs <- adj_i
        Exivs <- matrix(xivs, ncol=1)
        Txis[[s]][[i]] <- rbind(Exihs, Exivs)
        Txis[[s]][[i]] <- Txis[[s]][[i]][correct_order]
      }

      # Tf
      if ("Tfs" %in% to_calculate) {
        Tf[[i]] <- t(Phi[obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs)
      }

      # Tl
      if ("Tls" %in% to_calculate) {
        Tl[[i]] <- t(Lambda_s[[s]][obs.idx.i, ]) %*% Sig_vvs_inv %*% (xivs - mu_vs)
      }

    }

    # add up the contribution from incomplete observations
    incomplete_Txisxis_contribution <- Reduce('+', Filter(Negate(is.null), Txisxis))
    if (!is.null(incomplete_Txisxis_contribution)) Txsxs[[s]] <- incomplete_Txisxis_contribution

    # add up the contribution from incomplete observations
    incomplete_Txisf_contribution <- Reduce('+', Filter(Negate(is.null), Txisf))
    if (!is.null(incomplete_Txisf_contribution)) Txsfs[[s]] <- incomplete_Txisf_contribution

    # add up the contribution from incomplete observations
    incomplete_Txisl_contribution <- Reduce('+', Filter(Negate(is.null), Txisl))
    if (!is.null(incomplete_Txisl_contribution)) Txsls[[s]] <- incomplete_Txisl_contribution

    # add up the contribution from incomplete observations
    incomplete_Tff_contribution <- Reduce('+', Filter(Negate(is.null), Tff))
    if (!is.null(incomplete_Tff_contribution)) Tfsfs[[s]] <- incomplete_Tff_contribution
    
    # add up the contribution from incomplete observations
    incomplete_Tll_contribution <- Reduce('+', Filter(Negate(is.null), Tll))
    if (!is.null(incomplete_Tll_contribution)) Tlsls[[s]] <- incomplete_Tll_contribution

    # add up the contribution from incomplete observations
    incomplete_Tfl_contribution <- Reduce('+', Filter(Negate(is.null), Tfl))
    if (!is.null(incomplete_Tfl_contribution)) Tfsls[[s]] <- incomplete_Tfl_contribution

    # add up the contribution from incomplete observations
    incomplete_Tfs_contribution <- Reduce('+', Filter(Negate(is.null), Tf))
    if (!is.null(incomplete_Tfs_contribution)) Tfs[[s]] <- incomplete_Tfs_contribution

    # add up the contribution from incomplete observations
    incomplete_Tl_contribution <- Reduce('+', Filter(Negate(is.null), Tl))
    if (!is.null(incomplete_Tl_contribution)) Tls[[s]] <- incomplete_Tl_contribution

    # for observations that have all variables observed,
    # the only latent variables are f, l
    if (num_complete_cases > 0) {

      if ("Txsxs" %in% to_calculate) {
        
        if (is.null(incomplete_Txisxis_contribution)) {
          Txsxs[[s]] <- matrix(0, p, p)
        }

        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)
        Txsxs[[s]] <- Txsxs[[s]] + t(observed_vals_from_incomplete_obs) %*% observed_vals_from_incomplete_obs
      }

      if ("Txsfs" %in% to_calculate) {

        if (is.null(incomplete_Txisf_contribution)) {
          Txsfs[[s]] <- matrix(0, p, k)
        }
        
        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)
        
        Txsfs[[s]] <- Txsfs[[s]] + t(observed_vals_from_incomplete_obs) %*%
          t(t(Phi) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) -
                                          mu_s[[s]]))
      }

      if ("Txsls" %in% to_calculate) {

        if (is.null(incomplete_Txisl_contribution)) {
          Txsls[[s]] <- matrix(0, p, j_s[s])
        }

        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)
        
        Txsls[[s]] <- Txsls[[s]] + t(observed_vals_from_incomplete_obs) %*%
          t(t(Lambda_s[[s]]) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]]))
      }

      if ("Tfsfs" %in% to_calculate) {
        
        if (is.null(incomplete_Tff_contribution)) {
          Tfsfs[[s]] <- matrix(0, k, k)
        } 

        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)
        
        Tfsfs[[s]] <- Tfsfs[[s]] +
          num_complete_cases*(diag(k) - t(Phi) %*% Sig_s1[[s]] %*% Phi) +
          (t(Phi) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]])) %*%
          t(t(Phi) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]]))

      }


      if ("Tlsls" %in% to_calculate) {
        
        if (is.null(incomplete_Tll_contribution)) {
          Tlsls[[s]] <- matrix(0, j_s[s], j_s[s])
        }

        factor1 <- t(Lambda_s[[s]]) %*% Sig_s1[[s]]
        factor2 <- factor1 %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]])
        Tlsls[[s]] <- Tlsls[[s]] + num_complete_cases*(diag(j_s[[s]]) -
                                                         factor1 %*% Lambda_s[[s]]) +
          factor2 %*% t(factor2)

      }


      if ("Tfsls" %in% to_calculate) {
        
        incomplete_Tfl_contribution <- Reduce('+', Filter(Negate(is.null), Tfl)) # / n_s[s]
        
        if (is.null(incomplete_Tfl_contribution)) {
          Tfsls[[s]] <- matrix(0, k, j_s[s])
        }

        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)

        Tfsls[[s]] <- Tfsls[[s]] + num_complete_cases*(-t(Phi) %*% Sig_s1[[s]] %*% Lambda_s[[s]]) +
          (t(Phi) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]])) %*%
          t(t(Lambda_s[[s]]) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]]))

      }

      if (("Txs" %in% to_calculate) | ("Txis" %in% to_calculate)) {

        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)

        Txis[[s]][incomplete_obs[[s]][["complete"]]] <- lapply(split(observed_vals_from_incomplete_obs, 1:num_complete_cases), function(x) x)

      }

      if ("Tfs" %in% to_calculate) {
        
        if (is.null(incomplete_Tfs_contribution)) {
          Tfs[[s]] <- matrix(0, k, 1)
        }

        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)
        

        Tfs[[s]] <- Tfs[[s]] + matrix(rowSums(t(Phi) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]])), ncol=1)

      }

      if ("Tls" %in% to_calculate) {
        incomplete_Tl_contribution <- Reduce('+', Filter(Negate(is.null), Tl)) # / n_s[s]
        if (is.null(incomplete_Tl_contribution)) {
          Tls[[s]] <- matrix(0, j_s[s], 1)
        }
        
        observed_vals_from_incomplete_obs <- matrix(dat[incomplete_obs[[s]][["complete"]], ], ncol=p)
        
        Tls[[s]] <- Tls[[s]] + matrix(rowSums(t(Lambda_s[[s]]) %*% Sig_s1[[s]] %*% (t(observed_vals_from_incomplete_obs) - mu_s[[s]])), ncol=1)

      }
    }

    Txsxs[[s]] <- Txsxs[[s]] / n_s[s]
    Txsfs[[s]] <- Txsfs[[s]] / n_s[s]
    Txsls[[s]] <- Txsls[[s]] / n_s[s]
    Tfsfs[[s]] <- Tfsfs[[s]] / n_s[s]
    Tlsls[[s]] <- Tlsls[[s]] / n_s[s]
    Tfsls[[s]] <- Tfsls[[s]] / n_s[s]

    Txs[[s]] <- Reduce('+', Txis[[s]]) / n_s[s]
    
    # collapse Txis[[s]] into dataframe
    if ("Txis" %in% to_calculate) {
      Txis[[s]] <- Reduce('rbind', Txis[[s]])  # n_s x p
    }

    Tfs[[s]] <- Tfs[[s]] / n_s[s]
    Tls[[s]] <- Tls[[s]] / n_s[s]

  }

  return(list(Txsxs = Txsxs,
              Txsfs = Txsfs,
              Txsls = Txsls,
              Tfsfs = Tfsfs,
              Tlsls = Tlsls,
              Tfsls = Tfsls,
              Txis = Txis,
              Txs = Txs,
              Tfs = Tfs,
              Tls = Tls,
              ds_s = ds_s,
              Sig_s1 = Sig_s1, 
              Sig_s = Sig_s))
}

# log-likelihood
loglik_ecm_msfam <- function(Sig_s1, ds_s, n_s, mu_s, Txis)
{
  S <- length(n_s)
  
  # log-likelihood value for each study
  val_s <- c()
  for(s in 1:S){
    val_s[s] <- -(n_s[s]/2) * log(ds_s[[s]]) -
      (1/2) * sum(sapply(1:n_s[s], function(i) {
        (Txis[[s]][i, ] - mu_s[[s]]) %*% Sig_s1[[s]] %*% (Txis[[s]][i, ] - mu_s[[s]])
      }))
  }
  
  # sum of each study-likelihood
  val_tot <- sum(val_s)
  
  return(val_tot)
}

# ECM algorithm
ecm_msfam <- function(X_s,
                      start,
                      nIt=50000,
                      tol=10^-7,
                      constraint="block_lower2",
                      trace=TRUE) {
  
  S <- length(X_s)
  j_s <- n_s <- numeric(S)
  Phi <- start$Phi
  p <- dim(Phi)[[1]]
  k <- dim(Phi)[[2]]
  Lambda_s <- start$Lambda_s
  psi_s <- start$psi_s

  # dictionary of missingness
  incomplete_cols <- list()
  incomplete_obs <- list()

  # flattened vector of phi, lambda1, ..., psi1, ...
  theta <- param2vect(start, constraint)

  # define objects
  Psi_s1 <- Psi_s <- mu_s <- list()
  L_s <- list()

  # 1st round of cycle
  for(s in 1:S){
    n_s[s] <-  dim(X_s[[s]])[[1]]
    j_s[s] <-  dim(Lambda_s[[s]])[[2]]
    Psi_s[[s]] <- diag(psi_s[[s]])
    Psi_s1[[s]] <-  diag(1/psi_s[[s]])

    mu_s[[s]] <- colMeans(X_s[[s]], na.rm=TRUE)

    # get missing indices
    incomplete_cols[[s]] <- list()
    if (sum(is.na(X_s[[s]])) > 0) {
      incomplete_cols[[s]][["missing"]] <- apply(X_s[[s]], 1, function(ind) which(is.na(ind)))
      incomplete_cols[[s]][["observed"]] <- apply(X_s[[s]], 1, function(ind) which(!is.na(ind)))
    } else {
      incomplete_cols[[s]][["missing"]] <- incomplete_cols[[s]][["observed"]] <- vector(mode="list", length=n_s[s])
    }

    # incomplete_obs[[s]] is a vector containing the observations that are incomplete
    incomplete_obs[[s]] <- list()
    incomplete_obs[[s]][["incomplete"]] <- which(sapply(incomplete_cols[[s]][["missing"]], function(cols) length(cols) > 0))
    incomplete_obs[[s]][["complete"]] <- which(sapply(incomplete_cols[[s]][["missing"]], function(cols) length(cols) == 0))
  }
  
  # E-Step
  out <- exp_values_msfam(Phi, Lambda_s, Psi_s, Psi_s1, mu_s, n_s, X_s,
                          incomplete_cols, incomplete_obs, getdet = FALSE,
                          to_calculate=c("Txis")
  )
  
  Sig_s1 <- out$Sig_s1
  Sig_s <- out$Sig_s
  Txis <- out$Txis
  ds_s <- out$ds_s
  l_stop0 <- 0
  lm1 <- 0
  l0 <- loglik_ecm_msfam(Sig_s1,  ds_s, n_s, mu_s, Txis)
  loglik_chain <- rep(NA, nIt)   # save log-likelihood chain
  
  for (i in (1:nIt))
  {
    ###########CM0 ---------------------------------------------------------------------------------------

    # expected values
    out <- exp_values_msfam(Phi, Lambda_s, Psi_s, Psi_s1, mu_s, n_s, X_s,
                            incomplete_cols, incomplete_obs, getdet = FALSE,
                            to_calculate=c("Txs")
    )
    
    Txs <- out$Txs

    # update mu_s
    for (s in 1:S) {
      mu_s[[s]] <- Txs[[s]]
    }

    ###########CM1 ---------------------------------------------------------------------------------------

    # expected values
    out <- exp_values_msfam(Phi, Lambda_s, Psi_s, Psi_s1, mu_s, n_s, X_s,
                            incomplete_cols, incomplete_obs, getdet = FALSE,
                            to_calculate=c("Txsxs",
                                           "Tfsfs",
                                           "Tlsls",
                                           "Txsfs",
                                           "Txs",
                                           "Txsls",
                                           "Tfs",
                                           "Tls",
                                           "Tfsls")
    )
    
    Txsxs <- out$Txsxs;
    Tfsfs <- out$Tfsfs;
    Tlsls <- out$Tlsls;
    Txsfs <- out$Txsfs;
    Txs <- out$Txs;
    Txsls <- out$Txsls;
    Tfs <- out$Tfs;
    Tls <- out$Tls;
    Tfsls <- out$Tfsls;
    ds_s <- out$ds_s;
    Sig_s1 <- out$Sig_s1;

    # update Psi_s
    Psi_new <- list()
    Psi_new1 <- list()
    psi_new <- list()

    for(s in 1:S){
      psi_new[[s]] <- diag(Txsxs[[s]] +
                             mu_s[[s]] %*% t(mu_s[[s]]) +
                             Phi %*% Tfsfs[[s]] %*% t(Phi) +
                             Lambda_s[[s]] %*% Tlsls[[s]] %*% t(Lambda_s[[s]]) -
                             2 * Txs[[s]] %*% t(mu_s[[s]]) -
                             2 * Txsfs[[s]] %*% t(Phi) -
                             2 * Txsls[[s]] %*% t(Lambda_s[[s]]) +
                             2 * mu_s[[s]] %*% t(Tfs[[s]]) %*% t(Phi) +
                             2 * mu_s[[s]] %*% t(Tls[[s]]) %*% t(Lambda_s[[s]]) +
                             2 * Phi %*% Tfsls[[s]] %*% t(Lambda_s[[s]]))
      Psi_new[[s]] <- diag(psi_new[[s]])
      
      # invert
      Psi_new1[[s]] <- diag(1/diag(Psi_new[[s]]))
    }

    ###########CM2 ---------------------------------------------------------------------------------------

    # expected values
    out <- exp_values_msfam(Phi, Lambda_s, Psi_new, Psi_new1, mu_s, n_s, X_s,
                            incomplete_cols, incomplete_obs, getdet = FALSE,
                            to_calculate=c("Txsfs", "Tfs", "Tfsls", "Tfsfs")
    )
    
    Txsfs <- out$Txsfs;
    Tfs <- out$Tfs;
    Tfsls <- out$Tfsls;
    Tfsfs <- out$Tfsfs;

    # update Phi
    C_s <- list()
    kron_s <- list()
    for(s in 1:S){
      C_s[[s]] <- n_s[s]*Psi_new1[[s]] %*%
        (Txsfs[[s]] - mu_s[[s]]%*%t(Tfs[[s]]) - Lambda_s[[s]]%*%t(Tfsls[[s]]))
      kron_s[[s]] <- kronecker(Tfsfs[[s]], n_s[s]*Psi_new1[[s]])
    }
    
    C <- Reduce('+', C_s)
    kron <- Reduce('+', kron_s)
    Phi_vec <- solve(kron) %*% matrix(as.vector(C))
    Phi_new <- matrix(Phi_vec, p, k)


    ########CM3 ---------------------------------------------------------------------------------------

    # expected values
    out <- exp_values_msfam(Phi_new, Lambda_s, Psi_new, Psi_new1, mu_s, n_s, X_s,
                            incomplete_cols, incomplete_obs, getdet = FALSE,
                            to_calculate=c("Txsls", "Tls", "Tfsls", "Tlsls")
    )
    
    Txsls <- out$Txsls;
    Tlsls <-out$Tlsls;
    Tfsls <- out$Tfsls;
    Tls <- out$Tls;

    # update Lambda_s
    Lambda_new <- list()
    for(s in 1:S){
      Lambda_new[[s]] <- matrix(((Txsls[[s]] - mu_s[[s]] %*% t(Tls[[s]]) - Phi_new %*% Tfsls[[s]]) %*% solve(Tlsls[[s]])), p, j_s[s])
    }

    # apply identifiability constraint

    if (constraint == "null")  {
      Phi_new <- Phi_new
      for (s in 1:S) Lambda_new[[s]] <- Lambda_new[[s]]
    }
    
    if (constraint == "block_lower1")  {
      # see supplemental information in 
      # De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
      # paper for complete data
      # for information on this constraint 
      
      Phi_new[upper.tri(Phi_new)] <- 0
      for (s in 1:S){
        L_s[[s]] <- cbind(Phi_new, Lambda_new[[s]])
        L_s[[s]][upper.tri(L_s[[s]])] <- 0
        Phi_new <- matrix(L_s[[s]][,1:k], nrow=p, ncol = k)
        Lambda_new[[s]] <- L_s[[s]][,(k+1):(k+j_s[s])]
      }
    }

    # the following block ensures the full rank condition holds
    if (constraint == "block_lower2")  {
      lambda_vals <- c()
      psi_vals <- psi_new <- c()
      Phi_new[upper.tri(Phi_new)] <- 0
      phi_val <- as.vector(Phi_new[lower.tri(Phi_new, diag = TRUE)])
      for (s in 1:S){
        Lambda_new[[s]][upper.tri(Lambda_new[[s]])] <- 0
        lambda_vals <- c(lambda_vals, as.vector(Lambda_new[[s]][lower.tri(Lambda_new[[s]], diag = TRUE)]))
        psi_new[[s]] <- diag(Psi_new[[s]])
        psi_vals <- c(psi_vals, psi_new[[s]])
      }
      L_sTOT <- Reduce('cbind', Lambda_new)
      Omega <- cbind(Phi_new, L_sTOT)
      rank_tot <-  qr(Omega)$rank
      theta_new <- c(phi_val, lambda_vals, psi_vals)
      param.struct <- list(Phi = Phi_new, Lambda_s = Lambda_new, psi_s=psi_new)
      Delta <- theta_new - theta
      sh <- 0   # no more than 20 step-halving rounds

      while( (rank_tot < k + sum(j_s)) & (sh<20))
      {
        Delta <- Delta / 2
        sh <- sh + 1
        theta_new <- theta + Delta
        param <- vect2param(theta_new, param.struct, constraint, p, k, j_s)
        Lambda_new <- c()
        psi_new <- param$psi_new
        for(s in 1:S)
        {
          Lambda_new[[s]] <- param$Lambda_s[[s]]
          Psi_new[[s]] <- diag(psi_new[[s]])
          Psi1_new[[s]] <- diag(1 / psi_new[[s]])
        }
        L_sTOT <- Reduce('cbind', Lambda_new)
        Phi_new <- param$Phi
        Omega <- cbind(Phi_new, L_sTOT)
        rank_tot <-  qr(Omega)$rank
      }

      if(sh==20) stop("The full rank condition does not hold\n")
    }

    # stopping rule
    out <- exp_values_msfam(Phi_new, Lambda_new, Psi_new, Psi_new1, mu_s, n_s, X_s,
                            incomplete_cols, incomplete_obs, getdet = FALSE)
    Txis <- out$Txis; 
    Txsxs <- out$Txsxs; 
    Txsfs <- out$Txsfs; 
    Txsls <- out$Txsls; 
    Tfsfs <- out$Tfsfs; 
    Tlsls <- out$Tlsls; 
    Tfsls <- out$Tfsls; 
    Txs <- out$Txs; 
    Tfs <- out$Tfs; 
    Tls <- out$Tls;
    Txis <- out$Txis

    Sig_s1 <- out$Sig_s1
    Sig_s <- out$Sig_s
    ds_s <- out$ds_s
    l1 <- loglik_ecm_msfam(Sig_s1,  ds_s, n_s, mu_s, Txis)
    a <- (l1-l0)/(l0-lm1)
    # ^--ratio of current increase in likelihood to previous increase
    l_stop <- lm1 + (1/ (1-a)) * (l0-lm1)
    # ^--adjust stopping point
    l0 <- l1  # equivalent: loglik_ecm_msfam(Sig_s1,  ds_s, n_s, mu_s, Txis)
    if((trace) & (i %% 1000 == 0))  cat("i=", i, "Criterion for convergence ", abs(l_stop-l_stop0),  "\n")
    if((abs(l_stop-l_stop0)<tol) & i > 1 & l_stop != Inf) break
    Psi_s <- Psi_new
    psi_s <- psi_new
    Phi <- Phi_new
    Lambda_s <- Lambda_new
    if (constraint == "block_lower2") theta <- theta_new
    Psi_s1 <- Psi_new1
    lm1 <- l0
    l0 <- l1
    l_stop0 <- l_stop
    loglik_chain[i] <- l1

    # AIC and BIC computation
    if (constraint == "null") npar <- sum((j_s + k + 2*S) * p)  # total cells across lambda_s, phi, psi_s, mu_s
    if (constraint == "block_lower1") npar <- p * 2*S + k * (p - ( k - 1) / 2) +  sum(j_s * (p - k - (j_s - 1) / 2))
    if (constraint == "block_lower2")  npar <- p * 2*S + k * (p - ( k - 1) / 2) +  sum(j_s * (p  - (j_s - 1) / 2))
    n_tot <- sum(n_s)
    AIC <- -2 * l1 + npar * 2
    BIC <- -2 * l1 + npar * log(n_tot)
  }

  # return output
  res <- list(mu_s = mu_s,
              Phi = Phi, Lambda_s = Lambda_s, psi_s = psi_s, 
              loglik = l1,
              AIC = AIC, BIC = BIC, 
              npar=npar,
              iter = i,
              cov_s = Txis,
              n_s = n_s, 
              constraint = constraint, 
              loglik_chain = loglik_chain)
  return(res)
}

# returns starting values for parameters to intiate ECM algorithm
# this function is adapted from 
# De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
# to allow other factor analysis methods besides 'minres'
start_msfa <- function(X_s, k, j_s, 
                       fa.method = 'minres', 
                       constraint = "block_lower2", 
                       method = "adhoc", 
                       robust = FALSE, 
                       corr = FALSE, 
                       mcd = FALSE)
{
  X_used_s <- X_s
  S <- length(X_s)
  if(corr & !robust)
    for(s in 1:S)  X_used_s[[s]] <- scale(X_s[[s]])
  if(robust & corr & method=="adhoc"){
    for(s in 1:S){
      ogg_s <- if(mcd) covRob(X_s[[s]], estim = "mcd", quan = .75, ntrial = 1000) else covRob(X_s[[s]])
    }
    X_used_s[[s]] <- scale(X_s[[s]], center = ogg_s$center, scale = sqrt(diag(ogg_s$cov)))
  }
  p <- dim(X_s[[1]])[2]
  Phi <- matrix(0, nrow=p, ncol=k)
  Lambda_s <- psi_s <- list()
  if(method == "adhoc"){
    X <- Reduce(rbind, X_used_s)
    X.pcr <- prcomp(X)
    Phi <- matrix(X.pcr$rotation[,1:k], nrow=p, ncol=k, byrow=FALSE)
    
    if (constraint == "block_lower1") {
      Phi[upper.tri(Phi)] <- 0
      for(s in 1:S){
        iniLS <- array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        iniTot <- cbind(Phi, iniLS)
        iniTot[upper.tri(iniTot)] <- 0
        Lambda_s[[s]] <-  matrix(iniTot[,(k+1):(k+j_s[s])], p , j_s[s])
        psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k+j_s[s], fm=fa.method)$uniq
      }
    }
    
    if (constraint == "block_lower2") {
      Phi[upper.tri(Phi)] <- 0
      for(s in 1:S){
        Lambda_s[[s]] = array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        Lambda_s[[s]][upper.tri(Lambda_s[[s]])] <- 0
        psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k+j_s[s], fm=fa.method)$uniq
      }
    }
    
    if (constraint == "null") {
      Phi <- Phi
      for(s in 1:S){
        Lambda_s[[s]] = array(prcomp(X_used_s[[s]])$rotation, dim=c(p, j_s[s]))
        psi_s[[s]] <- fa(X_used_s[[s]], nfactors = k+j_s[s], fm=fa.method)$uniq
      }
    }
  }
  
  # it is important to post-process the output for avoiding sign changes
  if(method == "fa"){
    est <- ecm_fa(X_s, tot_s = k + j_s, robust = robust, mcd = mcd, corr = corr, tol = 10^-5, nIt = 5000, trace = FALSE)
    Phi <- est$Omega_s[[1]][,1:k] / S
    Lambda_s[[1]] <-  est$Omega_s[[1]][,(k+1):(k+j_s[1])]
    psi_s[[1]] <- est$psi_s[[1]]
    for(s in 2:S){
      Phi <- Phi + est$Omega_s[[s]][,1:k] / S * sign(Phi) * sign(est$Omega_s[[s]][,1:k]) ###to avoid sign changes
      Lambda_s[[s]] <-  est$Omega_s[[s]][,(k+1):(k+j_s[s])]
      psi_s[[s]] <- est$psi_s[[s]]
    }
  }
  out <- list(Phi=Phi, Lambda_s=Lambda_s, psi_s=psi_s)
  return(out)
}

# generate incomplete data as an example
# note, we use gene expression counts from microarray of ovarian cancer biopsies
# p=63 cells; n1=285; n2=140
# explained in De Vito, R., Bellio, R., Parmigiani, G. and Trippa, L. (2019). Multi-study Factor Analysis. Biometrics,  75, 337-346.
get_incomplete_data <- function(num_cols_incomplete,  # {1, ...., 63}
                                prop_col_incomplete,  # value between 0 and 1
                                missingness_type)  { # MNAR or MCAR

  dat <- MSFA::data_immune
  S <- length(dat)
  dims <- sapply(1:S, function(s) dim(dat[[s]])); n <- dims[1,]; p <- dims[2,1]
  
  incomplete_cols <- sample(1:p, size=num_cols_incomplete, replace=FALSE)
  
  if (missingness_type == "MCAR") {
    for (s in 1:S) {
      for (col in incomplete_cols) {
        dat[[s]][sample(1:n[s], 
                        size=floor(prop_col_incomplete*n[s]), 
                        replace=FALSE), col] <- NA
      }
    }
  } else if (missingness_type == "MNAR") {
    for (s in 1:S) {
      for (col in incomplete_cols) {
        
        # missingness more likely among low-valued datapoints
        col_values <- dat[[s]][, col]
        col_values_order <- order(col_values, decreasing=TRUE)
        weights <- exp(seq(0, 25, length.out=p)); weights <- weights[col_values_order]
        probs <- weights/sum(weights)
        
        smp <- sample(1:n[s], 
                      size=floor(prop_col_incomplete*n[s], replace=FALSE), 
                      prob=probs)
        
        dat[[s]][smp, col] <- NA
      }
    }
  }
  
  return(dat)
}





