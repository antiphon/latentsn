#' MCMC Estimate Shotnoise Kernel
#'
#' @param phi the data, tibble(x, y, xii) where the marks xii = 0,1 
#' @param prior_dbh generator dbh
#' @param prior_theta logistic regression coefficients
#' @param prior_log_tau shape parameter, normal
#' @param prior_lambda0 prior for intensity of generators (or first parameter if using Strauss)
#' @param log_tau optional initial value
#' @param theta optional initial value
#' @param lambda0 optional initial value
#' @param rho_scale adjustment parameter. Default 1000.
#' @param niter number of mcmc steps to run
#' @param print_every print freq
#' @param keep_every store generators only at this step iteration
#' @param loc_storage pre-sample this many uniform locations (faster to batch for polygon win)
#' @param win optional observation window, spatstat.geom::owin
#' @param ada_bi adaptive mcmc starts from this many interations (default: niter/2)
#' @param prior_psi strauss prior parameters c(gamma, range)
#' @param Rmax buffering / edge correction
#' @param cutoff_q kernel parameter
#' @param ... 
#' @param timeout max secs to use
#' 
#' @details all linear units are assumed to be the same, i.e. same units for locations and dbh and range-scales. 
#' Any covariates can be given as additional columns in phi, names starting with "X_", e.g. "X_elevation". Note that
#' prior_theta should have the corresponding prior added to between first (intercept) and the last (reserved for random field).
#' 
#' @import spatstat.geom mvnfast
#' @export
lsn_mcmc_estimate <- function(phi,             # data point pattern coordinate tibble(x,y,xii), xii = 0,1 infection status
                              prior_lambda0 = c(1, 1000), # Gamma(shape, rate)
                              prior_dbh     = c(1,3), # ~gamma. 
                              prior_theta   = NULL, # Depends on covariate amount 
                              prior_log_tau = c(0, 1), # ~ normal
                              lambda0 = NULL,
                              theta   = NULL,
                              log_tau = NULL, # initial values
                              rho_scale = 1000, # arbitrary scaler
                              niter = 1000,
                              print_every = -1, # when to print progress. if < 0, every 10%
                              keep_every  = 100,   # for storing psi
                              loc_storage = 10000, # for precaching sampling of new locations
                              win = NULL,
                              ada_bi = niter/2,
                              prior_psi = c(gamma = 0, range = 0), # 0 = poisson, < 0 = repulsion in psi
                              Rmax = NULL, # expansion for winodow of psi
                              cutoff_q = 0.997, # kernel cutoff 
                              ..., #ignored,
                              timeout = Inf # seconds to use max.
) {
  t0 <- Sys.time()
  
  # check time. Return TRUE if taken too long
  if(timeout < Inf) {
    timeout_check <- function() {
      secs <- as.numeric(Sys.time()-t0, "secs")
      if(timeout <= secs) stop(sprintf("Taking too long (>= %is)", round(timeout)))
    }
  }
  else timeout_check <- function() NULL
  
  
  # expanded window:
  bbox <- cbind(range(phi$x), range(phi$y))
  if(is.null(Rmax)) {
    # meaningful range?
    rvec <- seq(0, bbox[2,1]-bbox[1,1], l = 100)
    denv <- lsn_root_density_r_part(rvec, tau = 1)
    Rmax <- sum( denv > (0.01 * max(denv))) # 1% of maximum
  }
  win_ex  <- cbind(range(phi$x), range(phi$y)) + c(-1,1) * Rmax 
  bbox_ex <- win_ex
  vol     <- apply(win_ex, 2, diff) |> prod()
  # sample one point from the region.
  sample_loc1 <- function(...) c(runif(1, win_ex[1,1], win_ex[2,1]), runif(1, win_ex[1,2], win_ex[2,2]))
  # pre-simulate new locations, might be slow e.g. polygons
  if( is.null(win)) ready_locs <- sapply(1:loc_storage, sample_loc1) |> t()
  
  # if we have spatstat::owin
  if(!is.null(win)) {
    win_ex         <- spatstat.geom::dilation(win, Rmax)
    vol        <- spatstat.geom::area(win_ex)
    bbox_ex    <- spatstat.geom::boundingbox(win_ex) |> with(cbind(xrange, yrange))
    sample_loc1 <- function(..., n = 1) with(spatstat.random::runifpoint(n, win = win_ex)  , c(x,y))
    # pre-simulate new locations, might be slow e.g. polygons
    ready_locs <- sample_loc1(n=loc_storage) |> matrix(ncol = 2)
  }
  # the sampler of new locations.
  locs_next <- 1
  sample_loc <- function(...) {
    xout      <- ready_locs[locs_next, ]
    locs_next <<- locs_next + 1
    # create new set
    if(locs_next > loc_storage) {
      if(is.null(win)) ready_locs <<- sapply(1:loc_storage, sample_loc1) |> t()
      else ready_locs <<- sample_loc1(n=loc_storage) |> matrix(ncol = 2)
      locs_next <<- 1
    }
    xout
  }
  # prepare data location parts
  phi_xyt <- rbind(phi$x, phi$y)
  xii     <- phi$xii
  n       <- nrow(phi)
  #regression covariate structure.
  covnames <- grep("^X_", colnames(phi), value = TRUE)
  covX    <- phi[, covnames] |> as.matrix()
  fixefs <- c("Intercept", covnames)
  
  # Add intercept, and a placeholder for Us
  #covX    <- cbind(1, covX, 0) 
  covX    <- cbind(1, covX)  # omit U here for some speed gain
  n_theta <- ncol(covX) + 1 
  #
  # initial values 
  if(is.null(lambda0)) lambda0 <- rgamma(1, prior_lambda0[1], prior_lambda0[2])
  n0   <- max(2, round(lambda0 * vol)) # lazy here.
  # Uniform starting locations.
  psi  <- apply(t(1:n0), 2, sample_loc) |> t()
  # add initial dbh's 
  psi  <- cbind(psi, rgamma(n0, prior_dbh[1], prior_dbh[2]))
  
  # initial
  # other pars
  if(is.null(prior_theta)) {
    prior_theta <- list(m = rep(0, n_theta), S = diag(1e5, n_theta))
  }
  else if( length(prior_theta$m) != n_theta || ncol(prior_theta$S) != n_theta) {
    stop(sprintf("prior_theta should be a list with m = vect(%i), S=matrix(%i,%i)", n_theta, n_theta, n_theta))
  }
  if(is.null(theta)) theta     <- mvnfast::rmvn(1, prior_theta$m, prior_theta$S)
  if(is.null(log_tau)) log_tau <- rnorm(1, prior_log_tau[1], prior_log_tau[2])
  # The joint probability at the current state. For ratios
  # 
  # Prior functions
  l_prior_gamma    <- function(x, p)  p[1]*log(p[2]) - lgamma(p[1]) + log(x)*(p[1]-1) - p[2]*x
  l_prior_norm     <- function(x, p) -0.9189385 - log(p[2]) - 0.5 * (x - p[1])^2/p[2]^2 
  l_prior_mvnorm   <- function(x, p) mvnfast::dmvn(x, p$m, p$S, log = TRUE)
  l_prior_lnorm    <- function(x, p) l_prior_norm(x, p) - x
  #l_prior_lpoisson      <- function(x, p) x * log(p) - p - lfactorial(x) 
  #
  # Prior for the latent point location process.
  # Poisson 
  l_prior_psi_poi         <- function(x, p, ...) nrow(x) * log(p) + (1-p) * vol     
  # Strauss: lambda is the beta now. (gamma, R) fixed.
  # From approximate_strauss_constant_R:
  strauss_r2 <- prior_psi[['range']]^2
  dstrauss_constant_a <- -.5 * (1-exp(prior_psi[['gamma']])) * pi * strauss_r2/vol
  
  l_prior_psi_strauss <- function(x, p, npairs) {
    if(missing(npairs)) npairs <- sum(dist(x[,1:2,drop=FALSE]) < prior_psi[['range']])
    bee <- vol * p
    cee <- 0.5/bee - dstrauss_constant_a
    dee <- 0.5*(dstrauss_constant_a - 1) / cee
    ee  <- 0.5/cee
    lZ  <- -vol + bee + 0.5*log(ee/bee) - 0.5*bee + cee*dee^2
    nrow(x)*log(p) + npairs * prior_psi[['gamma']] - lZ
  }
  # choose, gamma=0 = poisson faster
  l_prior_psi <- if(prior_psi['gamma'] == 0) l_prior_psi_poi else l_prior_psi_strauss
  
  
  # likelihood 
  # For the observed binary infections.
  xii1 <- which(xii == 1)
  xii0 <- which(xii == 0)
  eps_low <- .Machine$double.eps # to avoid singularities
  eps_high <- 1 - eps_low
  l_lik_phi <- function(theta, Us) {
    # scale rho:
    eta    <- c(covX %*% theta[-n_theta] + Us * rho_scale * exp(theta[n_theta])) #
    p_xi_i <- 1/(1 + exp(-eta)) 
    p_xi_i[p_xi_i < eps_low] <-  eps_low
    p_xi_i[p_xi_i > eps_high] <- eps_high # pmax(eps_low, pmin(eps_high, p_xi_i) )
    sum( log(p_xi_i[xii1]) ) + sum( log(1-p_xi_i[xii0]))
  }
  # likelihood, full.
  lhood <- function(psi, lambda0, Us, log_tau, theta, npairs){
    # phi likelihood
    p_obs           <- l_lik_phi(theta, Us) 
    # Priors
    # Latent pattern: coords and sizes independent
    p_psi           <- l_prior_psi(psi, lambda0, npairs)
    p_prior_dbh     <- sum( l_prior_gamma(psi[,3], prior_dbh) )
    # Regression coefficients
    p_prior_theta      <- l_prior_mvnorm(theta, prior_theta) 
    # Hyper
    p_prior_lambda0 <- l_prior_gamma(lambda0, prior_lambda0)
    p_prior_log_tau <- l_prior_norm(log_tau, prior_log_tau) 
    # total
    p_obs + p_psi + p_prior_theta + p_prior_lambda0 + p_prior_log_tau +  
      p_prior_dbh 
  }
  # likelihood, as vector of parts
  lhood_parts <- function(psi, lambda0, Us, log_tau, theta, npairs) {
    c(
      obs           = l_lik_phi(theta, Us), 
      # Priors
      # Latent pattern: coords and sizes independent
      psi           = l_prior_psi(psi, lambda0, npairs),
      prior_dbh     = sum( l_prior_gamma(psi[,3], prior_dbh) ),
      # Regression coefficients
      prior_theta   =  l_prior_mvnorm(theta, prior_theta), 
      # Hyper
      prior_lambda0 = l_prior_gamma(lambda0, prior_lambda0),
      prior_log_tau = l_prior_norm(log_tau, prior_log_tau) 
    )
  }
  #
  # the interaction components.
  D     <- apply(psi[,1:2], 1, \(x) sqrt( colSums( (phi_xyt-x)^2 ) ) ) |> t()
  Us    <- colSums(  lsn_root_density_cutoff(D, psi[,3], tau = exp(log_tau), q = cutoff_q) )
  #
  lp_old <- lhood_parts(psi, lambda0, Us, log_tau, theta)
  # 
  # storage
  hist_n0      <- rep(0L, niter)
  psi_to_keep  <- ceiling(niter/keep_every)
  hist_psi     <- vector("list", psi_to_keep)
  hist_theta   <- array(0, dim = c(niter, n_theta) )
  hist_log_tau <- rep(0, niter)
  hist_lambda0 <- rep(0, niter)
  hist_acc     <- rep(0, niter)
  hist_psi_it  <- rep(0, psi_to_keep)
  acc <- c(b=0, d = 0, m = 0, resize = 0, theta = 0, lambda0=0, log_tau = 0)
  hpk <- 0 # counter for intermediate histories, e.g psi
  #
  # Proposal sd's when sampling from normal. adapt from this.
  proposal_qs <- c(lambda0 = 0.05, log_tau = 0.05)
  proposal_L_theta <- diag(0.1, n_theta) # Chol lower Tri of Proposal Sigma
  # Adaptive, bookkeeping
  ada_prop    <- 0 * c(proposal_qs , theta = 0)
  ada_rat     <- 0 * ada_prop
  #
  if(print_every < 0) print_every <- round(0.1 * niter, -1)
  #
  ch <- acc[1:3] # for the bdm part
  for(it in 1:niter) {
    # Modify locations:
    # Either move, or (b,d)
    p_bd <- 0.75 #  
    # what ?
    # birth or death
    if(runif(1) < p_bd ) {
      p_b  <- 0.5 # if empty, nothing to kill.
      if( runif(1)  <  p_b) { # birth
        ci      <- 1
        xnew    <- sample_loc()
        dbh_new <- rgamma(1, prior_dbh[1], prior_dbh[2]) 
        dnew    <- sqrt( colSums( (phi_xyt - xnew)^2 ) ) 
        Usnew   <- Us + lsn_root_density_cutoff(dnew, dbh_new, exp(log_tau), q = cutoff_q)
        psi_new <- rbind(psi, c(xnew, dbh_new) )
        lp_new        <- lp_old
        lp_new['obs'] <- l_lik_phi(theta, Usnew)
        lp_new['psi'] <- l_prior_psi(psi_new, lambda0)
        lp_new['prior_dbh'] <- sum( l_prior_gamma(psi_new[,3], prior_dbh) )
        q_prop        <-  (1-p_b)/p_b * (1 / (n0+1)) * (1/ (dgamma(dbh_new, prior_dbh[1], prior_dbh[2])/vol) )
        rat           <- exp( sum(lp_new - lp_old) ) * q_prop
        if(is.na(rat)) browser();         if(runif(1) < rat) {
          Us    <- Usnew
          psi   <- psi_new
          lp_old <- lp_new
          D     <- rbind(D, dnew)
          acc[1] <- acc[1] + 1
        }
      }
      # death
      else{
        ci      <- 2
        if(n0 > 0) {
          j       <- sample(1:n0, 1)
          psi_new <- psi[-j,,drop=FALSE]
          kold    <- lsn_root_density_r_part(D[j,], exp(log_tau))
          Usnew   <- Us - lsn_root_density_cutoff(D[j,], psi[j, 3], exp(log_tau), q = cutoff_q)
          lp_new        <- lp_old
          lp_new['obs'] <- l_lik_phi(theta, Usnew)
          lp_new['psi'] <- l_prior_psi(psi_new, lambda0)
          lp_new['prior_dbh'] <- sum( l_prior_gamma(psi_new[,3], prior_dbh) )
          q_prop        <- p_b/(1-p_b) * n0 / vol * (dgamma(psi[j,3], prior_dbh[1], prior_dbh[2]))
          rat           <- exp( sum(lp_new - lp_old) ) * q_prop
          if(runif(1) < rat) {
            Us    <- Usnew
            lp_old <- lp_new
            psi   <- psi_new
            D     <- D[-j,, drop=FALSE]
            acc[2] <- acc[2] + 1
          }
        }
      }
    }
    # move
    else{
      ci      <- 3
      if(n0>0){ # if nothing to move, do nothing
        j       <- sample(1:n0, 1)
        xnew    <- sample_loc()
        xold    <- psi[j, 1:2]
        psi[j, 1:2] <- xnew
        # strauss
        #browser()
        dnew    <- sqrt( colSums( (phi_xyt - xnew)^2 ) ) 
        Usnew   <- Us + (lsn_root_density_cutoff(dnew, psi[j, 3], exp(log_tau), q = cutoff_q) - 
                           lsn_root_density_cutoff(D[j,], psi[j, 3], exp(log_tau),  q = cutoff_q))
        lp_new        <- lp_old
        lp_new['obs'] <- l_lik_phi(theta, Usnew)
        lp_new['psi'] <- l_prior_psi(psi, lambda0)
        rat           <- exp( sum(lp_new - lp_old) )
        if(runif(1) < rat) {
          Us     <- Usnew
          lp_old  <- lp_new
          acc[3] <- acc[3] + 1
          D[j,]  <- dnew
        }else{psi[j,1:2] <- xold}
      }
    }
    npairs <- sum(dist(psi[,1:2,drop=FALSE]) < prior_psi['range'])
    ch[ci] <- ch[ci] + 1
    n0 <- nrow(psi)
    #
    # Update dbh's: log-gaussian scaling
    # Sweep? or Just some.
    if(n0 > 0){
      j       <- sample(n0, min(n0, 3))
      l_f     <- exp( rnorm(length(j), 0, .02) )
      dbh_old <- psi[j, 3]
      dbh_new <- dbh_old * l_f
      Kj       <- lsn_root_density_r_part(D[j,,drop=FALSE], exp(log_tau) )
      Usnew    <- Us + colSums(lsn_root_density_cutoff(D[j,,drop=FALSE], dbh_new, exp(log_tau),  q = cutoff_q) -
                                 lsn_root_density_cutoff(D[j,,drop=FALSE], dbh_old, exp(log_tau),  q = cutoff_q) )
      psi[j,3] <- dbh_new
      lp_new              <- lp_old
      lp_new['obs']       <- l_lik_phi(theta, Usnew)
      lp_new['prior_dbh'] <- sum( l_prior_gamma(psi[,3], prior_dbh) )
      rat                 <- exp( sum(lp_new - lp_old) ) * prod(l_f)
      if(runif(1) < rat) {
        lp_old <- lp_new
        Us    <- Usnew
        acc['resize'] <- acc['resize'] + 1
      }else {psi[j, 3] <- dbh_old}
    }
    
    # update other parameters:
    
    # theta symmetric proposition
    u_theta   <- rnorm(length(theta), 0, 1)
    theta_new <- theta + u_theta %*% proposal_L_theta
    lp_new        <- lp_old
    lp_new['obs'] <- l_lik_phi(theta_new, Us)
    lp_new['prior_theta'] <- l_prior_mvnorm(theta_new, prior_theta) 
    
    rat    <- exp( sum(lp_new - lp_old) )
    if(runif(1) < rat) {
      lp_old  <- lp_new
      theta  <- theta_new
      acc['theta'] <- acc['theta'] + 1
    }
    ada_rat['theta'] <- rat
    
    # lambda0: 
    # symmetric
    l_f      <- exp(rnorm(1, 0, proposal_qs['lambda0']))
    lambda0_new <- lambda0 * l_f
    lp_new                  <- lp_old
    lp_new['psi']           <- l_prior_psi(psi, lambda0_new, npairs)
    lp_new['prior_lambda0'] <- l_prior_gamma(lambda0_new, prior_lambda0)
    #
    rat    <- exp( sum(lp_new - lp_old) ) * l_f # log-gaussian proposal
    if(runif(1) < rat) {
      lp_old      <- lp_new
      lambda0    <- lambda0_new
      acc['lambda0'] <- acc['lambda0'] + 1
    }
    ada_rat['lambda0'] <- rat
    
    
    # Kernel parameters:
    l_f     <-  rnorm(1, 0, proposal_qs['log_tau'])
    log_tau_new <- log_tau + l_f
    # need to recompute all weights.
    Usnew  <- colSums( lsn_root_density_cutoff(D, psi[,3], exp(log_tau_new),  q = cutoff_q) )
    lp_new                  <- lp_old
    lp_new['obs']           <- l_lik_phi(theta, Usnew)
    lp_new['prior_log_tau'] <- l_prior_norm(log_tau_new, prior_log_tau)
    #
    rat    <- exp( sum(lp_new - lp_old) )
    if(runif(1) < rat) {
      Us         <- Usnew
      lp_old      <- lp_new
      log_tau    <- log_tau_new
      acc['log_tau'] <- acc['log_tau'] + 1
    }
    ada_rat['log_tau'] <- rat
    
    
    # remember
    hist_n0[it] <- n0
    hist_theta[it,] <- theta
    hist_lambda0[it] <- lambda0
    hist_acc[it] <- mean(acc/it)
    hist_log_tau[it] <- log_tau
    
    # adjust proposal variances
    if(ada_bi > it) {
      # # # just these for now
      ada_i <- 1:2
      # use the simplfilication for 1d no interaction
      change <- pmin(1, ada_rat[ada_i]) - 0.234
      uu     <- proposal_qs * sqrt(it^(-2/3) * abs(change))
      S      <- sqrt(proposal_qs^2 + sign(change) * uu^2)
      proposal_qs[ada_i] <- c(S)
      # then thetas
      proposal_L_theta  <- ramcmc::adapt_S( S = proposal_L_theta,
                                            u = u_theta,
                                            current = min(1, ada_rat['theta']),
                                            n = it)
      
      
      #if(it > 1000) browser()
    }
    #  if(it > 0) browser()
    #
    # psi storage
    if(it %% keep_every == 0) {
      hist_psi[[hpk <- hpk + 1]] <- psi
      hist_psi_it[hpk] <- it
    }
    
    if(it %% print_every == 0) 
      cat( sprintf("%i/%i (bdm:[%i/%i %i/%i %i/%i])(n0: %i) \n", 
                   it, niter, 
                   acc[1], ch[1],
                   acc[2], ch[2],
                   acc[3], ch[3], 
                   n0) )
    # timeout?
    timeout_check()
  } # eo MCMC loop
  
  colnames(hist_theta) <- c(fixefs, "U") 
  
  # store also parameter
  cfg <- list(priors = list(lambda0 = prior_lambda0,
                            theta   = prior_theta,
                            log_tau = prior_log_tau,
                            dbh = prior_dbh), 
              hyper = list(psi = prior_psi),
              niter = niter,
              keep_every = keep_every,
              rho_scale = rho_scale,
              Rmax = Rmax,
              cutoff_q = cutoff_q
  )
  # acceptance rates
  ch <- c(ch, acc[-(1:3)] * 0  + niter )
  # set nicer names
  #print(c(is=is(hist_theta), n=ncol(hist_theta)))
  
  # compile 
  list (rates  = cbind(tries = ch,
                       accept = acc,
                       prop   = acc/ch),
        last_psi = psi,
        hist_psi = hist_psi, 
        hist_psi_it = hist_psi_it,
        hist_par = data.frame(n0      = hist_n0,
                              lambda0 = hist_lambda0,
                              log_tau = hist_log_tau,
                              hist_theta),
        hist_acc = hist_acc,
        fixefs = fixefs,
        ranefs = "U",
        vol = vol,
        win_sim = win_ex,
        win_obs = bbox, # empirical bbox
        rho_scale = rho_scale,
        took = Sys.time() - t0,
        config = cfg,
        fitter_ver = "11.0",
        kernel_ver = "2.0"
  )
}
