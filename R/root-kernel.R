#' Root Kernel R-part
#' 
#' The range-dependent part of the kernel.
#' 
#' @param r distance argument
#' @param tau shape argument, >0
#' 
#' @export
lsn_root_density_r_part <- function(r, tau = 1) {
  # Mass preserving modification. the constant is 2/c.
  # if tau is massive, cut to a large number
  c_tau <- pmin(.Machine$double.xmax, (gamma(1/tau) * 1.435956 / tau )^tau)
  #exp(-r^tau * c_tau)  
  exp(-exp(tau * log(r)) * c_tau)
}

#' Root Kernel D-part
#' 
#' The size-dependent part of the kernel.
#' 
#' @param d dbh, in meters 
#' @export
lsn_root_density_d_part <- function(d) {
  #pi * exp(-14.63682 + 3.034176 * log(2 + 125*d))
  1.381846 * (2+125*d)^3.034 / 1e6
}

#' Root Kernel With Cutoff
#' 
#' The kernel with a cutoff for small values.
#' 
#' @param r distance argument
#' @param d dbh (in meters) 
#' @param tau shape argument, >0
#' @param q cutoff quantile of the kernel
#'
#' @export
lsn_root_density_cutoff <- function(r, d, tau = 1, q = 0.997) {
  e <- lsn_root_density_r_part(r, tau) * lsn_root_density_d_part(d)
  # cutoff, based on kernel shape only.
  if(length(e)) {
    qn <- (-log(1-q))^(1/tau) / (1.435956 * gamma(1/tau)/tau)
    e[ r > qn] <- 0
  }
  e
}


#' Root Kernel
#' 
#' Evaluate the kernel.
#' 
#' @param r distance argument
#' @param d dbh (in meters) 
#' @param tau shape argument, >0
#' 
#' @export
lsn_root_density <- function(r, d, tau = 1) {
  lsn_root_density_r_part(r, tau) * root_density_d_part(d)
}