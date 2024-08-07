#' Root Kernel R-part
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
#' @export
lsn_root_density_d_part <- function(d) {
  #pi * exp(-14.63682 + 3.034176 * log(2 + 125*d))
  1.381846 * (2+125*d)^3.034 / 1e6
}

#' Root Kernel With Cutoff
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
#' @export
lsn_root_density <- function(r, d, tau = 1) {
  lsn_root_density_r_part(r, tau) * root_density_d_part(d)
}