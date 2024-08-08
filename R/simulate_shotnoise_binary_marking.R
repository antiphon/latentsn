#' Default Simpars
#' 
#' @export
lsn_simpars_default <- function(bbox = rbind(0:1, 0:1) * 100,
                                Rmax = 10,
                                ndata = 100,
                                # shotnoise kernel
                                tau = 1, 
                                cutoff_q = 0.99,
                                # parent process parameters
                                psi_lambda0 = 0.01,
                                psi_gamma = 0, # strauss repulsion strength
                                psi_range = 5, # strauss range
                                psi_dbh = c(50, 200), # dbh in meters
                                # 
                                # logistic regression part
                                theta = c(qlogis(.15), log(2)),
                                rho_scale = 1000,
                                ...
                                ) {
  
  #formals()
  c(as.list(environment(), list(...)) )
}



#' Simulate Latent Shotnoise Dependent Binary Marking
#' 
#' @param simpars list of simulation parameters. See lsn_simpars_default
#' @param xy fixed locations of the trees to mark, ndata x 2 matrix of coordinates 
#' @param ... passed on to lsn_simpars_default
#' 
#' @details
#' Todo
#' @returns list(phi=tibble of data points with marks, psi = latent generator points, simpars = simpars)
#'
#' @import dplyr
#' @export
lsn_simulate <- function(simpars = lsn_simpars_default(), 
                         xy, # optional locations to mark 
                         psi, 
                         ...) {
  # shape of kernel
  tau <-  simpars$tau
  bbox <- simpars$bbox
  Rmax <- simpars$Rmax
  #gamma <- simpars$psi_gamma
  #range <- simpars$psi_range
  cutoff_q <- simpars$cutoffq
  #if(is.null(gamma)) gamma <- 0
  #if(is.null(range)) range <- 0
  if(is.null(Rmax)) Rmax <- 0
  if(is.null(cutoff_q)) cutoff_q <- 0.99
  
  # the generator point pattern
  if(!missing(psi)){
    stopifnot(is.data.frame(psi))
    stopifnot( all(c("x", "y", "d") %in% names(psi)) )
  }
  else{
    bbox_sim <- t(bbox) + c(-1,1) * Rmax
    psi_n <- rpois(1, simpars$psi_lambda0 * diff(bbox_sim[,1]) * diff(bbox_sim[,2]))
    psi   <- tibble(x = runif(psi_n, bbox_sim[1,1], bbox_sim[1,2]),
                    y = runif(psi_n, bbox_sim[2,1], bbox_sim[2,2]),
                    d = rgamma(psi_n, simpars$psi_dbh[1], simpars$psi_dbh[2]))
  }
  
  psi_x   <- psi[,c("x", "y")] |> as.matrix()
  psi_dbh <- psi[["d"]]
  psi_n <- nrow(psi_x)
  psi   <- tibble(
    x = psi_x[,1],
    y = psi_x[,2],
    d = psi_dbh) |> mutate(    tree = row_number(),
                               tau = tau # kernel shape
    ) # in meters
  # then shotnoise psi
  #U    <- shotnoise(psi, bbox, cutoff_q = cutoff_q)
  #
  # Then sample phi
  if(missing(xy)) {
    phi_n <- simpars$ndata
    # locations
    xy <- cbind(runif(phi_n, bbox[1,1], bbox[1,2]), runif(phi_n, bbox[2,1], bbox[2,2]))
  }
  phi_n <- nrow(xy)
  # compute exactly, like in the fitter
  D     <- apply(psi[,1:2] |> as.matrix(), 1, \(x) sqrt( colSums( (t(xy)-x)^2 ) ) ) |> t()
  u    <- colSums(  lsn_root_density_cutoff(D, psi$d, tau = tau, q = cutoff_q) )
  if(length(u) != phi_n) u <- 0
  # compile
  
  phi  <- tibble(
    tree = 1:phi_n,
    x = xy[,1],
    y = xy[,2],
    # And mark
    U =  u,
    # Coin flip
    xii = rbinom(phi_n, 1, 1/(1 + exp(-simpars$theta[1] - simpars$theta[2] * simpars$rho_scale * U))),
    infected = factor(xii, levels = 0:1, labels = c("no rot", "yes rot") )
  )
  # done.
  list(phi = phi,
       psi = psi,
       simpars = simpars)
}
