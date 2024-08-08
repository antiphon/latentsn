#' Shotnoise Field Evaluator
#' Create a shotnoise field from given Psi pattern (x,y,dbh,tau)
#' 
#' @param psi Dataframe with coordinates and dbh (in meters) in columns (x,y,d)
#' @param bbox Bounding box, used for creating a regular grid when argument `grid` is missing
#' @param resx Resolution of the regular grid
#' @param cutoff_q Parameter for the kernel
#' @param grid If given, assumed to be a matrix of coordinates (x,y) where to evaluate the field
#' 
#' @details Will compute the shot-noise field at locations (`grid`). 
#' 
#' @returns If grid is given, a vector of values. If grid not given, a spatstat::im object.
#'
#' @import spatstat.geom
#' @export
lsn_shotnoise <- function(psi, 
                          bbox, 
                          resx = 250, 
                          cutoff_q = .99, 
                          grid = NULL){
  #
  grid_given <- !is.null(grid)
  if(!grid_given) {
    resy <- round(  diff(bbox[2,])/diff(bbox[1,]) * resx )
    gx <- seq(bbox[1,1], bbox[1,2], l = resx)# - diff(bbox[1,])/2
    gy <- seq(bbox[2,1], bbox[2,2], l = resy)# - diff(bbox[1,])/2
    grid  <- expand.grid(gx, gy)
  }
  tgrid <- t(grid)
  # brute force.
  if(nrow(psi)>0) {
    D     <- vapply(1:nrow(psi), \(i) 
                    sqrt( (tgrid[1,] -psi$x[i])^2 + (tgrid[2,] - psi$y[i])^2 ) |> 
                      lsn_root_density_cutoff(d = psi$d[i], tau = psi$tau[i], q = cutoff_q), 
                    FUN.VALUE = rep(1.0, nrow(grid)) )
  } else D <- cbind(rep(0, nrow(grid))) # in case empty psi
  
  if(!grid_given) {  
    grid_r <- matrix(rowSums(D), nrow = length(gx))
    Y0     <- im(t(grid_r), gx, gy)
  }
  else {
    Y0 <- rowSums(D)
  }
  Y0
}
