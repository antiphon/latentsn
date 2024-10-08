% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/root-kernel.R
\name{lsn_root_density_cutoff}
\alias{lsn_root_density_cutoff}
\title{Root Kernel With Cutoff}
\usage{
lsn_root_density_cutoff(r, d, tau = 1, q = 0.997)
}
\arguments{
\item{r}{distance argument}

\item{d}{dbh (in meters)}

\item{tau}{shape argument, >0}

\item{q}{cutoff quantile of the kernel}
}
\description{
The kernel with a cutoff for small values.
}
