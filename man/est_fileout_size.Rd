% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{est_fileout_size}
\alias{est_fileout_size}
\title{Auxiliary for estimating the size of the written netcdf file}
\usage{
est_fileout_size(nloc, ntime, nvars)
}
\arguments{
\item{nloc}{number of loaded gridcells}

\item{ntime}{number of loaded time-slizes}

\item{nvars}{number of loaded variables}
}
\description{
Returns an estimate in Megabytes.
}
