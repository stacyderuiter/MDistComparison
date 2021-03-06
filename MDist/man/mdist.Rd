% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdist.R
\name{mdist}
\alias{mdist}
\title{calculate Mahalanobis distance for a multivariate time series.}
\usage{
mdist(data, fs = 1, smooth_dur = 0, overlap = 0, consec = FALSE,
  cum_sum = FALSE, baseline_start = 0,
  baseline_end = floor(nrow(data)/fs), baseline_cov = FALSE,
  parallel = FALSE)
}
\arguments{
\item{data}{A data frame or matrix with one row for each time point. Note that the Mahalanobis distance calculation should be carried out on continuous data only, so if your data contain logical, factor or character data, proceed at your own risk...errors (or at least meaningless results) will probably ensue.}

\item{fs}{The sampling rate of data in Hz (data should be regularly sampled). If not specified it will be assumed to be 1 Hz.}

\item{smooth_dur}{The length, in minutes, of the window to use for calculation of "comparison" values. If not specified or zero, there will be no smoothing (a distance will be calculated for each data observation).}

\item{overlap}{The amount of overlap, in minutes, between consecutive "comparison"  windows. smooth_dur - overlap will give the time resolution of the resulting distance time series. If not specified or zero,  there will be no overlap.  Overlap will also be set to zero if  smoothDur is unspecified or zero.}

\item{consec}{Logical. If consec=TRUE, then the calculated distances are between consecutive windows of duration smoothDur, sliding forward over  the data set by a time step of (smoothDur-overlap) minutes.   If TRUE, baselineStart and baselineEnd inputs will be used to define the period used to calculate the data covariance matrix. Default is consec=FALSE.}

\item{cum_sum}{Logical.  If cum_sum=TRUE, then output will be the cumulative sum of the calculated distances, rather than the distances themselves.  Default is cum_sum=FALSE.}

\item{baseline_start}{Start time (in seconds since start of the data set) of the baseline period  (the mean data values for this period will be used as the 'control' to which all   "comparison" data points (or windows) will be compared. if not specified,  it will be assumed to be 0 (start of record).}

\item{baseline_end}{End time (in seconds since start of the data set) of the baseline period.  If not specified, the entire data set will be used (baseline_end will  be the last sampled time-point in the data set).}

\item{baseline_cov}{Logical. if TRUE, the variance-covariance matrix for Mahalanobis distance calculation is estimated based on the baseline period only. If FALSE, estimate is based on the entire dataset. Default is FALSE. If using TRUE, be sure you are confident that the baseline period is long enough to provide adequate data for accurate estimation of the matrix.}

\item{parallel}{logical.  run in parallel?  NOT IMPLEMENTED YET.  would only help if one figured out how to do rollapply in parallel...}
}
\value{
a data frame with variables t (time in seconds since start of dataset at which distances are reported) and dist (Mahalanobis distances between the specified baseline period and the specified "comparison" periods)
}
\description{
Compute Mahalanobis distance comparing a baseline period to observations in a sliding window. For collapsing multivariate time-series prior to change-point analysis.
}
