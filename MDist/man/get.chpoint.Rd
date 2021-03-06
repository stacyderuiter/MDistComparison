% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.chpoint.R
\name{get.chpoint}
\alias{get.chpoint}
\title{carry out change-point analysis comparing mdist values to baseline period}
\usage{
get.chpoint(data, quantile_thr = 0.95, thr_type = "both", n = 10000,
  max_gap = 600, exp_times, baseline_start = 0,
  baseline_end = min(exp_times), md_params)
}
\arguments{
\item{data}{the output of the mdist calculation (a data frame with columns t (seconds since start of record) and dist (Mahalanobis distance))}

\item{quantile_thr}{the quantile to be used as a threshold for detection of a change-point. Default is 0.95. This quantile is computed from a set of n maxima from samples from the baseline period that are the same duration as the exposure period.}

\item{thr_type}{which type of threshold to use. options are 'quantile' (quantile_thr quantile of maxima from samples from baseline period that are the same duration as baseline), 'tukey_outlier' (median + 1.5 IQR of the same samples from baseline), or 'both'. Default is 'both'.}

\item{n}{number of iterations for resampling from the baseline period to find threshold}

\item{max_gap}{when determining the start and end times of responses, what is the maximum duration (in seconds) for which mdist can drop below threshold, before one long response is divided into two shorter ones. Default is 10 minutes (600 seconds).}

\item{exp_times}{is a two-element vector of times in seconds indicating the start and end of the experimental period. If baseline_start is not provided, start time of exposure period is used as the end of the baseline period for threshold estimation. Baseline period must be at least as long as the exposure.}

\item{baseline_start}{start time of baseline period for threshold determination. Defaults to start of record.}

\item{baseline_end}{end time of baseline period for threshold determination. Defaults to the start of the experimental period.}

\item{md_params}{a vector or list with information about the settings used for the Mahalanobis distance calculations. Defaults to NULL. This input is just passed unchanged to the result, to facilitate keeping track of experimental period times in simulations.}
}
\value{
a list with entries: resp_cst, resp_cst2, resp_dur=diff(resp_cst), resp_dur2=diff(resp_cst2), resp_ind (the indices in data of the detected response), resp_ind_nogap (first and last samples in data of the response),resp_ind2, resp_ind_nogap2, TH_percentile, TH_outlier, exp_times (passed through from input), mdist=data, md_params (passed through from input)
}
\description{
Do change-point analysis with mdist data. The threshold f
}
