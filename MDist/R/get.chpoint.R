#' carry out change-point analysis comparing mdist values to baseline period
#'
#' Do change-point analysis with mdist data. The threshold f
#' @param data the output of the mdist calculation (a data frame with columns t (seconds since start of record) and dist (Mahalanobis distance))
#' @param quantile_thr the quantile to be used as a threshold for detection of a change-point. Default is 0.95. This quantile is computed from a set of n maxima from samples from the baseline period that are the same duration as the exposure period.
#' @param thr_type which type of threshold to use. options are 'quantile' (quantile_thr quantile of maxima from samples from baseline period that are the same duration as baseline), 'tukey_outlier' (median + 1.5 IQR of the same samples from baseline), or 'both'. Default is 'both'.
#' @param n number of iterations for resampling from the baseline period to find threshold
#' @param max_gap when determining the start and end times of responses, what is the maximum duration (in seconds) for which mdist can drop below threshold, before one long response is divided into two shorter ones. Default is 10 minutes (600 seconds).
#' @param exp_times is a two-element vector of times in seconds indicating the start and end of the experimental period. If baseline_start is not provided, start time of exposure period is used as the end of the baseline period for threshold estimation. Baseline period must be at least as long as the exposure.
#' @param baseline_start start time of baseline period for threshold determination. Defaults to start of record.
#' @param baseline_end end time of baseline period for threshold determination. Defaults to the start of the experimental period.
#' @param md_params a vector or list with information about the settings used for the Mahalanobis distance calculations. Defaults to NULL. This input is just passed unchanged to the result, to facilitate keeping track of experimental period times in simulations.
#' @return a list with entries: resp_cst, resp_cst2, resp_dur=diff(resp_cst), resp_dur2=diff(resp_cst2), resp_ind (the indices in data of the detected response), resp_ind_nogap (first and last samples in data of the response),resp_ind2, resp_ind_nogap2, TH_percentile, TH_outlier, exp_times (passed through from input), mdist=data, md_params (passed through from input)
#' @importFrom magrittr "%>%"
#' @export

get.chpoint <- function(data,quantile_thr=0.95, thr_type='both', n=10000, 
                        max_gap=600, exp_times, baseline_start=0, 
                        baseline_end=min(exp_times), md_params){
  
  #resample control data to find threshold mdist for change point
  #extract data for times during the baseline period
  pre <- data %>% 
    dplyr::filter(t > baseline_start & t < baseline_end) %>% 
    dplyr::summarise(n()) %>%
    as.numeric()
  # find the duration of the baseline period in units of mdist samples
  edur <- data %>% dplyr::filter(t >= baseline_start & t <= baseline_end ) %>%  length(which(data$t >= exp_times[1] & data$t <= exp_times[2])) #duration of exposure in distance samples
  # find start indices () for n resampled periods the same duration as the real exposure, in samples (row numbers) within "pre"
  S <- sample(1:(length(pre)-edur),n,replace=TRUE) 
  # pre-allocate space for results
  THpre <- matrix(0,nrow=n,ncol=1)
  for (k in c(1:n)){
    THpre[k] <- max(pre[S[k]:(S[k]+edur)],na.rm=TRUE)
  }  
  TH <- stats::quantile(THpre,quantile_thr)
  TH2 <- 1.5*stats::IQR(THpre) + stats::median(THpre)
  
  if (thr_type %in% c('quantile', 'both')){
    #based on resampled periods from the baseline, determine the change threshold
    # find indices of observations where distance exceeds threshold
    ix <- which(data$dist > TH)
    # find indices of times when there is a gap longer than max_gap between above-threshold mdist values
    kk <- which(diff(data$t[ix])>max_gap)
    if (length(kk)>0){ # if there are gaps longer than max_gap
      ix <- ix[1:(kk[1]-1)] #take the indices of the first "change" that exceeds the threshold
    }
    #get the start and end time, in seconds since start of record, of the "response"
    resp_cst <- data$t[c(min(ix),max(ix))]
  }else{
    resp_cst<-NULL
    ix<-0}
  
  if (thr_type %in% c('tukey_outlier', 'both')){
    # second method of threshold estimation -- threshold is median + 1.5 IQR
    #the lines below are same as for the first method above - could make a small function instead but didn't bother
    ixb <- which(data$dist > TH2)
    kkb <- which(diff(data$t[ixb])>max_gap)
    if (length(kkb)>0){
      ixb <- ixb[1:(kkb[1]-1)]
    }
    #get the cst of the start and end times of "response"
    resp_cst2 <- data$t[c(min(ixb),max(ixb))]
  }else{
    resp_cst2<-NULL
    ixb <- 0}
  
  return(list(resp_cst=resp_cst, resp_cst2=resp_cst2,
              resp_dur=diff(resp_cst), resp_dur2=diff(resp_cst2),
              resp_ind=ix, resp_ind_nogap=seq(from=min(ix), to=max(ix)),
              resp_ind2=ixb, resp_ind_nogap2=seq(from=min(ixb), to=max(ixb)),
              TH_percentile=TH, TH_outlier=TH2,
              exp_times=exp_times, mdist=data,
              md_params=md_params))
}