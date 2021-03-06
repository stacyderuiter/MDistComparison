#' calculate Mahalanobis distance for a multivariate time series.
#'
#' Compute Mahalanobis distance comparing a baseline period to observations in a sliding window. For collapsing multivariate time-series prior to change-point analysis.
#' @param data A data frame or matrix with one row for each time point. Note that the Mahalanobis distance calculation should be carried out on continuous data only, so if your data contain logical, factor or character data, proceed at your own risk...errors (or at least meaningless results) will probably ensue.
#' @param fs The sampling rate of data in Hz (data should be regularly sampled). If not specified it will be assumed to be 1 Hz.
#' @param smooth_dur The length, in minutes, of the window to use for calculation of "comparison" values. If not specified or zero, there will be no smoothing (a distance will be calculated for each data observation).
#' @param overlap The amount of overlap, in minutes, between consecutive "comparison"  windows. smooth_dur - overlap will give the time resolution of the resulting distance time series. If not specified or zero,  there will be no overlap.  Overlap will also be set to zero if  smoothDur is unspecified or zero.
#' @param consec Logical. If consec=TRUE, then the calculated distances are between consecutive windows of duration smoothDur, sliding forward over  the data set by a time step of (smoothDur-overlap) minutes.   If TRUE, baselineStart and baselineEnd inputs will be used to define the period used to calculate the data covariance matrix. Default is consec=FALSE.  
#' @param cum_sum Logical.  If cum_sum=TRUE, then output will be the cumulative sum of the calculated distances, rather than the distances themselves.  Default is cum_sum=FALSE.
#' @param baseline_start Start time (in seconds since start of the data set) of the baseline period  (the mean data values for this period will be used as the 'control' to which all   "comparison" data points (or windows) will be compared. if not specified,  it will be assumed to be 0 (start of record).
#' @param baseline_end End time (in seconds since start of the data set) of the baseline period.  If not specified, the entire data set will be used (baseline_end will  be the last sampled time-point in the data set).
#' @param baseline_cov Logical. if TRUE, the variance-covariance matrix for Mahalanobis distance calculation is estimated based on the baseline period only. If FALSE, estimate is based on the entire dataset. Default is FALSE. If using TRUE, be sure you are confident that the baseline period is long enough to provide adequate data for accurate estimation of the matrix.
#' @param parallel       logical.  run in parallel?  NOT IMPLEMENTED YET.  would only help if one figured out how to do rollapply in parallel...
#' @return a data frame with variables t (time in seconds since start of dataset at which distances are reported) and dist (Mahalanobis distances between the specified baseline period and the specified "comparison" periods)
#' @export

mdist <- function(data,fs=1, smooth_dur=0, overlap=0, consec=FALSE, cum_sum=FALSE,  
                  baseline_start=0, baseline_end=floor(nrow(data)/fs), 
                  baseline_cov=FALSE, parallel=FALSE){
  
############################################################
# preliminaries - conversion, preallocate space, etc.
############################################################
bs <- floor(fs*baseline_start) + 1   #start of baseline period in samples
be <- min( ceiling(fs*baseline_end) , nrow(data) ) #end of baseline period in samples
W<-max(1,smooth_dur*fs*60)           #window length in samples
O<-overlap*fs*60                    #overlap between subsequent window, in samples
N<-ceiling(nrow(data)/(W-O))          #number of start points at which to position the window -- start points are W-O samples apart
k <- matrix(c(1:N),ncol=1)          #index vector
ss <- (k-1)*(W-O) + 1               #start times of comparison windows, in samples
ps <- ((k-1)*(W-O) + 1) + smooth_dur*fs*60/2             #mid points of comparison windows, in samples (times at which distances will be reported)
t <- ps/fs                          #mid-point times in seconds
ctr <- colMeans(data[bs:be,], na.rm=T)       #mean values during baseline period
if (baseline_cov){
  bcov <- stats::cov(data[bs:be,], use="complete.obs")           #covariance matrix using all data in baseline period
}else{
  bcov <- stats::cov(data, use="complete.obs")
}
  
############################################################
# Calculate distances!
############################################################
Ma <- function(d,Sx)#to use later...alternate way of calc Mdist
    #d is a row vector of pairwise differences between the things you're comparing
    #Sx is the inverse of the cov matrix
        { sum((d %*% Sx) %*% d)}


if(consec==FALSE){
  #compute rolling means, potentially with overlap, between baseline and individual time-periods
  comps <- zoo::rollapply(data, width = W, mean, by= W-O, by.column=TRUE, align="left",
                          fill=NULL, partial=TRUE, na.rm=T)
  d2 <- apply(comps, MARGIN=1, FUN=stats::mahalanobis, cov=bcov, center=ctr, inverted=FALSE)
  }else{# if consec = TRUE, compute mdists between consecutive windows
    i_bcov <- solve(bcov) #inverse of the baseline cov matrix
    #rolling means, potentially with overlap
    ctls <- zoo::rollapply(data, width = W, mean, by= W-O, by.column=TRUE, align="left", 
                           fill=NULL, partial=TRUE, na.rm=T)
    #compare a given control window with the following comparison window.
    comps <- rbind( ctls[2:nrow(ctls),] , NA*vector(mode="numeric",length=ncol(data)) ) 
    pair_diffs <- as.matrix(ctls-comps)
    d2 <- apply(pair_diffs, MARGIN=1, FUN=Ma, Sx=i_bcov)
    #first dist should be at midpoint of first comp window
    d2 <- c(NA, d2[1:(length(d2)-1)]) 
  }

#functions return squared Mahalanobis dist so take sqrt
dist<-sqrt(d2) 

#erase the values for partial windows and replace with NAs.  
#because the distances get bigger for partial windows, not b/c of change, but because of less averaging...
dist[t > (nrow(data)/fs - smooth_dur*60)] <- NA

# Calculate cumsum of distances if requested
# this is kind of silly. maybe it'll be more use having this in here if we decide to calculate 
# the cumsum after a specified start time, e.g. from start of exposure...
# or maybe just better to do later in plotting routines.

if(cum_sum==TRUE){
  dist<-cumsum(dist)
  }

#Ta-Da!
D <- data.frame(t,dist)
return(D)
}
