get.chpoint <- function(mdo,quantile.thr, CEEtimes, md.params){
  #mdo is the output of the mdist calculation,
  # data frame with columns t and dist
  # quantile.thr is threshold from baseline to use for chpoint.
  #CEEtimes is (CEEstarttime, CEEendtime)) in seconds since start
  
  #resample control data to find threshold mdist for change point
  pre <- mdo$dist[mdo$t < CEEtimes[1]]
  Nr <- 10000 #number of times to resample
  edur <- length(which(mdo$t >= CEEtimes[1] & mdo$t <= CEEtimes[2])) #duration of exposure in distance samples
  S <- sample(1:(length(pre)-edur),Nr/2,replace=TRUE) #start indices for resampled bits
  THpre <- matrix(0,nrow=Nr/2,ncol=1)
  for (k in 1:(Nr/2)){
    THpre[k] <- max(pre[S[k]:(S[k]+edur)],na.rm=TRUE)
   }           
  TH <- quantile(THpre,quantile.thr) 
  ix <- which(mdo$dist > TH)
  kk <- which(diff(mdo$t[ix])>600)#without a more than 10 minute gap of lower dist
  ix2 <- ix; ixs <- ix; resp.cst <- ix
  if (length(kk)>0){
    ix2 <- ix[1:(kk[1]-1)]
    ixs <- seq(from=min(ix2), to=max(ix2))
  }
  #get the cst of the "response"
  resp.cst <- mdo$t[c(min(ixs),max(ixs))]
  
  TH2 <- 1.6*sd(THpre) + median(THpre)
  ixb <- which(mdo$dist > TH2)
  kkb <- which(diff(mdo$t[ixb])>600)#without a more than 10 minute gap of lower dist
  ix2b <- ixb; ixsb <- ixb; resp.cst2 <- ixb
  if (length(kkb)>0){
    ix2b <- ixb[1:(kkb[1]-1)]
    ixsb <- seq(from=min(ix2b), to=max(ix2b))
  }
  #get the cst of the "response"
  resp.cst2 <- mdo$t[c(min(ixsb),max(ixsb))]
  
  return(list(resp.cst=resp.cst, resp.cst2=resp.cst2,
              resp.dur=diff(resp.cst), resp.dur2=diff(resp.cst2),
              resp.ind=ix, resp.ind.nogap=ixs,
              resp.ind2=ixb, resp.ind.nogap2=ixsb,
              TH.95prctile=TH, TH.1p6sd=TH2,
              CEEtimes=CEEtimes, mdist=mdo,
              md.params=md.params))
}