#' Simulate data from an animal-borne tag
#'
#' This function generates simulated data from an animal-borne tag. It was written for DTAG data collected on whales, but may be applicable in other contexts; more changes would be required if the animal doesn't do depth/altitude excursions (dives or flights) most of the time, or if the data are not high-resolution (DTAGs collect acceleration and magnetic field and depth data 5-500 times per second).
#' @param n number of dives (or flights, although this is untested) to simulate
#' @param G transition probability matrix indicating probability of transition between a set of dive (or perhaps flight) types
#' @param G_exp transition probability matrix during experimental exposure (if any). Defaults to G (i.e., no exposure or no change in TPM during exposure).
#' @param first integer; dive type for the first dive.  1 corresponds to the first row of G, 2 to the second, etc. If NULL (the default), the first dive type will be a random sample from the stationary distribution determined from G. 
#' @param phase_params parameters of multivariate normal distribution fit to (perhaps transformed) durations of different dive phases, by dive type. Assumed to be output from mclust::mvn, in terms of format.
#' @param phase_params_exp same as mvn_params, but during experimental exposure (if any). Defaults to mvn_est (i.e., no exposure or no change in MVN params during exposure)
#' @param fine_params parameters of multivariate AR(1) process fit to data for each observed tag data stream (e.g. acceleration, depth, pitch, roll...) within each phase. Assumed to be output from mAr::mAr.est, in terms of format.
#' @param fine_params_exp same as fine_params, but during experimental exposure (if any). Defaults to fine_params (i.e., no exposure or no change in within-phase fine-scale data during exposure).
#' @param transform_params parameters related to normalizing Box-Cox transformations of dive phase data (if used before estimating MVN parameters for dive phase durations). 
#' @param resp_start start time (as dive number) of simulated "response" to experimental exposure, if any. Defaults to 0 (no response).
#' @param resp_dur duration (in dives) of simulated "response" to experimental exposure. Defaults to 0. (NOTE - IT MAY BE OF INTEREST TO ALLOW FOR THIS and resp_start TO BE SPECIFIED IN UNITS OF EITHER DIVES, OR TIME)
#' @return a data.table object containing simulated tag data. Columns in the data.table, and their names, will correspond to inputs provided to tag_sim.
#' @importFrom magrittr "%>%"
#' @export

tag_sim <- function(n, G, G_exp=G, first=NULL, phase_params, phase_params_exp,
                    fine_params, fine_params_exp=fine_params, 
                    transform_params=NULL,
                    resp_start=0, resp_dur=0){
  
  # INPUT CHECKING AND DEFAULTS
  # -----------------------------------------------------------------------
  # if first dive type is not provided, find it as a random draw from the stationary distribution
  if (missing(first) || is.null(first)){
    # find first left eigenvector of dive transition matrix
    stationary <- eigen(t(G))$vectors[,1]
    # normalize eigenvector so entries sum to 1
    stationary <- first/sum(first)
    # draw a random integer from 1:n (where n is number of dive types) with probability [stationary] of each dive type
    first <- which.max(rmultinom(n=1, size=1, prob=stationary))
  }
  
  # SIMULATE DIVE TYPE TIME-SERIES
  # -----------------------------------------------------------------------
  
  # set up data frame for dive-level data: dive number, dive type, whether it's baseline or response, and duration
  sim_dives <- data.frame(dnum=seq(from=1, to=n)) %>%
    # add column indicating whether each dive will be "baseline" or "response" to experimental exposure.
    dplyr::mutate(resp = ifelse(dnum >= resp_start & dnum <= resp_start + resp_dur,
                         'response', 'baseline'))

  #simulate dive types based on G and G_exp
  sim.dives$clus.num <- MC.sim(n, G, G_exp, first)
  
  # SIMULATE DURATIONS OF PHASES WITHIN EACH DIVE
  # WORKING HERE 10/17
  # -----------------------------------------------------------------------
  
  #next simulate duration of each phase of each dive.
  colx <- c((ncol(sim.dives)+1):(ncol(sim.dives)+4))
  sim.dives[,colx] <- NA
  names(sim.dives)[colx] <- 
    c("b.dur", "d.dur", "a.dur", "s.dur")
  lmat.sim <- tp[sim.dives$clus.num]
  lmat.sim <- matrix(unlist(lmat.sim), ncol=4, byrow=T)
  
  k.best <- nlevels(as.factor(dives$clus.num))
  my.boxcox <- function(x, lmat) {
    #x should be a matrix of dive phase durations
    #lmat should be a matrix r of transform params
    (x^lmat-1)/lmat
  }
  my.boxcox.inv <- function(x, lmat){
    (x*lmat + 1)^(1/lmat)
  }
  for (dt in 1:k.best){
    while (sum(is.na(sim.dives[sim.dives$clus.num==dt , 
                               c("b.dur", "d.dur",
                                 "a.dur", "s.dur")])) > 0){ 
      btt <- which(sim.dives$clus.num==dt & sim.dives$resp.ind=="baseline")
      if (length(btt) > 0){
        sim.dives[btt , 
                  c("b.dur", "d.dur",
                    "a.dur", "s.dur")] <-
          my.boxcox.inv(rmvnorm(n=length(btt),
                                mean=mvn.est[[dt]]$parameters$mean , 
                                sigma=mvn.est[[dt]]$parameters$variance$Sigma,
                                method="svd") , lmat.sim[btt,])
      }
      rtt <- which(sim.dives$clus.num==dt & sim.dives$resp.ind=="response")
      if (length(rtt) > 0){
        sim.dives[rtt, 
                  c("b.dur", "d.dur",
                    "a.dur", "s.dur")] <-
          my.boxcox.inv(rmvnorm(n=length(rtt),
                                mean=mvn.est.resp[[dt]]$parameters$mean , 
                                sigma=mvn.est.resp[[dt]]$parameters$variance$Sigma,
                                method="svd") , lmat.sim[rtt,])
      }
    }
    #redo if there are any NA durations (AKA "your durations aren't really 
    #MVN...happens esp with short dive phases)
  }
  
  sim.dives$dur <- apply(sim.dives[,3:6], 1, sum)
  #set up the time-series data frame
  sim.data <- data.table(t=seq(from=(1/fs),by=(1/fs),
                               to=60*sum(sim.dives$dur)))
  setkey(sim.data,t)
  sim.dives$d.end <- 
    round((sim.dives$d.dur +
             cumsum(c(0,head(sim.dives$dur,-1)))) * fs*60)
  sim.dives$b.end <- 
    sim.dives$d.end + round(sim.dives$b.dur * fs*60)
  sim.dives$a.end <- 
    sim.dives$b.end + round(sim.dives$a.dur * fs*60)
  sim.dives$s.end <- 
    sim.dives$a.end + round(sim.dives$s.dur * fs*60)
  #add phase labels to sim.data
  #BEGIN section where I am super ashamed of the code.
  get.ind <- function(st, et){
    ix <- c(st:et)
  }
  sim.data$phase <- 
    factor("surface",
           levels=c("bottom", "descent",
                    "ascent", "surface"))
  di <- unlist(mapply(FUN=get.ind, 
                      st=c(1,head(sim.dives$s.end,-1)) ,
                      et=sim.dives$d.end))
  bi <- unlist(mapply(FUN=get.ind, 
                      st=sim.dives$d.end+1 ,
                      et=sim.dives$b.end))
  ai <- unlist(mapply(FUN=get.ind, 
                      st=sim.dives$b.end+1 ,
                      et=sim.dives$a.end))
  si <- unlist(mapply(FUN=get.ind, 
                      st=sim.dives$a.end+1 ,
                      et=sim.dives$s.end))
  di[di > nrow(sim.data)] <- nrow(sim.data)
  bi[bi > nrow(sim.data)] <- nrow(sim.data)
  ai[ai > nrow(sim.data)] <- nrow(sim.data)
  si[si > nrow(sim.data)] <- nrow(sim.data)
  sim.data$phase[di] <- "descent"
  sim.data$phase[bi] <- "bottom"
  sim.data$phase[ai] <- "ascent"
  sim.data$phase[si] <- "surface"
  #END shame. at least the extreme shame.
  
  s <- sim.data[phase=="surface",,which=T]
  sim.starts <- c(1,s[diff(sim.data[s,t])> 1]+1)
  sim.ends <- c(tail(sim.starts,-1), nrow(sim.data)+1) - 1
  #add identifier of dive-type (cluster number)
  sim.data[,clus.num:=
             unlist(mapply(FUN=function(dive.start, dive.end, clus.id)
               rep(clus.id, times=(dive.end-dive.start+1)), sim.starts, sim.ends, 
               sim.dives$clus.num))]
  #add dive number (in sequence - first, second, 3rd dive etc) 
  sim.data[,dive.num:=
             unlist(mapply(FUN=function(dive.start, dive.end, dive.id)
               rep(dive.id, times=(dive.end-dive.start+1)), sim.starts, sim.ends, 
               c(1:nrow(sim.dives))))]
  #fill in the simulated DTAG data within dive phases.
  #find the places where dive phase changes
  #this gives indices of the last entry within phases:
  setkey(sim.data,t)
  pe <- c(which(diff(unclass(sim.data$phase)) != 0), nrow(sim.data))
  #start indices of phases
  ps <- c(1, head(pe, -1) + 1) 
  #indicator of whether each time point is "response" or not
  sim.data[,resp.ind:=sim.dives[sim.data$dive.num,"resp.ind"]]
  
  #not needed for data.table:
  #sim.data[, keep2] <- vector("numeric", nrow(sim.data))
  #keep2 (phases to simulate):
  # list(depth, Ax, Ay, Az,cpitch, spitch, croll,sroll, shead, chead,nmsa, nodba, verticalvelocity, tip)
  for (p in 1:length(ps)) {
    pix <- seq(from = ps[p], to = pe[p])
    np <- length(pix)
    dt <- sim.data[pix[1],clus.num]
    # ph <- unclass(sim.data$phase[pix[1]])
    if (sim.data[pix[1], resp.ind]=="baseline"){
      mod <- inPhase[[dt]][[sim.data[pix[1],phase]]]
    }else{
      mod <- inPhase.resp[[dt]][[sim.data[pix[1],phase]]]
    }
    fillin <- try(mAr.sim(w = mod$wHat, A = mod$AHat,
                               C = mod$CHat, N = np),
                       silent=TRUE)
   if(class(fillin)=="try-error"){
      fillin <- NA
        #try(mAr.sim(w = mod$wHat, A = mod$AHat,
                #            C = mod$CHat, N = np),
                #    silent=TRUE)
    }
    sim.data[pix, (keep2):=fillin]
   rm(fillin)
  }
  
  #verify odba, msa are +ve
  sim.data[,nodba:=abs(nodba)]
  sim.data[,nmsa:=abs(nmsa)]
  
  # convert sin/cos heading to heading angle:
  sincos2angle <- function(sin.data, cos.data) {
    i <- complex(real = 0, imaginary = 1)
    ang <- Arg(cos.data + i * sin.data)
    return(ang)
    # three cheers for Euler!
  }
  ## if heading is modelled as sin(heading) & cos(heading)
  #sim.data[,heading:=sincos2angle(sim.data$shead,sim.data$chead)] 
  ##if heading is modelled as diff(heading)
  head0 <- runif(1, min = -pi, max = pi) + cumsum(sim.data$dhead)
  sim.data[,shead:=sin(head0)]
  sim.data[,chead:=cos(head0)]
  sim.data[,heading:=sincos2angle(sin.data = shead, cos.data = chead)]
  
  sim.data[,pitch:=sincos2angle(sim.data$spitch, sim.data$cpitch)]
  sim.data[,roll:=sincos2angle(sim.data$sroll, sim.data$croll)]
  # also generate a depth profile based on: 1. z=0 at start of
  # each descent phase 2. then use vertical velocity Note: this
  # will never work quite right because the whale won't end up
  # back at the surface at the end of the dive.
  setkey(sim.data,t)
  ds <- intersect(ps, sim.data[phase == "descent",,which=T])
  sim.data[,vd:=verticalvelocity * (1/fs)]
  # get 'dead reckoned' depth track
  vd2z <- function(dive.start, dive.end, vd) {
    cumsum(vd[c(dive.start:dive.end)]) - vd[dive.start]
  }
  sim.data[,z0:=unlist(mapply(FUN = vd2z, dive.start = sim.starts,
                              dive.end = sim.ends, MoreArgs = list(vd = sim.data$vd)))]
  # adjust ascent vv so the whale ends at the surface.
  as <- intersect(ps, sim.data[phase == "ascent",,which=T])
  ss <- intersect(ps, sim.data[phase == "surface",,which=T])
  adj.asc.vd <- function(as, ae, z0, vd) {
    # adjust ascent depth profile and vv so whale ends up at the
    # surface @ end of dive
    zb <- z0[as - 1] #depth at end of bottom phase
    za <- -sum(vd[c(as:ae)]) #net distance covered during the ascent
    vd.new <- zb/za * vd[c(as:ae)]
    return(vd.new)
  }
  sim.data[,vd.adj:=vd]
  sim.data[phase == "ascent" ,vd.adj:=unlist(mapply(FUN = adj.asc.vd,
                                                    as = as, ae = ss - 1, MoreArgs = list(vd = sim.data$vd, z0 = sim.data$z0)))]
  sim.data[,z:=unlist(mapply(FUN = vd2z, dive.start = sim.starts,
                             dive.end = sim.ends, MoreArgs = list(vd = sim.data$vd.adj)))]
  sim.data[,z:=ifelse(z < -5, 0, z)]
  return(sim.data)
  
}
