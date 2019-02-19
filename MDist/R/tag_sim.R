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
#' @param sampling_rate sampling rate in Hz (samples per second) for fine-scale time-series data. Defaults to 1 Hz.
#' @param transform_params parameters related to normalizing Box-Cox transformations of dive phase data (if used before estimating MVN parameters for dive phase durations). 
#' @param resp_start start time (as dive number) of simulated "response" to experimental exposure, if any. Defaults to 0 (no response).
#' @param resp_dur duration (in dives) of simulated "response" to experimental exposure. Defaults to 0. (NOTE - IT MAY BE OF INTEREST TO ALLOW FOR THIS and resp_start TO BE SPECIFIED IN UNITS OF EITHER DIVES, OR TIME)
#' @return a data.table object containing simulated tag data. Columns in the data.table, and their names, will correspond to inputs provided to tag_sim.
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":="
#' @export

tag_sim <- function(n, G, G_exp=G, first=NULL, phase_params, phase_params_exp = phase_params,
                    fine_params, fine_params_exp=fine_params, 
                    transform_params=NULL, sampling_rate=1,
                    resp_start=0, resp_dur=0){
  
  # INPUT CHECKING AND DEFAULTS
  # -----------------------------------------------------------------------
  # if first dive type is not provided, find it as a random draw from the stationary distribution
  if (missing(first) || is.null(first)){
    # find first left eigenvector of dive transition matrix
    stationary <- eigen(t(G))$vectors[,1]
    # normalize eigenvector so entries sum to 1
    stationary <- stationary/sum(stationary)
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
  sim_dives <- dplyr::mutate(sim_dives, dive_type = mc_sim(sim_dives, n, G, G_exp, first))

  # SIMULATE DURATIONS OF PHASES WITHIN EACH DIVE
 
  # next simulate duration of each phase of each dive.
  # simulation uses multivariate normal distribution fitted to dive durations by phase,
  # params in inputs phase_params and phase_params_exp
  # possibly transformed via Box-Cox power transforms; if so params are in input transform_params
  # make a matrix (4 columns, same number of rows as sim_dives)
  
  phase_durs  <-  function(dive_type, nn=1, phase_params, transform_params=NULL){
    # get the right box-cox transformation parameters lambda for this dive type
    tp <- transform_params[[dive_type]]
    # generate random dive phase durations from MVN distribution
    durs <- mvtnorm::rmvnorm(nn, mean=phase_params[[dive_type]]$parameters$mean,
                                   sigma=phase_params[[dive_type]]$parameters$variance$Sigma,
                                   method="svd") %>%
      # make the results into a data frame
      data.frame() %>%
      # change the names of the durations columns to match the dive phases (fetch names from phase_params)
      rlang::set_names(row.names(phase_params[[1]]$parameters$mean))
    # undo the box-cox transformation, if necessary
    if (!is.null(transform_params)){
      return(durs)
    }else{
      return(box_cox(durs, lambda=tp, inverse=TRUE))
    }
  }
  
  # apply the phase_durs function to get durations for the dives in sim_dives
  sim_dives <- sim_dives %>%
    mutate(dur = purrr::map(dive_type, function(x) phase_durs(x, nn=1, phase_params=phase_params, transform_params= transform_params))) %>%
    tidyr::unnest(dur) %>%
    #Note: it may happen that there are some NA phase durations (AKA "your durations aren't really 
    #MVN...happens esp with short dive phases)
    # get the duration of the full dive
    mutate(dur = a.dur + b.dur + s.dur + d.dur)
# ---------------------------------------------------- end of "dive metadata" table (sim_dives) creation ---------------------------
  
# ---------------------------------------------create sim_data (fine-scale time-series data.table) ---------------------------------
# Note : Code was modernized up to here but :( no further yet.
    #set up the time-series data frame
  sim_data <- data.table::data.table(t=seq(from=(1/sampling_rate),by=(1/sampling_rate),
                               to=60*sum(sim_dives$dur)))
  #organize sim_data by time
  data.table::setkey(sim_data,t)
  # create phase end times using phase durations 
  sim_dives$d.end <- 
    round((sim_dives$d.dur +
             cumsum(c(0,head(sim_dives$dur,-1)))) * sampling_rate*60)
  sim_dives$b.end <- 
    sim_dives$d.end + round(sim_dives$b.dur * sampling_rate*60)
  sim_dives$a.end <- 
    sim_dives$b.end + round(sim_dives$a.dur * sampling_rate*60)
  sim_dives$s.end <- 
    sim_dives$a.end + round(sim_dives$s.dur * sampling_rate*60)
  #add phase labels to sim_data
  #BEGIN section where I am super ashamed of the code. this is old and slow.
  get.ind <- function(st, et){
    ix <- c(st:et)
  }
  # allocate space for phase data
  sim_data$phase <- 
    factor("surface",
           levels=c("bottom", "descent",
                    "ascent", "surface"))
  # find indices of  descent...
  di <- unlist(mapply(FUN=get.ind, 
                      st=c(1,head(sim_dives$s.end,-1)) ,
                      et=sim_dives$d.end))
  # bottom
  bi <- unlist(mapply(FUN=get.ind, 
                      st=sim_dives$d.end+1 ,
                      et=sim_dives$b.end))
  # ascent
  ai <- unlist(mapply(FUN=get.ind, 
                      st=sim_dives$b.end+1 ,
                      et=sim_dives$a.end))
  # and surface phases
  si <- unlist(mapply(FUN=get.ind, 
                      st=sim_dives$a.end+1 ,
                      et=sim_dives$s.end))
  # match lengths of labels and dataset (sometimes rounding error hurts)
  di[di > nrow(sim_data)] <- nrow(sim_data)
  bi[bi > nrow(sim_data)] <- nrow(sim_data)
  ai[ai > nrow(sim_data)] <- nrow(sim_data)
  si[si > nrow(sim_data)] <- nrow(sim_data)
  # add phase labels to phase column of sim_data
  sim_data$phase[di] <- "descent"
  sim_data$phase[bi] <- "bottom"
  sim_data$phase[ai] <- "ascent"
  sim_data$phase[si] <- "surface"
  #END shame. at least the extreme shame.
  
  # find dive start and end times
  s <- sim_data[phase=="surface",,which=T]
  sim.starts <- c(1,s[diff(sim_data[s,t])> 1]+1)
  sim.ends <- c(tail(sim.starts,-1), nrow(sim_data)+1) - 1
  #add identifier of dive-type (cluster number)
  sim_data$dive_type=
             unlist(mapply(FUN=function(dive.start, dive.end, clus.id)
               rep(clus.id, times=(dive.end-dive.start+1)), 
               sim.starts, sim.ends, sim_dives$dive_type))
  #add dive number (in sequence - first, second, 3rd dive etc) 
  sim_data$dive_num=
             unlist(mapply(FUN=function(dive.start, dive.end, dive.id)
               rep(dive.id, times=(dive.end-dive.start+1)),
               sim.starts, sim.ends, c(1:nrow(sim_dives))))
  #fill in the simulated DTAG data within dive phases.
  #find the places where dive phase changes
  #this gives indices of the last entry within phases:
  data.table::setkey(sim_data,t)
  # indices of phase ends
  pe <- c(which(diff(unclass(sim_data$phase)) != 0), nrow(sim_data))
  #start indices of phases
  ps <- c(1, head(pe, -1) + 1) 
  #indicator of whether each time point is "response" or not
  sim_data$resp=sim_dives[sim_data$dive_num,"resp"]
  # add new columns for new Dtag within phase sims
  sim_data[,as.character(keep2)] <- as.numeric(0)
  
  # loop over all dive phases in simulated dataset
  for (p in 1:length(ps)) {
    # indices of this phase
    pix <- seq(from = ps[p], to = pe[p])
    # duration of the phase in samples
    np <- length(pix)
    # which dive type this dive is
    dt <- sim_data[pix[1],dive_type]
    # ph <- unclass(sim_data$phase[pix[1]])
    
    if (sim_data[pix[1], resp]=="baseline"){
      # simulation params for this phase, this dive type, baseline
      mod <- fine_params[[dt]][[sim_data[pix[1],phase]]]
    }else{
      # simulation params for this phase, this dive type, response
      mod <- fine_params_exp[[dt]][[sim_data[pix[1],phase]]]
    }
    
    # try to simulate data
    fillin <- try(mAr.sim(w = mod$wHat, A = mod$AHat,
                               C = mod$CHat, N = np),
                       silent=TRUE)
    # if there was an error then fill in missing data
   if(class(fillin)=="try-error"){
      fillin <- NA
        #try(mAr.sim(w = mod$wHat, A = mod$AHat,
                #            C = mod$CHat, N = np),
                #    silent=TRUE)
   }
    # put the data into the sim_data data table
    sim_data[pix, 6:18:= fillin]
    # remove the phase data so it doesn't gum up the next iteration
    rm(fillin)
  }
  
  #verify odba, msa are +ve.
  # this is no good, no good at all. I seriously wonder if it would be better to just generate these from simulated Acc data
  # need to try this when I have time...
  sim_data[,nodba:=abs(nodba)]
  sim_data[,nmsa:=abs(nmsa)]
  
  ## if heading is modelled as sin(heading) & cos(heading)
  #sim_data[,heading:=sincos2angle(sim_data$shead,sim_data$chead)] 
  ##if heading is modelled as diff(heading) (this seems to work better?)
  head0 <- runif(1, min = -pi, max = pi) + cumsum(sim_data$dhead)
  # this switching back and forth is to get the headings into the 0-2pi range
  sim_data[,shead:=sin(head0)]
  sim_data[,chead:=cos(head0)]
  sim_data[,heading:=sincos2angle(sin = shead, cos = chead)]
  
  # fill in pitch and roll from simulated sin(pitch) cos(pitch) sin(roll) cos(roll)
  sim_data[,pitch:=sincos2angle(sim_data$spitch, sim_data$cpitch)]
  sim_data[,roll:=sincos2angle(sim_data$sroll, sim_data$croll)]
  # also generate a depth profile based on: 1. z=0 at start of
  # each descent phase 2. then use vertical velocity Note: this
  # will never work quite right because the whale won't end up
  # back at the surface at the end of the dive.
  data.table::setkey(sim_data,t)
  ds <- intersect(ps, sim_data[phase == "descent",,which=T])
  sim_data[,vd:=verticalvelocity * (1/sampling_rate)]
  # get 'dead reckoned' depth track
  vd2z <- function(dive.start, dive.end, vd) {
    cumsum(vd[c(dive.start:dive.end)]) - vd[dive.start]
  }
  # put initial simulated depth into the data table
  sim_data[,z0:=unlist(mapply(FUN = vd2z, dive.start = sim.starts,
                              dive.end = sim.ends, MoreArgs = list(vd = sim_data$vd)))]
  # adjust ascent vv so the whale ends at the surface.
  as <- intersect(ps, sim_data[phase == "ascent",,which=T])
  ss <- intersect(ps, sim_data[phase == "surface",,which=T])
  adj.asc.vd <- function(as, ae, z0, vd) {
    # adjust ascent depth profile and vv so whale ends up at the
    # surface @ end of dive
    zb <- z0[as - 1] #depth at end of bottom phase
    za <- -sum(vd[c(as:ae)]) #net distance covered during the ascent
    vd.new <- zb/za * vd[c(as:ae)]
    return(vd.new)
  }
  # put adjusted vertical velocity into the data table
  sim_data[,vd.adj:=vd]
  sim_data[phase == "ascent" ,vd.adj:=unlist(mapply(FUN = adj.asc.vd,
                                                    as = as, ae = ss - 1, MoreArgs = list(vd = sim_data$vd, z0 = sim_data$z0)))]
  sim_data$z = unlist(mapply(FUN = vd2z, dive.start = sim.starts,
                             dive.end = sim.ends, MoreArgs = list(vd = sim_data$vd.adj)))
  # don't let the whales fly even after all this adjustment they try to.
  sim_data[,z:=ifelse(z < -5, 0, z)]
  return(sim_data)
  
}
