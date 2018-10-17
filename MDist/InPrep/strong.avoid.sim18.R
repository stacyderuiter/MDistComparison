source("mdist.R")
source("get.chpoint.R")
source('dtag.sim.R')
require(mvtnorm)
require(data.table)
require(mAr)
require(CircStats)
#number of whales to simulate
nw <- 500#150
#number of dives to simulate, per whale
nd <- 20
md.params <- expand.grid(smoothDur=2^seq(from=1, to=5, by=2) , 
                         overlap.proportion=c(0.1, 0.9),
                         crit=0.95, consec=c(TRUE,FALSE), 
                         cumSum=FALSE)

load("C:/Users/sld33/Dropbox/MDistComparison/start_data_4sims.RData")
#for Ziphius here
#dive type 1 is shallow, 2 is deep foraging
#no change in dive types during response
G.resp <- G
# G.resp[1,1] <- G[1,1] + 0.8*(1-G[1,1])
# G.resp[1,2] <- G[1,2] - 0.8*(1-G[1,1])
# G.resp[2,1] <- G[2,1] + 0.8*(1-G[2,1])
# G.resp[2,2] <- G[2,2] - 0.8*(1-G[2,1])

# no changes to mvn.est (phase durations) during response
mvn.est.resp <- mvn.est


#make longer bottom and ascent phases for deep dives

######################################
# changes to inPhase during response:

# decrease heading and pitch variability and increase ODBA mean
######################################
inPhase.resp <- inPhase

#make heading, pitch less variable
#keep2 is a vector of variable names in inPhase
pitch.vv <- is.element(keep2, c("cpitch", "spitch", "dhead"))
var.decr.factor <- 2 #new var is old/var.decr.factor
for (dtype in c(1:2)){
  for (phase in c(1:length(inPhase.resp[[2]]))){
  diag(inPhase.resp[[dtype]][[phase]]$Chat)[pitch.vv] <- diag(inPhase.resp[[dtype]][[phase]]$Chat)[pitch.vv]*var.decr.factor
  }}

#make ODBA higher overall
mean.incr.factor <- 2
od <- is.element(keep2, "nodba")
for (dtype in c(1:2)){
  #loop over all dive phases, desc bott asc surf
  for (ph in c(1:4)){
    inPhase.resp[[dtype]][[ph]]$wHat[od] <- inPhase[[dtype]][[ph]]$CHat[od]*mean.incr.factor
  }
}


tryCatch.W.E <- function(expr)
  {
      W <- NULL
      w.handler <- function(w){ # warning handler
  	W <<- w
       	invokeRestart("muffleWarning")
           }
         list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                            				     warning = w.handler),
               	 warning = W)
     }

run.sim <- function(nw,nd,md.params){
  avoid.results <- list()
  feed.results <- list()
  lit.results <- list()
  parset.counter <- 0
  

    for (whale in c(1:nw)){
      
      rm(wdat)
      #if (whale/20 == round(whale/20)){
      print(paste("Whale ", whale,
                  " of ", nw, "...", date(),sep="" ))
      #}
      wdat <- dtag.sim(n.dives=nd, dives, G, G.resp, mvn.est, 
                       mvn.est.resp, inPhase, 
                       inPhase.resp=inPhase.resp, resp.st=6, resp.dur=2, tp)
      #variability of heading
      fs <-  1/mean(diff(wdat$t))
      #averaging period for calculations (var.head, etc) in seconds:
      varsec <- 60
      #var.head <- numeric(nrow(wdat))
      #     circ.var <- function(X) { circ.disp(X)$var }
      #     #ROLLAPPLY IS SO SLOW TOO!
      #     var.head <- rollapply(data=wdat$heading, width=varsec*fs,
      #                           FUN=circ.var, fill=NULL, partial=TRUE, 
      #                           align="center")
      get.var <- function(ix, data, win.size, circ){
        if (circ==TRUE){
          v <- circ.disp(data[seq(from=max(c(1, ix-win.size/2)),
                                  to=min(c(length(data),ix+win.size/2)),
                                  by=1)])$var
        }else{
          v <- var(data[seq(from=max(c(1, ix-win.size/2)),
                            to=min(c(length(data),ix+win.size/2)),
                            by=1)])
        }
        return(v)
      }
      wdat[,var.head := vapply(X=c(1:length(wdat$heading)) , FUN=get.var,
                               FUN.VALUE=1, data=wdat$heading, 
                               win.size=ceiling(varsec*fs),
                               circ=TRUE)]
      wdat[,var.pitch := vapply(X=c(1:length(wdat$pitch)) , FUN=get.var,
                                FUN.VALUE=1, data=wdat$pitch, 
                                win.size=ceiling(varsec*fs),
                                circ=TRUE)]
      wdat[,var.odba := vapply(X=c(1:length(wdat$nodba)) , FUN=get.var,
                               FUN.VALUE=1, data=wdat$nodba, 
                               win.size=ceiling(varsec*fs),
                               circ=FALSE)]
      wdat[,var.vv := vapply(X=c(1:length(wdat$verticalvelocity)) , FUN=get.var,
                               FUN.VALUE=1, data=wdat$verticalvelocity, 
                               win.size=ceiling(varsec*fs),
                               circ=FALSE)]
      #this loop is WAY TOO SLOW -- maybe in RCpp would be ok...
      #may want to try it if time because the vapply is still slow.
      #     var.pitch <- numeric(nrow(wdat))
      #     var.odba <- numeric(nrow(wdat))#
      #     for (j0 in c(1:nrow(wdat))){
      #       if (j0==1){
      #         var.head[j0] <- 0
      #         var.pitch[j0] <- 0
      #         var.odba[j0] <- 0
      #       }else{ 
      #         var.head[j0] = circ.disp(wdat$heading[seq(from=max(c(1, j0-(varsec/2*fs))),
      #                                               to=min(c(nrow(wdat),(j0+varsec/2*fs))),
      #                                               by=1)])$var
      #         var.pitch[j0] = circ.disp(wdat$pitch[seq(from=max(c(1, j0-(varsec/2*fs))),
      #                                               to=min(c(nrow(wdat),(j0+varsec/2*fs))),
      #                                               by=1)])$var
      #         var.odba[j0] = var(wdat$nodba[seq(from=max(c(1, j0-(varsec/2*fs))),
      #                                               to=min(c(nrow(wdat),(j0+varsec/2*fs))),
      #                                               by=1)])
      #       }
      #     }
      #     wdat[,var.head:=var.head]
      #     wdat[,var.pitch:=var.pitch]
      #     wdat[,var.odba:=var.odba]
      for (parset in c(1:nrow(md.params))){
        list.id <- (whale-1)*nrow(md.params) + parset 
        print(paste("Mdist params ", parset,
                    " of ", nrow(md.params), "...", date(), sep="" ))
        N <- md.params[parset,"smoothDur"]
      O <- md.params[parset, "smoothDur"] * md.params[parset,"overlap.proportion"]
      lit.vars <- c("nodba", "var.head", "var.pitch", "shead", "chead")
      avoid.vars <- c("var.head", "var.pitch", "nodba") 
      feed.vars <- c("var.odba", "var.pitch", "verticalvelocity", "var.vv")
      resp.start.sec <- match("response", wdat$resp.ind)/fs
      CEEstart <- min(resp.start.sec, 60*60*1.5)
      CEEtimes <- c(CEEstart, CEEstart+0.5*60*60 )
      md.avoid <- tryCatch.W.E(mdist(data=wdat[,avoid.vars,with=FALSE],
                        fs=fs, smoothDur=N, 
                        overlap=O, consec=md.params[parset,"consec"], 
                        cumSum=md.params[parset,"cumSum"], 
                        expStart=CEEstart, expEnd= CEEstart+0.5*60*60, 
                        baselineStart=0, baselineEnd=CEEtimes[1]-1))
      if(is.data.frame(md.avoid$value)){
        avoid.results[[list.id]] <- get.chpoint(mdo=md.avoid$value,
                                                quantile.thr=md.params[parset,"crit"],
                                                CEEtimes=CEEtimes,
                                                md.params=md.params[parset,])
      }else{whale <- whale - 1} #redo if error
      md.feed <- tryCatch.W.E(mdist(data=wdat[,feed.vars,with=FALSE],
                       fs=fs, smoothDur=N, 
                       overlap=O, consec=md.params[parset,"consec"], 
                       cumSum=md.params[parset,"cumSum"], 
                       expStart=CEEstart, expEnd= CEEstart+0.5*60*60,  
                       baselineStart=0, baselineEnd=CEEtimes[1]-1))
      if(is.data.frame(md.feed$value)){
        feed.results[[list.id]] <- get.chpoint(mdo=md.feed$value,
                                               quantile.thr=md.params[parset,"crit"],
                                               CEEtimes=CEEtimes,
                                               md.params=md.params[parset,])
        
      }
      md.lit <- tryCatch.W.E(mdist(data=wdat[,lit.vars,with=FALSE],
                                    fs=fs, smoothDur=N, 
                                    overlap=O, consec=md.params[parset,"consec"], 
                                    cumSum=md.params[parset,"cumSum"], 
                                    expStart=CEEstart, expEnd= CEEstart+0.5*60*60,  
                                    baselineStart=0, baselineEnd=CEEtimes[1]-1))
      if(is.data.frame(md.lit$value)){
        lit.results[[list.id]] <- get.chpoint(mdo=md.lit$value,
                                               quantile.thr=md.params[parset,"crit"],
                                               CEEtimes=CEEtimes,
                                               md.params=md.params[parset,])
        
      }
      

      
      
      
      save(avoid.results, feed.results, lit.results, file='strong.avoid.recover18.RData')
    }
  }
 
  return(list(avoid.results=avoid.results, feed.results=feed.results, lit.results=lit.results))
}

strong.avoid.sim.out <- run.sim(nw,nd,md.params)
