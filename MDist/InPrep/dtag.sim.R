dtag.sim <- function(n.dives, dives, G, G.resp, mvn.est, mvn.est.resp, inPhase, 
                     inPhase.resp, resp.st, resp.dur, tp){
  #simulate dtag data
  #dives is a table of real data
  #G is TPM for dives in dives
  #mvn.est is dive phase parameters
  #response is the type of CEE response observed, 
  #resp.st is start time in DIVES
  #and resp.dur is response duration in DIVES
  #tp are parameters for box-cox transforms
  
  ################################################
  rm(sim.dives,sim.data, sim.starts, sim.ends)
  sim.dives <- data.frame(dnum=c(1:n.dives))
  #randomly select the first dive in the record
  d1 <- sample(dives$clus.num , 1)
  sim.dives$resp.ind <- factor(rep("baseline", times=n.dives), 
                               levels=c("baseline", "response"))
  if (resp.dur > 0){
    sim.dives$resp.ind[resp.st:(resp.st+resp.dur-1)] <- "response"
  }
  
  MC.sim <- function(n,G,G.resp,d1) {
    sim <- numeric(n)
    m <- ncol(G)
    if (missing(d1)) {
      sim[1] <- sample(1:m,1) # random start
    } else { sim[1] <- d1 }
    for (i in 2:n) {
      if (sim.dives$resp.ind[i]=="baseline"){
        newstate <- sample(1:m,1,prob=G[sim[i-1],])
      }else{
        #if this is a "response" dive
        newstate <- sample(1:m,1,prob=G.resp[sim[i-1],])
      }
      sim[i] <- newstate
    }
    return(sim)
  }
  
  sim.dives$clus.num <- MC.sim(n.dives, G, G.resp, d1)
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