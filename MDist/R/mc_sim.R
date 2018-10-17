
mc_sim <- function(n,G,G.resp,d1) {
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