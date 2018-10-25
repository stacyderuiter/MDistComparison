#' simulate from a Markov chain
#'
#' @param n number of time-steps to simulate
#' @param G transition probability matrix (baseline)
#' @param G_exp transition probability matrix (exposure)
#' @param first l
#' @return a vector of numeric values (1 up to number of rows/columns in G) indicating the simulated state sequence
#' @export



mc_sim <- function(n,G,G_resp,first) {
  sim <- numeric(n)
  m <- ncol(G)
sim[1] <- first
  for (i in 2:n) {
    if (sim_dives$resp[i]=="baseline"){
      newstate <- sample(1:m,1,prob=G[sim[i-1],])
    }else{
      #if this is a "response" dive
      newstate <- sample(1:m,1,prob=G_exp[sim[i-1],])
    }
    sim[i] <- newstate
  }
  return(sim)
}