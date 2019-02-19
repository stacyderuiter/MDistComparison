# example whale dive simulation ######################
# author: katherine whyte ############################

# load Stacy example start parameters
load("~/GitHub/mdist_comparison/data/start_data_4sims.RData")

# source all scripts
source("code/test_code/box_cox.R")
source("code/test_code/mc_sim.R")
source("code/test_code/sincos2angle.R")
source("code/test_code/tag_sim.R")

# load libraries
require(mvtnorm)
require(data.table)
require(mAr)
require(CircStats)
require(magrittr)
require(dplyr)



#number of whales to simulate
nw <- 1
#number of dives to simulate, per whale
nd <- 50

# OPTION A ##
# no change in dive type sequence during response
G.resp <- G
# no changes to mvn.est (phase durations) during response
mvn.est.resp <- mvn.est
# no change in within phase measurements during response
inPhase.resp <- inPhase

# generate data
test <- tag_sim(n=nd, G=G, G_exp = G.resp, 
                first = 1, # try changing to NULL if this works
                phase_params = mvn.est, phase_params_exp=mvn.est.resp, 
                fine_params = inPhase, fine_params_exp=inPhase.resp, 
                transform_params = tp,
                sampling_rate = 1, 
                resp_start = 30, resp_dur = 10)

# plot simulated dive profiles
plot(test$t, -test$z, xlab="Time (seconds)", ylab="Depth (m)")
