library(coda)
library(ggmcmc)
library(tidyverse)

mclist <- ggs(as.mcmc.list(map(chain_output, ~.x$samples$samples)))
mclist <- ggs(as.mcmc(out_v1_State1$samples))
ggs_traceplot(mclist,family = 's.B.int')
ggs_traceplot(mclist,family = 'r.B.int')
ggs_traceplot(mclist,family = 'f.B.int')
ggs_traceplot(mclist,family = 'f.B.int\\[1\\]')
ggs_traceplot(mclist,family = 'sd.yr')
ggs_traceplot(mclist,family = 'Sigma.yr')
ggs_traceplot(mclist,family = 'Sigma.farm')
chain_output[[1]]$dur
