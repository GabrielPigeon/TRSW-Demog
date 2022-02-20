#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# args:  [name] [array] [startState.rds]
# startState not yet implemented
if(length(args)>0) curModelName=args[1]
if(length(args)==0) curModelName='v1'


library(parallel)
library(doParallel)
library(nimble)
library(coda)
library(tidyverse)
# library(nimbleEcology)
source('R/999_MyFunc.R')

# source('R/1_cleanMultiState.R')
load('cache/cleanMultiState.Rdata')
source(paste0('R/3_CMR_',curModelName,'.R'))


# checks   ------------------
#mydat$nbFledge[mydat$state==2] %>% table(useNA = 'a')
#mydat$nbFledge[mydat$state==3] %>% table(useNA = 'a')
#mydat$nbFledge[mydat$state==1] %>% table(useNA = 'a')
#mydat$nbFledge[mydat$obs==2] %>% table(useNA = 'a')
#mydat$state[mydat$obs==2] %>% table(useNA = 'a')
#mydat$nbFledge[mydat$state==1] <- 0
# mydat$obs %>% table(useNA = 'a')
# tmp <- sapply(1:nrow(mydat$obs), function(i) mydat$obs[i,myconst$first[i]]  )
# tmp%>% table()
# mydat$state[which(tmp==3),]
# myconst$first[which(tmp==3)]
#
# sum(mydat$state %in% 1:2)-nrow(mydat$state)  # effective sample size for prob. to transit to breeding/non-breeding
#
# mydat$state %>% table(useNA = 'a')
# myconst$first %>% table(useNA = 'a')
#
# sapply(1:nrow(mydat$obs), function(i) mydat$state[i,myconst$first[i]]  ) %>% table(useNA = 'a')
# sapply(1:nrow(mydat$obs), function(i) mydat$nbFledge[i,myconst$first[i]]  ) %>% table()
# sapply(1:nrow(mydat$obs), function(i) myconst$age[i,myconst$first[i]]  ) %>% table(useNA = 'a' )
# sapply(1:nrow(mydat$obs), function(i) mydat$farm[i,myconst$first[i]]  ) %>% table(useNA = 'a' )
# sapply(1:3 , function(x ) mydat$state[mydat$obs==x] %>% table(useNA = 'a') )
# sapply(1:3 , function(x ) mydat$obs[mydat$state==x] %>% table(useNA = 'a') )



paraNimble <- function(seed,curCode,curConst,curDat,curInits=myInits,nburn=1,ntin=1,nkeep=1000, modName='curMod' ,checkpt=NULL){
    #     curConst=microConst ;curDat=microDat   ;  seed=1 ;modName='m0ffarm'  ; checkpt=4; nkeep=500; ntin=1; nburn=200
    X=seed
    strt=Sys.time()
    library(nimble)
    library(coda)
    source('R/999_MyFunc.R')
    # library(MCMCpack)
    # nb.t=curConst$nb.t
    # nb.id=curConst$nb.id
    # state=curDat$state
    # nb.mat=curConst$nb.mat
    # nbFledge=curDat$nbFledge
    # RS=curDat$RS
    # obs=curDat$obs
    # first=curConst$first
    # farm=curDat$farm
    # nE.s=dim(curDat$x.farmYrEnv)[3]
    # nE.r=dim(curDat$x.farmYrEnv)[3]
    # nE.f=dim(curDat$x.farmYrEnv)[3]
    
    
    myMod <- nimbleModel(code = curCode,
                         constants = curConst,data = curDat,
                         inits = curInits(curDat,curConst)
    )
    # myMod$getDependencies('sig')
    # myMod$getDependencies('f.B.int[2]')
    # myMod$getNodeNames() %>% length
    
    
    vars=c('s.B.int','r.B.int','f.B.int',
           'Sigma.id','rho.id', 'eps.id', 'xi.id',
           'Sigma.yr','rho.yr',  'xi.yr','eps.yr' ,
           'Sigma.farm', 'rho.farm','xi.farm', 'eps.farm'  ,
           's.B.Env','r.B.Env','f.B.Env',
           's.ranef.yr','r.ranef.yr','f.ranef.yr',
           'sig',
           # 'Pim','Nm',
           # 'Notzero',
           # 'obs.hat','nbFledge.hat',
           # 'state', 'nff',
           # 'obs',
           # 'farm',
           'mu.p','p1','p2','sd.p')
    
    Confmcmc <- configureMCMC(myMod,monitors=vars, enableWAIC = T)
    # Confmcmc$removeSamplers(c('sig','omsig'))
    # Confmcmc$addSampler(target = c('sig','omsig'),
    #            type = "RW_block",
    #            control = list(adaptInterval = 20, tries = 2))
    mymcmc <- buildMCMC(Confmcmc)
    cModel <- compileNimble(myMod)
    CmyMCMC <- compileNimble(mymcmc,project = cModel)
    dur.Compile=Sys.time()-strt
    
    CmyMCMC$run(nburn)
    
    onenit <- runMCMC(CmyMCMC,samplesAsCodaMCMC = F,
                      nburnin = 0,niter = 1,thin = 1, WAIC = F)
    output <- list(samples=matrix(NA,nrow=nkeep,
                                  ncol = ncol(onenit),
                                  dimnames = list(NULL,colnames(onenit))),
                   stateList=list(),
                   dur=list(),
                   WAIC=list())
    rm(onenit)
    output$stateList <- list(modelState = getModelState(cModel),
                             mcmcState = getMCMCstate(Confmcmc, CmyMCMC))
    output$dur=dur.Compile
    write_rds(output, file = paste0('cache/out_',modName,'_State',X,'.rds'),compress = 'gz')
    
    if(!is.null(checkpt) & checkpt>1){
        nkeepPercp <- nkeep/checkpt
        for(ii in 1:checkpt){
            output$samples[(ii-1)*nkeepPercp+ 1:nkeepPercp,] <- runMCMC(CmyMCMC,samplesAsCodaMCMC = F,
                                                                        nburnin = 0,niter = nkeepPercp*ntin,thin = ntin, WAIC = F)
            # CmyMCMC$run(pnit,thin =ntin, reset=F)
            output$stateList=list(modelState = getModelState(cModel),
                                  mcmcState = getMCMCstate(Confmcmc, CmyMCMC),
                                  code=curCode)
            output$dur=list(dur.Compile=dur.Compile,dur.tot=Sys.time()-strt)
            write_rds(output, file = paste0('cache/out_',modName,'_State',X,'.rds'),compress = 'gz')
            print(paste('saved checkpoint nb.',ii))
        }
    }else{
        # CmyMCMC$run(nkeep*ntin,thin = ntin, reset=T)
        output$samples <- runMCMC(CmyMCMC,samplesAsCodaMCMC = T,summary = T,
                                  nburnin = 0,niter = nkeep*ntin,thin = ntin, WAIC = F)
    }
    # library(nimble)
    output$stateList=list(modelState = getModelState(cModel),
                          mcmcState = getMCMCstate(Confmcmc, CmyMCMC),
                          code=curCode)
    output$dur=list(dur.Compile=dur.Compile,dur.tot=Sys.time()-strt)
    output$WAIC <- CmyMCMC$getWAIC()
    write_rds(output, file = paste0('cache/out_',modName,'_State',X,'.rds'),compress = 'gz')
    
    return(output)
}



# run the models  -----------------------
if(length(args)>1) {
    chain_output=paraNimble(seed = as.numeric(args[2]),
                            # nburn=50000,ntin=5,nkeep=2000, # waic Comparaison
                            # nburn=80000,ntin=10,nkeep=2000, # longuer run
                            # curConst=myconst,curDat=mydat,
                            # nburn=2,ntin=1,
                            curConst=microConst,curDat=microDat,
                            # curConst=miniConst,curDat=miniDat,
                            curCode=myCode,curInits = myInits,
                            modName=curModelName,checkpt=5
    )
}else{
    start.t <- Sys.time()
    this_cluster <- makeCluster(2)
    chain_output <- parLapply(cl = this_cluster, X = 1:2,
                              fun = paraNimble,
                              nburn=50000,ntin=20,
                              # nburn=50000,ntin=5,nkeep=2000, # waic Comparaison
                              # nburn=80000,ntin=10,nkeep=2000, # longuer run
                              # curConst=curConst,curDat=curDat,
                              # nburn=10000,ntin=10,
                              curConst=microConst,curDat=microDat,
                              # curConst=miniConst,curDat=miniDat,
                              curCode=myCode,curInits = myInits ,
                              modName=curModelName,checkpt=0
    )
    # It's good practice to close the cluster when you're done with it.
    stopCluster(this_cluster)
    dur=Sys.time()-start.t
    save(chain_output,myCode,file = paste0('cache/out_',curModelName,'.Rdata'))
}



# chain_output <- list(
#   read_rds(paste0('cache/out_',curModelName,'_State',1,'.rds')),
#   read_rds(paste0('cache/out_',curModelName,'_State',2,'.rds')),
#   read_rds(paste0('cache/out_',curModelName,'_State',3,'.rds'))
# )
# save(chain_output,myCode,file = paste0('cache/out_',curModelName,'.Rdata'))


