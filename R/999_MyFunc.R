library(tidyverse)
library(stringr)
library(lubridate)
library(magrittr)
library(nimble)

theme_gab <- function(base_size=12,base_family = "Sans"){
    theme_bw(base_size,base_family)+
        theme(panel.border = element_rect(linetype = "solid", colour = "black",size=0.5, fill=NA),
              plot.margin= unit(c(0.25,0.25,0.25,0.25), "lines"),
              axis.text	=element_text(colour = "black")

        )
}



printMeanCI <- function(x,quant=c(0.025,0.975),nsmall=3,sep=", ",par="["){
    if(par %in% c(")",'(','()')) {par=c(' (',')') } else{par=c('[',']')}

    mean=format(round(mean(x),nsmall),nsmall = nsmall)
    q=format(round(quantile(x,probs = quant),nsmall),nsmall = nsmall)
    paste0(mean,par[1],q[1],sep,q[2],par[2])
}

ArchiveSave <- function(...,file,list=NULL){
    if(length(file)>1 | !is.character(file)) stop('incompatible file names')

	file <- str_split(file,pattern = '/')[[1]]
	if(length(file)>1) {
	    folder=paste(file[-length(file)],sep='/')
	    file=file[length(file)]
	} else {
	    folder=NULL
	    file=file
	}

	file1=paste(folder,file,sep="/")
	save(...,list = list, file = file1,compress = 'xz')

	date <- str_replace_all(today(), fixed("-"), "")
	if(!dir.exists(paste0(folder,'/Archive'))) dir.create(paste0(folder,'/Archive'))
	file2=paste0(folder,'/Archive/',date,'_',file)
	save(...,list = list, file = file2,compress = 'xz')

}

getPace <- function(coda,dur,param){
  # coda <- as.mcmc(samples$samples)[,param]
ef <-effectiveSize(coda[,param])
print('effective sample / hours')
ef/as.numeric(dur,units='hours')
}


dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(),
                 notZeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(notZeroProb))
      ## or the probability if log = FALSE
      else return((notZeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- (1-notZeroProb) + (notZeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })
rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), notZeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = (1-notZeroProb), size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })
registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, notZeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'notZeroProb = double()')
  )))


g.summaryBayes <- function(mcmc,null=0){
  library(coda)
  library(purrr)
  library(dplyr)

  if(class(mcmc)=="mcmc.list"){
    tmp <- summary(mcmc)
    tmp=data.frame(var=rownames(tmp[[1]]),
                   tmp[[1]][,c(1,2)],
                   effSize=effectiveSize(mcmc),
                   tmp[[2]][,c(1,5)],
                   row.names = NULL
    )
    colnames(tmp)[c(5,6)] <- c('cil','cih')
    pdir <- mcmc %>% map_dfr(as.data.frame) %>% apply(2,function(x){
      max(
            c(
              length(x[x > null]) / length(x), # pd positive
              length(x[x < null]) / length(x) # pd negative
            )
          )

  })

    out <- tmp %>% mutate(Pdir=pdir,
                   S=ifelse((cil* cih)>0,'*','')
                   ) %>%
      as_tibble()

  }

  if(class(mcmc) %in% c('numeric','matrix',"mcmc")){
    mcmc <- as.data.frame(mcmc)
  }

  if(class(mcmc)=='data.frame'){

    tmp=data.frame(var=colnames(mcmc),
                   Mean=colMeans(mcmc),
                   SD=apply(mcmc, 2, sd),
                   effSize=apply(mcmc, 2, function(x) coda::effectiveSize(x)),
                   cil=apply(mcmc, 2, quantile,0.025),
                  cih=apply(mcmc, 2, quantile,0.975),
                   row.names = NULL

    )
    pdir <- mcmc %>% map_dfr(as.data.frame) %>% apply(2,function(x){
      max(
        c(
          length(x[x > null]) / length(x), # pd positive
          length(x[x < null]) / length(x) # pd negative
        )
      )

    })

    out <- tmp %>% mutate(Pdir=pdir,
                   S=ifelse((cil* cih)>0,'*','')
    ) %>%
      as_tibble()

  }
return(out)
}



getStateVariableNames <- function(samplerDef) {
    resetMethod <- body(samplerDef$reset)
    stateVars <- character()
    if(resetMethod[[1]] != '{') stop('something wrong')
    numLines <- length(resetMethod)
    for(i in 1:numLines) {
        if(i == 1) next
        thisLine <- resetMethod[[i]]
        if(thisLine[[1]] == '<<-') {
            LHS <- thisLine[[2]]
            if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
            stateVars <- c(stateVars, as.character(LHS))
        }
        if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
            stateVars <- c(stateVars, 'my_calcAdaptationFactor')
        }
    }
    setupMethod <- body(samplerDef$setup)
    if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
    return(stateVars)
}

getModelState <- function(model) {
    modelVarNames <- model$getVarNames()
    modelVarValuesList <- vector('list', length(modelVarNames))
    names(modelVarValuesList) <- modelVarNames
    for(var in modelVarNames) {
        modelVarValuesList[[var]] <- model[[var]]
    }
    return(modelVarValuesList)
}

getMCMCstate <- function(conf, mcmc) {
    stateVarNamesList <- vector('list', length(conf$samplerConfs))
    mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
    for(i in seq_along(conf$samplerConfs)) {
        samplerDef <- conf$getSamplerDefinition(i)
        theseStateNames <- getStateVariableNames(samplerDef)
        theseStateValuesList <- vector('list', length(theseStateNames))
        names(theseStateValuesList) <- theseStateNames
        for(j in seq_along(theseStateNames)) {
            if(is.nf(mcmc)) {
                if(theseStateNames[j] == 'my_calcAdaptationFactor') {
                    theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                                      gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
                } else
                    theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
            }
            if(is.Cnf(mcmc)) {
                if(theseStateNames[j] == 'my_calcAdaptationFactor') {
                    theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                                      gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
                } else
                    theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
            }
        }
        mcmcStateValuesList[[i]] <- theseStateValuesList
    }
    return(mcmcStateValuesList)
}

setModelState <- function(model, modelState) {
    modelVarNames <- model$getVarNames()
    if(!identical(sort(modelVarNames), sort(names(modelState)))) stop('saved model variables don\'t agree')
    for(var in modelVarNames) {
        model[[var]] <- modelState[[var]]
    }
    invisible(model$calculate())
}

setMCMCstate <- function(conf, mcmc, mcmcState) {
    if(length(mcmcState) != length(conf$samplerConfs)) stop('saved mcmc samplers don\'t agree')
    for(i in seq_along(conf$samplerConfs)) {
        theseStateValuesList <- mcmcState[[i]]
        for(j in seq_along(theseStateValuesList)) {
            samplerStateName <- names(theseStateValuesList)[j]
            if(is.nf(mcmc)) {
                if(samplerStateName == 'my_calcAdaptationFactor') {
                    mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$timesAdapted <- theseStateValuesList[[j]]$timesAdapted
                    mcmc$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$gamma1 <- theseStateValuesList[[j]]$gamma1
                } else {
                    mcmc$samplerFunctions$contentsList[[i]][[samplerStateName]] <- theseStateValuesList[[samplerStateName]]
                }
            }
            if(is.Cnf(mcmc)) {
                if(samplerStateName == 'my_calcAdaptationFactor') {
                    invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted', theseStateValuesList[[j]]$timesAdapted))
                    invisible(valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1', theseStateValuesList[[j]]$gamma1))
                } else {
                    invisible(valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], samplerStateName, theseStateValuesList[[samplerStateName]]))
                }
            }
        }
    }
}