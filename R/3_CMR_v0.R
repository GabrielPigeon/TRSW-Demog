# CMR mutistale joints
# survival + RS  / nbFledge
# multi-state  /  zero-inflated poisson

# states:
# 1=alive & 0 nestling fledged
# 2=alive & >0 nestling fledged
# 3=dead
# observations
# 1=seen with no nestling fledged
# 2=seen with >0 nestling fledged
# 3=not seen 


library(nimble)


# model -------------------------------------------------------------------


myCode <- nimbleCode({
    # Priors 
    for(a in 1:3){ # loop throught age class to get age specific parameter estimate (ALL is in interaction with age)
        s.B.int[a]~dlogis(0,1)   # intercept for survival
        r.B.int[a]~dlogis(0,1)   # intercept for reproductive success (>0 fledgelings)
        f.B.int[a]~dnorm(0,0.001) # intercept for poisson part of zip on nb of fledgelings produced
        
        # slopes for the effects of env covariates
        s.B.Env[1,a] <- 0#  ~dnorm(0,0.001)            # temp
        s.B.Env[2,a] <- 0#  ~dnorm(0,0.001) # prec
        s.B.Env[3,a] <- 0#  ~dnorm(0,0.001)# coldSnap
        s.B.Env[4,a] ~dnorm(0,0.001)            # hosp
        s.B.Env[5,a] <- 0# ~dnorm(0,0.001)# IA
        s.B.Env[6,a] <- 0#  ~dnorm(0,0.001)# prec:cold
        s.B.Env[7,a] <- 0#  ~dnorm(0,0.001)# prec:IA
        s.B.Env[8,a] <- 0#  ~dnorm(0,0.001)# IA:colf
        s.B.Env[9,a] <- 0#  ~dnorm(0,0.001)# triple
        
        r.B.Env[1,a] <- 0#  ~dnorm(0,0.001)            # temp
        r.B.Env[2,a] <- 0#  ~dnorm(0,0.001) # prec
        r.B.Env[3,a] <- 0#  ~dnorm(0,0.001)# coldSnap
        r.B.Env[4,a] ~dnorm(0,0.001)            # hosp
        r.B.Env[5,a] <- 0# ~dnorm(0,0.001)# IA
        r.B.Env[6,a] <- 0#  ~dnorm(0,0.001)# prec:cold
        r.B.Env[7,a] <- 0#  ~dnorm(0,0.001)# prec:IA
        r.B.Env[8,a] <- 0#  ~dnorm(0,0.001)# IA:colf
        r.B.Env[9,a] <- 0#  ~dnorm(0,0.001)# triple
        
        f.B.Env[1,a] <- 0#  ~dnorm(0,0.001)            # temp
        f.B.Env[2,a] <- 0#  ~dnorm(0,0.001) # prec
        f.B.Env[3,a] <- 0#  ~dnorm(0,0.001)# coldSnap
        f.B.Env[4,a] ~dnorm(0,0.001)            # hosp
        f.B.Env[5,a] <- 0# ~dnorm(0,0.001)# IA
        f.B.Env[6,a] <- 0#  ~dnorm(0,0.001)# prec:cold
        f.B.Env[7,a] <- 0#  ~dnorm(0,0.001)# prec:IA
        f.B.Env[8,a] <- 0#  ~dnorm(0,0.001)# IA:cold
        f.B.Env[9,a] <- 0#  ~dnorm(0,0.001)# triple
        
    }
    
    # ######  random effects   -------------------------
    # modeled using multivariate normal distribution to get age:yr random effects. and account for annual covariance in survival between traits and ages. 
    # use of a scaling parameter to help convergeance of random effect with very different scale. see Gelman et al., 
    
    for (a in 1:9){
        xi.yr[a] ~ dunif(0,10) # Scaling
    } #i
    
    for (t in 1:(nb.t-1)) {
        eps.yr[t,1:9]  ~ dmnorm(zero[1:9], Tau.yr[1:9, 1:9])
        for (a in 1:3){
            s.ranef.yr[t,a] <- xi.yr[a] * eps.yr[t,a] 
            r.ranef.yr[t,a] <- xi.yr[a+3] * eps.yr[t,a+3]
            f.ranef.yr[t,a] <- xi.yr[a+6] * eps.yr[t,a+6]
        } #a
    } #t
    
    ## Prior for precision matrix
    Tau.yr[1:9, 1:9]  ~ dwish(W[1:9, 1:9], 10)
    Sigma.yr.raw[1:9, 1:9] <- inverse(Tau.yr[1:9, 1:9])
    for (k in 1:9){ # get the posterior estimation of sigma and correlation
        for (k.prime in 1:9){
            rho.yr[k,k.prime] <- Sigma.yr.raw[k,k.prime]/
                sqrt(Sigma.yr.raw[k,k]*Sigma.yr.raw[k.prime,k.prime])
        }
        Sigma.yr[k] <- abs(xi.yr[k])*sqrt(Sigma.yr.raw[k,k])
    }
    
    # [id, Z]  &  [farm,z] random effect. also from multivariate normal distribution to acount for cov between traits. but no age interaction to avoid too much over parametrization
    xi.id[1] ~ dunif(0,10)
    xi.id[2] ~ dunif(0,10)
    
    for (i in 1:(nb.id)){
        eps.id[i,1:2]~ dmnorm(zero[1:2], Tau.id[1:2, 1:2])
        r.ranef.id[i] <-  xi.id[1] *eps.id[i,1]
        f.ranef.id[i] <- xi.id[2] *eps.id[i,2]
    } #i
    Tau.id[1:2, 1:2]  ~ dwish(W[1:2, 1:2], 3)
    Sigma.id.raw[1:2, 1:2] <- inverse(Tau.id[1:2, 1:2])
    for (k in 1:2){# get the posterior estimation of sigma and correlation
        for (k.prime in 1:2){
            rho.id[k,k.prime] <- Sigma.id.raw[k,k.prime]/
                sqrt(Sigma.id.raw[k,k]*Sigma.id.raw[k.prime,k.prime])
        }
        Sigma.id[k] <- abs(xi.id[k])*sqrt(Sigma.id.raw[k,k])
    }
    
    
    for(z in 1:3){
        xi.farm[z] ~ dunif(0,4)
    }
    for (f in 1:nb.farm){
        eps.farm[f,1:3]  ~ dmnorm(zero[1:3], Tau.farm[1:3, 1:3])
        s.ranef.farm[f] <- xi.farm[1] *eps.farm[f,1]
        r.ranef.farm[f] <- xi.farm[2] *eps.farm[f,2]
        f.ranef.farm[f] <- xi.farm[3] *eps.farm[f,3]
    } #f
    Tau.farm[1:3, 1:3]  ~ dwish(W[1:3, 1:3], 4)
    Sigma.farm.raw[1:3, 1:3] <- inverse(Tau.farm[1:3, 1:3])
    for (k in 1:3){# get the posterior estimation of sigma and correlation
        for (k.prime in 1:3){
            rho.farm[k,k.prime] <- Sigma.farm.raw[k,k.prime]/
                sqrt(Sigma.farm.raw[k,k]*Sigma.farm.raw[k.prime,k.prime])
        }
        Sigma.farm[k] <- abs(xi.farm[k])*sqrt(Sigma.farm.raw[k,k])
    }
    
    # capture probs
    mu.p[1] ~ dlogis(0,1)
    mu.p[2] ~ T(dlogis(0,1),-5,10) # constrained it a bit because it sometimes went crazy high
    
    tau.p <- 1/(sd.p*sd.p)
    sd.p~ dunif(0,5) # sd for the random yr effect on capture prob
    
    
    # Farm location and movement
    sig~dunif(0,50) 
    farm.probs[1:nb.farm,1:nb.farm] <-exp(- D[1:nb.farm,1:nb.farm] / sig)  # exponential decline in mvt prob with distance. 
    
    for(i in 1:nb.id){
        for(t in (first[i]+1):nb.t){
            farm[i,t]~dcat(farm.probs[farm[i,t-1],1:nb.farm])
        }
    }
    
    
    # Define state-transition and observation matrices      -------------------
    for (i in 1:nb.id){
        # Define probabilities of state S(t+1) given S(t)
        for (t in first[i]:(nb.t-1)){
            logit(phi[i,t]) <- s.B.int[age[i,t]] +
                s.B.Env[1, age[i,t]]*  x.farmYrEnv[farm[i,t],t,1]+ # temp
                s.B.Env[2, age[i,t]]*  x.farmYrEnv[farm[i,t],t,2]+ # prec
                s.B.Env[3, age[i,t]]*  x.farmYrEnv[farm[i,t],t,3]+ # coldSnap
                s.B.Env[4, age[i,t]]*  x.farmYrEnv[farm[i,t],t,4]+ # hosp
                s.B.Env[5, age[i,t]]*  x.farmYrEnv[farm[i,t],t,5]+ # IA
                s.B.Env[6, age[i,t]]*  x.farmYrEnv[farm[i,t],t,6]+# prec:cold
                s.B.Env[7, age[i,t]]* x.farmYrEnv[farm[i,t],t,7]+# prec:IA
                s.B.Env[8, age[i,t]]* x.farmYrEnv[farm[i,t],t,8]+# IA:colf
                s.B.Env[9, age[i,t]]*  x.farmYrEnv[farm[i,t],t,9]+ # triple
                s.ranef.farm[farm[i,t]]+
                s.ranef.yr[t,age[i,t]]
            
            logit(r[i,t]) <-  r.B.int[age[i,t]] +
                r.B.Env[1, age[i,t]]* x.farmYrEnv[farm[i,t],t+1,1]+ # temp
                r.B.Env[2, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,2]+ # prec
                r.B.Env[3, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,3]+ # coldrnap
                r.B.Env[4, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,4]+ # horp
                r.B.Env[5, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,5]+ # IA
                r.B.Env[6, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,6]+ # prec:cold
                r.B.Env[7, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,7]+ # prec:IA
                r.B.Env[8, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,8]+ # IA:colf
                r.B.Env[9, age[i,t]]*  x.farmYrEnv[farm[i,t],t+1,9]+ # triple
                r.ranef.farm[farm[i,t]]+
                r.ranef.yr[t,age[i,t]]+
                r.ranef.id[i]
            
            # ps[from,i,t,to]
            ps[1,i,t,1] <- phi[i,t]*(1-r[i,t])
            ps[1,i,t,2] <- phi[i,t]*r[i,t]
            ps[1,i,t,3] <- 1-phi[i,t]
            
            ps[2,i,t,1] <- phi[i,t]*(1-r[i,t])
            ps[2,i,t,2] <- phi[i,t]*r[i,t]
            ps[2,i,t,3] <- 1-phi[i,t]
            
            ps[3,i,t,1] <- 0
            ps[3,i,t,2] <- 0
            ps[3,i,t,3] <- 1
        }
    }
    for (t in 1:(nb.t-1)){
        ranef.pa[t] ~ dnorm(0,tau.p)
        logit(p1[t]) <- mu.p[1] +ranef.pa[t]
        logit(p2[t]) <-  mu.p[2] +ranef.pa[t]
        
        po[1,t,1] <- p1[t]
        po[1,t,2] <- 0
        po[1,t,3] <- 1-p1[t]
        
        po[2,t,1] <- 0
        po[2,t,2] <- p2[t]
        po[2,t,3] <- 1-p2[t]
        
        po[3,t,1] <- 0
        po[3,t,2] <- 0
        po[3,t,3] <- 1
    } #t
    
    for(i in 1:nb.mat){
        for( t in first[iMat[i]]:(nb.t-1)){
            log(lambda[i,t]) <- f.B.int[age[iMat[i],t]] +
                f.B.Env[1, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,1]+ # temp
                f.B.Env[2, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,2]+ # prec
                f.B.Env[3, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,3]+ # coldfnap
                f.B.Env[4, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,4]+ # hofp
                f.B.Env[5, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,5]+ # IA
                f.B.Env[6, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,6]+ # prec:cold
                f.B.Env[7, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,7]+ # prec:IA
                f.B.Env[8, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,8]+ # IA:colf
                f.B.Env[9, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t],t+1,9]+ # triple
                f.ranef.farm[farm[iMat[i],t]]+
                f.ranef.yr[t,(age[iMat[i],t])]+
                f.ranef.id[iMat[i]]
        }
    }
    
    # Likelihood   -----------------------------------
    
    for (i in 1:nb.id){
        for (t in (first[i]+1):nb.t){
            # State process: draw S(t) given S(t-1)
            state[i,t] ~ dcat(ps[state[i,t-1], i, t-1,1:3])
            # Observation process: draw O(t) given S(t)
            obs[i,t] ~ dcat(po[state[i,t], t-1,1:3])
            # obs.hat[i,t] ~ dcat(po[state[i,t], t-1,1:3])  # for posterior predictive checks keep commented if not used
        } #t
    }
    # for nb of fledging produced. zero-inflated poisson
    for(i in 1:nb.mat){
        for( t in (first[iMat[i]]+1):nb.t){
            nff[i,t] <- lambda[i,t-1]* (state[iMat[i],t]==2)+0.00001
            nbFledge[iMat[i],t]~dpois(nff[i,t])
            #   nbFledge.hat[i,t] ~ dpois(nff[i,t])   # for posterior predictive checks keep commented if not used
        }
    } #i
    
    
    # Calculate derived population parameters  -------------------
    
    # get nb of id of each age class in the pop	
     for(t in 1:nb.t){
       for(i in 1:nb.id) {
           st1[i,t] <- (state[i,t]>0) * (state[i,t]<3)* age[i,t]
       }
       Nm[1,t] <- sum(st1[1:nb.id,t]==1) # nb marked ois
       Nm[2,t] <- sum(st1[1:nb.id,t]==2) # nb marked SY
       Nm[3,t] <- sum(st1[1:nb.id,t]==3) # nb marked SY
     }
    # #
    # 
    # # immigration
    # # similar to Taylor et al. 2018  & (Schaub and Fletcher 2015).
     for(t in 2:nb.t){
       imia[t] <-  round(imi[t,2] *0.48)  # *48 = prop of ASY which come from SY according to Esther
       imib[t] <- imi[t,2]-imia[t]
       Pim[1,t] <- (imi[t,1])/Nm[1,t-1]
       Pim[2,t] <- (imia[t])/Nm[2,t-1]
       Pim[3,t] <- (imib[t])/Nm[3,t-1]
     }
}
)


myInits <- function(curDat,curConst){
    l=list(
        s.B.int=rnorm(3,c(-3,0,1),0.05),
        r.B.int=rnorm(3,0, 0.05),
        f.B.int=rnorm(3,1.5, 0.05),
        s.B.Env=matrix(rnorm(9*3,0, 0.25),ncol = 3),
        r.B.Env=matrix(rnorm(9*3,0, 0.25),ncol = 3),
        f.B.Env=matrix(rnorm(9*3,0, 0.25),ncol = 3),
        mu.p= rnorm(2,c(0,3), 0.12),
        ranef.pa=rnorm(14-1,0,0.1),
        tau.p=runif(1,4,8),
        Tau.yr=inverse(diag(c(0.2,0.3,0.2,.8,.5,.5,.4,.1,.2))),
        Tau.farm=inverse(diag(c(0.2,0.5,0.1))),
        Tau.id=inverse(diag(c(1.5,0.1))),
        sig=rnorm(1,3,.5),
        xi.farm=runif(3,0.95,1.05),
        xi.id=runif(2,0.95,1.05),
        xi.yr=runif(9,0.95,1.05),
        eps.yr=matrix(rnorm((14-1)*9,0, 0.05),ncol = 9),
        farm=apply(curDat$farm,1:2,function(ff) ifelse(is.na(ff),sample(1:40,1),NA)),
        state=matrix(NA,nrow=nrow(curDat$state),ncol=ncol(curDat$state))
    )
    nbFledge.tmp <- matrix(sample(1:6,replace = T,size = curConst$nb.id*curConst$nb.t),nrow = curConst$nb.id,ncol = curConst$nb.t)
    for(i in 1:curConst$nb.id){
        t2 <- max(which(curDat$state[i,] %in% 1:2))
        l$state[i,1:t2]=sample(1:2,t2,replace = T)
        if(t2<curConst$nb.t){
            l$state[i,t2:curConst$nb.t]=sample(1:3,length(t2:curConst$nb.t),replace = T)
            if(sum(l$state[i,]==3)>0){
                t3<- min(which(l$state[i,] %in% 3))
                if(t3<=curConst$nb.t) l$state[i,t3:curConst$nb.t] <- 3
            }
        }
    }
    nbFledge.tmp <- nbFledge.tmp * l$state==2
    l$state <- ifelse(is.na(curDat$state),l$state,NA)
    l$nbFledge <- ifelse(is.na(curDat$nbFledge),nbFledge.tmp,NA)
    # l$state=ms.init.z(obs, first)
    
    return(l)
}

MyVars=c('s.B.int','r.B.int','f.B.int',
         'Sigma.id','rho.id', 'xi.id', 'eps.id',
         'Sigma.yr','rho.yr',  'xi.yr','eps.yr' ,
         'Sigma.farm', 'rho.farm','xi.farm', 'eps.farm'  ,
         "cor.id","sd.id", 'cor.yr','sd.yr',
         'cor.farm','sd.farm',
         's.B.Env','r.B.Env','f.B.Env',
         's.ranef.yr','r.ranef.yr','f.ranef.yr',
         's.ranef.farm','r.ranef.farm','f.ranef.farm',
         'sig',
         'Pim','Nm',
         # 'Notzero',
         # 'obs.hat','nbFledge.hat',
         # 'state', 'nff',
         # 'obs',
         # 'farm',
         'mu.p','p1','p2','sd.p')
