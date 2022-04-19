# CMR mutistale joints
# survival + RS  / nbFledge
# multi-state  /  zero-inflated poisson
# choleski decomposition for mvnorm rnadom effect
# both yr and farm random effect are age depedent
# no E effect

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
        s.B.Env[4,a] <- 0# ~dnorm(0,0.001)            # hosp
        s.B.Env[5,a] <- 0# ~dnorm(0,0.001)# IA
        s.B.Env[6,a] <- 0#  ~dnorm(0,0.001)# prec:cold
        s.B.Env[7,a] <- 0#  ~dnorm(0,0.001)# prec:IA
        s.B.Env[8,a] <- 0#  ~dnorm(0,0.001)# IA:colf
        s.B.Env[9,a] <- 0#  ~dnorm(0,0.001)# triple
        
        r.B.Env[1,a] <- 0#  ~dnorm(0,0.001)            # temp
        r.B.Env[2,a] <- 0#  ~dnorm(0,0.001) # prec
        r.B.Env[3,a] <- 0#  ~dnorm(0,0.001)# coldSnap
        r.B.Env[4,a] <- 0# ~dnorm(0,0.001)            # hosp
        r.B.Env[5,a] <- 0# ~dnorm(0,0.001)# IA
        r.B.Env[6,a] <- 0#  ~dnorm(0,0.001)# prec:cold
        r.B.Env[7,a] <- 0#  ~dnorm(0,0.001)# prec:IA
        r.B.Env[8,a] <- 0#  ~dnorm(0,0.001)# IA:colf
        r.B.Env[9,a] <- 0#  ~dnorm(0,0.001)# triple
        
        f.B.Env[1,a] <- 0#  ~dnorm(0,0.001)            # temp
        f.B.Env[2,a] <- 0#  ~dnorm(0,0.001) # prec
        f.B.Env[3,a] <- 0#  ~dnorm(0,0.001)# coldSnap
        f.B.Env[4,a] <- 0# ~dnorm(0,0.001)            # hosp
        f.B.Env[5,a] <- 0# ~dnorm(0,0.001)# IA
        f.B.Env[6,a] <- 0#  ~dnorm(0,0.001)# prec:cold
        f.B.Env[7,a] <- 0#  ~dnorm(0,0.001)# prec:IA
        f.B.Env[8,a] <- 0#  ~dnorm(0,0.001)# IA:cold
        f.B.Env[9,a] <- 0#  ~dnorm(0,0.001)# triple
        
    }
    
    # ######  random yr effects   -------------------------
    for (j in 1:9) {
        A[j, j] ~ T(dnorm(0.0, 0.4444444),0.0,)
        Delta[j, j] <- 1/tau[j] ; tau[j] ~ dgamma(1.5, 1.5)
        l[j, j] <- 1.0
    }
    for (j in 1:(9-1)) {
        for (k in (j+1):9) {
            l[j, k] <- 0.0; A[j, k] <- 0.0; Delta[j, k] <- 0.0
            l[k, j] ~ dnorm(0.0, 4.0); A[k, j] <- 0.0; Delta[k, j] <- 0.0
        }
    }
    Omega.yr[1:9,1:9] <- A[1:9,1:9]%*%l[1:9,1:9]%*%Delta[1:9,1:9]%*%t(l[1:9,1:9])%*%A[1:9,1:9]
    for(i in 1:nb.t){
        for(j in 1:9){
            xi[i, j] ~ dnorm(0.0, tau[j])
        }
        s.ranef.yr[i, 1] <- A[1, 1]*(l[1, 1]*xi[i, 1]);
        s.ranef.yr[i, 2] <- A[2, 2]*(l[2, 1]*xi[i, 1] +
                            l[2, 2]*xi[i, 2])
        s.ranef.yr[i, 3] <- A[3, 3]*(l[3, 1]*xi[i, 1] +
                                    l[3, 2]*xi[i, 2] +
                                    l[3, 3]*xi[i, 3])
        r.ranef.yr[i, 1] <- A[4, 4]*(l[4, 1]*xi[i, 1] +
                                l[5, 2]*xi[i, 2] +
                                l[5, 4]*xi[i, 3] +
                                l[5, 5]*xi[i, 4])
        r.ranef.yr[i, 2] <- A[5, 5]*(l[5, 1]*xi[i, 1] +
                                l[5, 2]*xi[i, 2] +
                                l[5, 3]*xi[i, 3] +
                                l[5, 4]*xi[i, 4] +
                                l[5, 5]*xi[i, 5])
        r.ranef.yr[i, 3] <- A[6, 6]*(l[6, 1]*xi[i, 1] +
                                l[6, 2]*xi[i, 2] +
                                l[6, 3]*xi[i, 3] +
                                l[6, 4]*xi[i, 4] +
                                l[6, 5]*xi[i, 5] +
                                l[6, 6]*xi[i, 6])
        f.ranef.yr[i, 1] <- A[7, 7]*(l[7, 1]*xi[i, 1] +
                                l[7, 2]*xi[i, 2] +
                                l[7, 3]*xi[i, 3] +
                                l[7, 4]*xi[i, 4] +
                                l[7, 5]*xi[i, 5] +
                                l[7, 6]*xi[i, 6] +
                                l[7, 7]*xi[i, 7])
        f.ranef.yr[i, 2] <- A[8, 8]*(l[8, 1]*xi[i, 1] +
                                l[8, 2]*xi[i, 2] +
                                l[8, 3]*xi[i, 3] +
                                l[8, 4]*xi[i, 4] +
                                l[8, 5]*xi[i, 5] +
                                l[8, 6]*xi[i, 6] +
                                l[8, 7]*xi[i, 7] +
                                l[8, 8]*xi[i, 8])
        f.ranef.yr[i, 3] <- A[9, 9]*(l[9, 1]*xi[i, 1] +
                                l[9, 2]*xi[i, 2] +
                                l[9, 3]*xi[i, 3] +
                                l[9, 4]*xi[i, 4] +
                                l[9, 5]*xi[i, 5] +
                                l[9, 6]*xi[i, 6] +
                                l[9, 7]*xi[i, 7] +
                                l[9, 8]*xi[i, 8] +
                                l[9, 9]*xi[i, 9])
    }
    
    # random age:farm effect 
    for (j in 1:9) {
        A.farm[j, j] ~ T(dnorm(0.0, 0.4444444),0.0,)
        Delta.farm[j, j] <- 1/tau.farm[j] ; tau.farm[j] ~ dgamma(1.5, 1.5)
        l.farm[j, j] <- 1.0
    }
    for (j in 1:(9-1)) {
        for (k in (j+1):9) {
            l.farm[j, k] <- 0.0; A.farm[j, k] <- 0.0; Delta.farm[j, k] <- 0.0
            l.farm[k, j] ~ dnorm(0.0, 4.0); A.farm[k, j] <- 0.0; Delta.farm[k, j] <- 0.0
        }
    }
    Omega.farm[1:9,1:9] <- A.farm[1:9,1:9]%*%l.farm[1:9,1:9]%*%Delta.farm[1:9,1:9]%*%t(l.farm[1:9,1:9])%*%A.farm[1:9,1:9]
    for(i in 1:nb.farm){
        for(j in 1:9){
            xi.farm[i, j] ~ dnorm(0.0, tau.farm[j])
        }
        s.ranef.farm[i, 1] <- A.farm[1, 1]*(l.farm[1, 1]*xi.farm[i, 1]);
        s.ranef.farm[i, 2] <- A.farm[2, 2]*(l.farm[2, 1]*xi.farm[i, 1] +
                                         l.farm[2, 2]*xi.farm[i, 2])
        s.ranef.farm[i, 3] <- A.farm[3, 3]*(l.farm[3, 1]*xi.farm[i, 1] +
                                         l.farm[3, 2]*xi.farm[i, 2] +
                                         l.farm[3, 3]*xi.farm[i, 3])
        r.ranef.farm[i, 1] <- A.farm[4, 4]*(l.farm[4, 1]*xi.farm[i, 1] +
                                         l.farm[5, 2]*xi.farm[i, 2] +
                                         l.farm[5, 4]*xi.farm[i, 3] +
                                         l.farm[5, 5]*xi.farm[i, 4])
        r.ranef.farm[i, 2] <- A.farm[5, 5]*(l.farm[5, 1]*xi.farm[i, 1] +
                                         l.farm[5, 2]*xi.farm[i, 2] +
                                         l.farm[5, 3]*xi.farm[i, 3] +
                                         l.farm[5, 4]*xi.farm[i, 4] +
                                         l.farm[5, 5]*xi.farm[i, 5])
        r.ranef.farm[i, 3] <- A.farm[6, 6]*(l.farm[6, 1]*xi.farm[i, 1] +
                                         l.farm[6, 2]*xi.farm[i, 2] +
                                         l.farm[6, 3]*xi.farm[i, 3] +
                                         l.farm[6, 4]*xi.farm[i, 4] +
                                         l.farm[6, 5]*xi.farm[i, 5] +
                                         l.farm[6, 6]*xi.farm[i, 6])
        f.ranef.farm[i, 1] <- A.farm[7, 7]*(l.farm[7, 1]*xi.farm[i, 1] +
                                         l.farm[7, 2]*xi.farm[i, 2] +
                                         l.farm[7, 3]*xi.farm[i, 3] +
                                         l.farm[7, 4]*xi.farm[i, 4] +
                                         l.farm[7, 5]*xi.farm[i, 5] +
                                         l.farm[7, 6]*xi.farm[i, 6] +
                                         l.farm[7, 7]*xi.farm[i, 7])
        f.ranef.farm[i, 2] <- A.farm[8, 8]*(l.farm[8, 1]*xi.farm[i, 1] +
                                         l.farm[8, 2]*xi.farm[i, 2] +
                                         l.farm[8, 3]*xi.farm[i, 3] +
                                         l.farm[8, 4]*xi.farm[i, 4] +
                                         l.farm[8, 5]*xi.farm[i, 5] +
                                         l.farm[8, 6]*xi.farm[i, 6] +
                                         l.farm[8, 7]*xi.farm[i, 7] +
                                         l.farm[8, 8]*xi.farm[i, 8])
        f.ranef.farm[i, 3] <- A.farm[9, 9]*(l.farm[9, 1]*xi.farm[i, 1] +
                                         l.farm[9, 2]*xi.farm[i, 2] +
                                         l.farm[9, 3]*xi.farm[i, 3] +
                                         l.farm[9, 4]*xi.farm[i, 4] +
                                         l.farm[9, 5]*xi.farm[i, 5] +
                                         l.farm[9, 6]*xi.farm[i, 6] +
                                         l.farm[9, 7]*xi.farm[i, 7] +
                                         l.farm[9, 8]*xi.farm[i, 8] +
                                         l.farm[9, 9]*xi.farm[i, 9])
    }
    
    # ######  random effects id 
    for (j in 1:2) {
        A.id[j, j] ~ T(dnorm(0.0, 0.4444444),0.0,)
        Delta.id[j,j] <- 1/tau.id[j] ; tau.id[j] ~ dgamma(1.5, 1.5);
        l.id[j, j] <- 1.0;
    }
    l.id[1, 2] <- 0.0; A.id[1, 2] <- 0.0; Delta.id[1, 2] <- 0.0;
    l.id[2, 1] ~ dnorm(0.0, 4.0); A.id[2, 1] <- 0.0; Delta.id[2, 1] <- 0.0;
    
    Omega.id[1:2,1:2] <- A[1:2,1:2]%*%l[1:2,1:2]%*%Delta[1:2,1:2]%*%t(l[1:2,1:2])%*%A[1:2,1:2]
    for(i in 1:nb.id){
        for(j in 1:2){
            xi.id[i, j] ~ dnorm(0.0, tau.id[j])
        }
        r.ranef.id[i] <- A.id[1, 1]*(l.id[1, 1]*xi.id[i, 1]);
        f.ranef.id[i] <- A.id[2, 2]*(l.id[2, 1]*xi.id[i, 1] +
                                         l[2, 2]*xi.id[i, 2])
    }
    
    
    s.cost~dnorm(0,0.001)
    r.cost~dnorm(0,0.001)
    
    # capture probs
    mu.p[1] ~ dlogis(0,1)
    mu.p[2] ~ T(dlogis(0,1),-5,10) # constrained it a bit because it sometimes went crazy high
    
    tau.p <- 1/(sd.p*sd.p)
    sd.p~ dunif(0,5) # sd for the random yr effect on capture prob
    
    
    # Farm location and movement
    # sig[1]~dunif(0,50) 
    # sig[2]~dunif(0,50) 
    # sig[3] <- 0.1
    # for(s in 1:3){
    #     farm.probs[1:nb.farm,1:nb.farm,s] <-exp(- D[1:nb.farm,1:nb.farm] / sig[s])  # exponential decline in mvt prob with distance. 
    # }
    
    sig~dunif(0,50) 
    farm.probs[1:nb.farm,1:nb.farm] <-exp(- D[1:nb.farm,1:nb.farm] / sig)  # exponential decline in mvt prob with distance. 
    for(i in 1:nb.id){
        for(t in (first[i]+1):nb.t){
            # farm[i,t]~dcat(farm.probs[farm[i,t-1],1:nb.farm,state[i,t-1]]) # for state dependant dispersion
            farm[i,t]~dcat(farm.probs[farm[i,t-1],1:nb.farm])
        }
    }
    
    
    # Define state-transition and observation matrices      -------------------
    for (i in 1:nb.id){
        # Define probabilities of state S(t+1) given S(t)
        for (t in first[i]:(nb.t-1)){
            logit(phi1[i,t]) <- s.B.int[age[i,t]] +
                s.B.Env[1, age[i,t]]*  x.farmYrEnv[farm[i,t],t,1]+ # temp
                s.B.Env[2, age[i,t]]*  x.farmYrEnv[farm[i,t],t,2]+ # prec
                s.B.Env[3, age[i,t]]*  x.farmYrEnv[farm[i,t],t,3]+ # coldSnap
                s.B.Env[4, age[i,t]]*  x.farmYrEnv[farm[i,t],t,4]+ # hosp
                s.B.Env[5, age[i,t]]*  x.farmYrEnv[farm[i,t],t,5]+ # IA
                s.B.Env[6, age[i,t]]*  x.farmYrEnv[farm[i,t],t,6]+# prec:cold
                s.B.Env[7, age[i,t]]* x.farmYrEnv[farm[i,t],t,7]+# prec:IA
                s.B.Env[8, age[i,t]]* x.farmYrEnv[farm[i,t],t,8]+# IA:colf
                s.B.Env[9, age[i,t]]*  x.farmYrEnv[farm[i,t],t,9]+ # triple
                s.ranef.farm[farm[i,t],age[i,t]]+
                s.ranef.yr[t,age[i,t]]
            
            logit(phi2[i,t]) <- s.B.int[age[i,t]] + s.cost +
                s.B.Env[1, age[i,t]]*  x.farmYrEnv[farm[i,t],t,1]+ # temp
                s.B.Env[2, age[i,t]]*  x.farmYrEnv[farm[i,t],t,2]+ # prec
                s.B.Env[3, age[i,t]]*  x.farmYrEnv[farm[i,t],t,3]+ # coldSnap
                s.B.Env[4, age[i,t]]*  x.farmYrEnv[farm[i,t],t,4]+ # hosp
                s.B.Env[5, age[i,t]]*  x.farmYrEnv[farm[i,t],t,5]+ # IA
                s.B.Env[6, age[i,t]]*  x.farmYrEnv[farm[i,t],t,6]+# prec:cold
                s.B.Env[7, age[i,t]]* x.farmYrEnv[farm[i,t],t,7]+# prec:IA
                s.B.Env[8, age[i,t]]* x.farmYrEnv[farm[i,t],t,8]+# IA:colf
                s.B.Env[9, age[i,t]]*  x.farmYrEnv[farm[i,t],t,9]+ # triple
                s.ranef.farm[farm[i,t],age[i,t]]+
                s.ranef.yr[t,age[i,t]]

            
            logit(r1[i,t]) <-  r.B.int[age[i,t]] +
                r.B.Env[1, age[i,t]]* x.farmYrEnv[farm[i,t+1],t+1,1]+ # temp
                r.B.Env[2, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,2]+ # prec
                r.B.Env[3, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,3]+ # coldrnap
                r.B.Env[4, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,4]+ # horp
                r.B.Env[5, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,5]+ # IA
                r.B.Env[6, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,6]+ # prec:cold
                r.B.Env[7, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,7]+ # prec:IA
                r.B.Env[8, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,8]+ # IA:colf
                r.B.Env[9, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,9]+ # triple
                r.ranef.farm[farm[i,t+1],age[i,t]]+
                r.ranef.yr[t+1,age[i,t]]+
                r.ranef.id[i]

            logit(r2[i,t]) <-  r.B.int[age[i,t]] + r.cost +
                r.B.Env[1, age[i,t]]* x.farmYrEnv[farm[i,t+1],t+1,1]+ # temp
                r.B.Env[2, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,2]+ # prec
                r.B.Env[3, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,3]+ # coldrnap
                r.B.Env[4, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,4]+ # horp
                r.B.Env[5, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,5]+ # IA
                r.B.Env[6, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,6]+ # prec:cold
                r.B.Env[7, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,7]+ # prec:IA
                r.B.Env[8, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,8]+ # IA:colf
                r.B.Env[9, age[i,t]]*  x.farmYrEnv[farm[i,t+1],t+1,9]+ # triple
                r.ranef.farm[farm[i,t+1],age[i,t]]+
                r.ranef.yr[t+1,age[i,t]]+
                r.ranef.id[i]
            
            # ps[from,i,t,to]
            ps[1,i,t,1] <- phi1[i,t]*(1-r1[i,t])
            ps[1,i,t,2] <- phi1[i,t]*r1[i,t]
            ps[1,i,t,3] <- 1-phi1[i,t]
            
            ps[2,i,t,1] <- phi2[i,t]*(1-r2[i,t])
            ps[2,i,t,2] <- phi2[i,t]*r2[i,t]
            ps[2,i,t,3] <- 1-phi2[i,t]
            
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
                f.B.Env[1, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,1]+ # temp
                f.B.Env[2, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,2]+ # prec
                f.B.Env[3, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,3]+ # coldfnap
                f.B.Env[4, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,4]+ # hofp
                f.B.Env[5, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,5]+ # IA
                f.B.Env[6, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,6]+ # prec:cold
                f.B.Env[7, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,7]+ # prec:IA
                f.B.Env[8, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,8]+ # IA:colf
                f.B.Env[9, age[iMat[i],t]]*  x.farmYrEnv[farm[iMat[i],t+1],t+1,9]+ # triple
                f.ranef.farm[farm[iMat[i],t+1],age[i,t]]+
                f.ranef.yr[t+1,(age[iMat[i],t])]+
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
        r.cost=rnorm(1,0,0.5),s.cost=rnorm(1,0,0.5),
        mu.p= rnorm(2,c(0,3), 0.12),
        ranef.pa=rnorm(14-1,0,0.1),
        xi=matrix(0.01,nrow = curConst$nb.t,ncol = 9),
        xi.farm=matrix(0.01,nrow = 40,ncol = 9),
        xi.id=matrix(0.01,nrow = curConst$nb.id,ncol = 2),
        tau.p=runif(1,4,8),
        sig=rnorm(1,3,.5),
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
         's.B.Env','r.B.Env','f.B.Env',
         'Omega.id','Omega.farm','Omega.yr',
         's.ranef.yr','r.ranef.yr','f.ranef.yr',
         's.ranef.farm','r.ranef.farm','f.ranef.farm',
         'sig','s.cost','r.cost',
         'Pim','Nm',
         # 'xi', 'A','Delta'   , 'A.farm','Delta.farm','xi.farm','A.id','Delta.id','xi.id',
         'mu.p','p1','p2','sd.p')




## testing ---------------------

library(coda)
library(tidyverse)
source('R/999_MyFunc.R')
load('cache/cleanMultiState.Rdata')
#
start <- now()
nimbleOut <- nimbleMCMC(myCode,constants = miniConst,data = miniDat,
           # niter = 2000,nburnin = 1000,nchains = 3,
            # niter = (50000+1000*5),nburnin = 50000,nchains = 3,
           niter = (50000+1000*1),nburnin = 2000,nchains = 2,
           monitors = MyVars,
           summary = T,samplesAsCodaMCMC = T,
           inits = myInits(curDat = miniDat,curConst = miniConst)
        )
dur=now()
save(nimbleOut,dur,file = 'cache/v4CostMini1c.Rdata')

# 
# 
# plot(nimbleOut$samples[,grepl(".B.int",colnames(nimbleOut$samples[[1]]))])
# plot(nimbleOut$samples[,"Omega.yr[3, 3]"])
# plot(nimbleOut$samples[,"Omega.yr[3, 2]"])
# plot(nimbleOut$samples[,"sig"])
# plot(nimbleOut$samples[,"mu.p[1]"])
# plot(nimbleOut$samples[,"mu.p[2]"])
# 
nimbleOut$summary[,3] %>% hist

