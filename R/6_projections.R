library(tidyverse)
library(coda)
library(popbio)
library(boot)
library(scales)
library(corrplot)

source('R/999_MyFunc.R')
# source('R/110_cleanMultiState.R')
load('cache/cleanMultiState.Rdata')

load('../14_Trsw_CMR/cache/313_m6f_0916.Rdata')
betazzP <- chain_output %>% #map(~.x[[1]]) %>%
  map(~as.data.frame(.x)) %>% bind_rows() %>%
  select(contains('Pim'))

load('cache/out_v0.Rdata')
betazz <- chain_output %>% map(~.x[[1]]) %>%
  map(~as.data.frame(.x[[1]])) %>% bind_rows() %>%
  select(!contains('hat') & !contains('nff') )

ltmp=max(nrow(betazz),nrow(betazzP))
betazz <- cbind(betazz[1:ltmp,],betazzP[1:ltmp,])



mclist <- chain_output %>% map(~as.mcmc(.x[[1]][[1]]) )%>% as.mcmc.list()
# mclist <- chain_output %>% map(~as.mcmc(.x) )%>% as.mcmc.list()

tmp <- summary(mclist[,grepl('.B.Env',colnames(mclist[[1]]))])
tmp <- tibble(var=row.names(tmp[[1]]),
              mean=tmp[[1]][,'Mean'],
              sd=tmp[[1]][,'SD'],
              cil=tmp[[2]][,1],
              cih=tmp[[2]][,5])
tmp %>% filter(sd>0) %>% filter((cil*cih)>0)
tmp %>% filter(sd>0)


plot(mclist[,'s.B.Env[5, 2]'])
plot(mclist[,'s.B.Env[5, 1]'])
plot(mclist[,'s.B.Env[4, 2]'])
plot(mclist[,'s.B.int[2]'])

plot(mclist[,'mu.p[1]'])
plot(mclist[,'mu.p[2]'])

rm(betazzP,chain_output,mclist)

betazz %>% select(contains('.B.')) %>%
  select(., where(function(x){sd(x)>0} ))%>%
  cor %>% corrplot()


nb1=colSums((mydat$obs<3)*(myconst$age==1))
nb2=colSums((mydat$obs<3)*(myconst$age==2))
nb3=colSums((mydat$obs<3)*(myconst$age==3))
nb1['2019'] <- round(mean(nb1/(nb1+nb2+nb3))* (nb2['2019']+nb3['2019']))

obsPim <- matrix(NA,14,3)
for(t in 2:14){
  obsPim[t,1] <- mydat$imi[t,1]/nb1[t-1]
  obsPim[t,2] <- (mydat$imi[t,2]*.48)/nb2[t-1]
  obsPim[t,3] <- (mydat$imi[t,2]*.52)/nb3[t-1]
}

yrMeanEnvs <- apply(mydat$x.farmYrEnv,c(2:3) , mean) %>% as_tibble() %>%
  rename(prec.cold=`prec:cold`,
         prec.IA=`prec:IA`,
         IA.cold=`IA:cold`)

g.Pred.LMat <- function(newd,betaz=betazz,t=NA,f=NA,pim=NA){
  #     newd = 0;betaz = betazz; t = 3 ; f=13; pim=NA ;it=1:5

  it= betaz %>% as.data.frame() %>% as.matrix() %>% is.na() %>% rowSums()
  it= as.numeric(which(it==0))
  nb.t=13
  
  if(is.na(t)) t <- 999;  if(is.na(f)) f <- 999;if(is.na(pim)) pim <- 999

  if(sum(grepl('.B.Env',colnames(betaz)))<1){
    be=matrix(0,ncol = 9*3*3,nrow=nrow(betaz))
    colnames(be) <- c(paste0('s.B.Env[',rep(1:9,each=3),', ',1:3,']'),
                      paste0('r.B.Env[',rep(1:9,each=3),', ',1:3,']'),
                      paste0('f.B.Env[',rep(1:9,each=3),', ',1:3,']'))
    betaz=cbind(betaz,be)
  }

  if (pim>0){
      Pim1 <- betaz[,paste0("Pim[1, ",1:14,"]")]
      Pim2 <- betaz[,paste0("Pim[2, ",1:14,"]")]
      Pim3 <- betaz[,paste0("Pim[3, ",1:14,"]")]
      Pim1[1,1] <-Pim2[1,1] <-Pim3[1,1] <- 0   # was NA for some resaon. probably error in inits
  }


  s.ranef.yr<- array(0,c(3,length(it)))
  r.ranef.yr<-array(0,c(3,length(it)))
  f.ranef.yr<-array(0,c(3,length(it)))
  s.ranef.farm<- rep(0,c(length(it)))
  r.ranef.farm<-rep(0,c(length(it)))
  f.ranef.farm<-rep(0,c(length(it)))

  pim.ranef <-array(0,c(3,length(it)))



  if(t %in% (1:(nb.t))) {
    # for (a in 1:3){
    # s.ranef.yr[a,] <- betaz[it,paste0('s.ranef.yr[',t,', ',a,']')]
    # r.ranef.yr[a,] <- betaz[it,paste0('r.ranef.yr[',t,', ',a,']')]
    # f.ranef.yr[a,] <- betaz[it,paste0('f.ranef.yr[',t,', ',a,']')]    # ~dnorm(0,f1.tau.yr)
    # } #i
      s.ranef.yr[1,] <- betaz$`xi.yr[1]`[it] * betaz[[paste0('eps.yr[',t,', 1]')]][it]
      s.ranef.yr[2,] <- betaz$`xi.yr[2]`[it] * betaz[[paste0('eps.yr[',t,', 2]')]][it]
      s.ranef.yr[3,] <- betaz$`xi.yr[3]`[it] * betaz[[paste0('eps.yr[',t,', 3]')]][it]
      r.ranef.yr[1,] <- betaz$`xi.yr[4]`[it] * betaz[[paste0('eps.yr[',t,', 4]')]][it]
      r.ranef.yr[2,] <- betaz$`xi.yr[5]`[it] * betaz[[paste0('eps.yr[',t,', 5]')]][it]
      r.ranef.yr[3,] <- betaz$`xi.yr[6]`[it] * betaz[[paste0('eps.yr[',t,', 6]')]][it]
      f.ranef.yr[1,] <- betaz$`xi.yr[7]`[it] * betaz[[paste0('eps.yr[',t,', 7]')]][it]
      f.ranef.yr[2,] <- betaz$`xi.yr[8]`[it] * betaz[[paste0('eps.yr[',t,', 8]')]][it]
      f.ranef.yr[3,] <- betaz$`xi.yr[9]`[it] * betaz[[paste0('eps.yr[',t,', 9]')]][it]
  }
  if(t=='ran'){
    s.ranef.yr<- array(0,c(3,length(it)))
    r.ranef.yr<-array(0,c(3,length(it)))
    f.ranef.yr<-array(0,c(3,length(it)))
    eps.yr <- sapply(it,function(i){
      corr.farm <- matrix(as.numeric(betaz[i,grepl('rho.yr',colnames(betaz))]),9,9)
      # corr.farm <-diag(9)
      sd.farm <- as.numeric(betaz[i,grepl('Sigma.yr',colnames(betaz))])
      if(length(sd.farm)==0) sd.farm=rep(0,ncol(corr.farm))
      VarMat <- lavaan::cor2cov(R = corr.farm,sds = sd.farm)
      eps.farm <- MASS::mvrnorm(1,rep(0,length(sd.farm)),VarMat)
      return(eps.farm)
    })
    s.ranef.yr[1,] <- betaz$`xi.yr[1]`[it] *eps.yr[1,]
    s.ranef.yr[2,] <- betaz$`xi.yr[2]`[it] *eps.yr[2,]
    s.ranef.yr[3,] <- betaz$`xi.yr[3]`[it] *eps.yr[3,]
    r.ranef.yr[1,] <- betaz$`xi.yr[4]`[it] *eps.yr[4,]
    r.ranef.yr[2,] <- betaz$`xi.yr[5]`[it] *eps.yr[5,]
    r.ranef.yr[3,] <- betaz$`xi.yr[6]`[it] *eps.yr[6,]
    f.ranef.yr[1,] <- betaz$`xi.yr[7]`[it] *eps.yr[7,]
    f.ranef.yr[2,] <- betaz$`xi.yr[8]`[it] *eps.yr[8,]
    f.ranef.yr[3,] <- betaz$`xi.yr[9]`[it] *eps.yr[9,]
  }

  if(pim %in% 1:(nb.t)){
    pim.ranef[1,] <- (Pim1[it,pim+1])
    pim.ranef[2,] <- (Pim2[it,pim+1])
    pim.ranef[3,] <- (Pim3[it,pim+1])
  }
  if(pim>nb.t){
    pim.ranef[1,] <- rowMeans(Pim1[,-1],na.rm = T)[it]
    pim.ranef[2,] <- rowMeans(Pim2[,-1],na.rm = T)[it]
    pim.ranef[3,] <- rowMeans(Pim3[,-1],na.rm = T)[it]
  }else{
    pim.ranef[1,] <- 0
    pim.ranef[2,] <- 0
    pim.ranef[3,] <- 0
  }

  if(f %in% 1:40){
    s.ranef.farm <- betaz$`xi.farm[1]`[it] *betaz[[paste0('eps.farm[',f,', 1]')]][it]
    r.ranef.farm <- betaz$`xi.farm[2]`[it] *betaz[[paste0('eps.farm[',f,', 2]')]][it]
    f.ranef.farm <- betaz$`xi.farm[3]`[it] *betaz[[paste0('eps.farm[',f,', 3]')]][it]
  } else{
    s.ranef.farm=0
    r.ranef.farm=0
    f.ranef.farm=0
  }

  if(f=='ran'){
    eps.farm <- sapply(it,function(i){
      corr.farm <- matrix(as.numeric(betaz[i,grepl('rho.farm',colnames(betaz))]),3,3)
      sd.farm <- as.numeric(betaz[i,grepl('Sigma.farm',colnames(betaz))])
      if(length(sd.farm)==0) sd.farm=rep(0,ncol(corr.farm))

      VarMat <- lavaan::cor2cov(R = corr.farm,sds = sd.farm)
      eps.farm <- MASS::mvrnorm(1,rep(0,length(sd.farm)),VarMat)
      return(eps.farm)
    })
    s.ranef.farm <- betaz$`xi.farm[1]`[it] *eps.farm[1,]
    r.ranef.farm <- betaz$`xi.farm[1]`[it] *eps.farm[2,]
    f.ranef.farm <- betaz$`xi.farm[1]`[it] *eps.farm[3,]

  }



  mu.s1 <- betaz[it,'s.B.int[1]']+
    betaz[it,"s.B.Env[1, 1]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"s.B.Env[2, 1]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"s.B.Env[3, 1]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"s.B.Env[4, 1]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"s.B.Env[5, 1]" ]*  newd$IA[1]+ # IA
    betaz[it,"s.B.Env[6, 1]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"s.B.Env[7, 1]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"s.B.Env[8, 1]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"s.B.Env[9, 1]" ]*  newd$triple[1] # triple

  mu.s2 <- betaz[it,'s.B.int[2]']+
    betaz[it,"s.B.Env[1, 2]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"s.B.Env[2, 2]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"s.B.Env[3, 2]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"s.B.Env[4, 2]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"s.B.Env[5, 2]" ]*  newd$IA[1]+ # IA
    betaz[it,"s.B.Env[6, 2]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"s.B.Env[7, 2]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"s.B.Env[8, 2]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"s.B.Env[9, 2]" ]*  newd$triple[1] # triple

  mu.s3 <- betaz[it,'s.B.int[3]']+
    betaz[it,"s.B.Env[1, 3]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"s.B.Env[2, 3]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"s.B.Env[3, 3]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"s.B.Env[4, 3]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"s.B.Env[5, 3]" ]*  newd$IA[1]+ # IA
    betaz[it,"s.B.Env[6, 3]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"s.B.Env[7, 3]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"s.B.Env[8, 3]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"s.B.Env[9, 3]" ]*  newd$triple[1] # triple

  mu.r1 <- betaz[it,'r.B.int[1]']+
    betaz[it,"r.B.Env[1, 1]" ]*  newd$meanTemp.rearing[t]+
    betaz[it,"r.B.Env[2, 1]" ]*  newd$precip.rearing[t]+ # prec
    betaz[it,"r.B.Env[3, 1]" ]*  newd$coldSnap.rearing[t]+ # coldfnap
    betaz[it,"r.B.Env[4, 1]" ]*  newd$hosp[t]+ # hofp
    betaz[it,"r.B.Env[5, 1]" ]*  newd$IA[t]+ # IA
    betaz[it,"r.B.Env[6, 1]" ]*  newd$prec.cold[t]+ # prec:cold
    betaz[it,"r.B.Env[7, 1]" ]*  newd$prec.IA[t]+ # prec:IA
    betaz[it,"r.B.Env[8, 1]" ]*  newd$IA.cold[t]+ # IA:colf
    betaz[it,"r.B.Env[9, 1]" ]*  newd$triple[t] # triple

  mu.r2 <- betaz[it,'r.B.int[2]']+
    betaz[it,"r.B.Env[1, 2]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"r.B.Env[2, 2]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"r.B.Env[3, 2]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"r.B.Env[4, 2]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"r.B.Env[5, 2]" ]*  newd$IA[1]+ # IA
    betaz[it,"r.B.Env[6, 2]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"r.B.Env[7, 2]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"r.B.Env[8, 2]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"r.B.Env[9, 2]" ]*  newd$triple[1] # triple

  mu.r3 <- betaz[it,'r.B.int[3]']+
    betaz[it,"r.B.Env[1, 3]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"r.B.Env[2, 3]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"r.B.Env[3, 3]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"r.B.Env[4, 3]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"r.B.Env[5, 3]" ]*  newd$IA[1]+ # IA
    betaz[it,"r.B.Env[6, 3]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"r.B.Env[7, 3]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"r.B.Env[8, 3]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"r.B.Env[9, 3]" ]*  newd$triple[1] # triple

  mu.f1 <-  betaz[it,'f.B.int[1]']+
    betaz[it,"f.B.Env[1, 1]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"f.B.Env[2, 1]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"f.B.Env[3, 1]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"f.B.Env[4, 1]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"f.B.Env[5, 1]" ]*  newd$IA[1]+ # IA
    betaz[it,"f.B.Env[6, 1]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"f.B.Env[7, 1]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"f.B.Env[8, 1]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"f.B.Env[9, 1]" ]*  newd$triple[1] # triple
  mu.f2 <- betaz[it,'f.B.int[2]']+
    betaz[it,"f.B.Env[1, 2]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"f.B.Env[2, 2]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"f.B.Env[3, 2]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"f.B.Env[4, 2]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"f.B.Env[5, 2]" ]*  newd$IA[1]+ # IA
    betaz[it,"f.B.Env[6, 2]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"f.B.Env[7, 2]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"f.B.Env[8, 2]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"f.B.Env[9, 2]" ]*  newd$triple[1] # triple

  mu.f3 <- betaz[it,'f.B.int[3]']+
    betaz[it,"f.B.Env[1, 3]" ]*  newd$meanTemp.rearing[1]+
    betaz[it,"f.B.Env[2, 3]" ]*  newd$precip.rearing[1]+ # prec
    betaz[it,"f.B.Env[3, 3]" ]*  newd$coldSnap.rearing[1]+ # coldfnap
    betaz[it,"f.B.Env[4, 3]" ]*  newd$hosp[1]+ # hofp
    betaz[it,"f.B.Env[5, 3]" ]*  newd$IA[1]+ # IA
    betaz[it,"f.B.Env[6, 3]" ]*  newd$prec.cold[1]+ # prec:cold
    betaz[it,"f.B.Env[7, 3]" ]*  newd$prec.IA[1]+ # prec:IA
    betaz[it,"f.B.Env[8, 3]" ]*  newd$IA.cold[1]+ # IA:colf
    betaz[it,"f.B.Env[9, 3]" ]*  newd$triple[1] # triple


  s1 <- inv.logit(mu.s1 + s.ranef.yr[1,] + s.ranef.farm )
  s2 <- inv.logit(mu.s2 + s.ranef.yr[2,] + s.ranef.farm)
  s3 <- inv.logit(mu.s3 + s.ranef.yr[3,] + s.ranef.farm)
  r1 <- inv.logit(mu.r2 + r.ranef.yr[1,] + r.ranef.farm)
  r2 <- inv.logit(mu.r2 + r.ranef.yr[2,] + r.ranef.farm)
  r3 <- inv.logit(mu.r3 + r.ranef.yr[3,] + r.ranef.farm)
  f1 <- exp(mu.f1+ f.ranef.yr[1,]+ f.ranef.farm)
  f2 <- exp(mu.f2+ f.ranef.yr[2,]+ f.ranef.farm)
  f3 <- exp(mu.f3+ f.ranef.yr[3,]+ f.ranef.farm)

  i1 <-    pim.ranef[1,]
  i2 <-   pim.ranef[2,]
  i3 <-     pim.ranef[3,]
  # i1 <-  rep(obsPim[t+1,1],length(it))
  # i2 <-   rep(obsPim[t+1,2] ,length(it))
  # i3 <-   rep(obsPim[t+1,3] ,length(it))

  #post breeding transition matrix
  LMatrix <- sapply(1:length(it), function(i){
    LMatrix <- matrix(0,3,3)
    LMatrix[2,1] <-  s1[i] + (i1[i]) # ois -> SY
    LMatrix[3,2] <-  s2[i] + (i2[i])  # SY -> ASY
    LMatrix[3,3] <-  s3[i] + (i3[i])  # ASY -> ASY

    LMatrix[1,1] <- s1[i] * r1[i] * f1[i] * 0.5
    LMatrix[1,2] <- s2[i] * r2[i] * f2[i] * 0.5
    LMatrix[1,3] <- s3[i] * r3[i] * f3[i] * 0.5

    return(LMatrix)
  },simplify = 'array')

  tibble(s1,s2,s3,
         r1,r2,r3,
         f1,f2,f3,
         i1,i2,i3)
  out <- list(LMatrix=LMatrix,
              VR= tibble(s1,s2,s3,
                         r1,r2,r3,
                         f1,f2,f3,
                         i1,i2,i3))

  return(out)

}


# plot pred hosp ----------------------------------------------------------
newd <- yrMeanEnvs[1,] %>% mutate_all(~0) %>% slice_sample(n=11,replace = T) %>%
  mutate(hosp=unique(as.numeric(mydat$x.farmYrEnv[,,'hosp']))) %>% arrange(hosp)


allYrMat <- lapply(1:nrow(newd),function(t){
  g.Pred.LMat(newd = newd[t,],betaz = betazz,t=NA,f=NA,pim=NA)
})

meanMatz <- allYrMat %>% map(~ apply(.x[[1]],c(1,2),mean))
lbds <- map_dfr(1:nrow(newd),function(i){
  tibble(hosp=newd$hosp[i],
         lambda=apply(allYrMat[[i]]$LMatrix, 3, function(it) eigen.analysis(it)$lambda1),
         r3=allYrMat[[i]]$VR$r3,
         s3=allYrMat[[i]]$VR$s3,
         f3=allYrMat[[i]]$VR$f3,

         hosppc=(i-1)/10*100)
})

meanlambds_hosp <- allYrMat %>% map_dfr(function(t){
  lbdz <- apply(t$LMatrix, 3, function(it) eigen.analysis(it)$lambda1)

  tibble(lambda=median(lbdz),
         lbd.q25=quantile(lbdz,0.25),
         lbd.q75=quantile(lbdz,0.75),
         lbd.cil=quantile(lbdz,0.025),
         lbd.cih=quantile(lbdz,0.975))
}) %>% cbind(.,newd)

ggplot(data=meanlambds_hosp,aes(x = as.factor(hosp),y=lambda,ymin=lbd.cil,ymax=lbd.cih))+
  geom_pointrange(aes(ymin=lbd.cil,ymax=lbd.cih))+
  # geom_path()+
  scale_x_discrete(labels=paste0(seq(0,100,by=10),'%'))+
  labs(y=expression(paste("Taux de croissance assymptotique (",lambda,")")),
       x="Taux d'occupation par le moineau")+
  theme_gab() +  ylim(0.43,.95)


hospD <- tibble(hosp=0:10,
                n=as.numeric(table(yrFarmEnv[,,'hosp'])))
g.hosp <- ggplot(data=meanlambds_hosp,aes(x = as.factor(hosp),y=lambda))+
  geom_violin(data=lbds,fill='#50938a')+
  geom_errorbar(aes(ymin=lbd.cil,ymax=lbd.cih),width=.25)+
  geom_linerange(aes(ymin=lbd.q25,ymax=lbd.q75),size=2,color='#7eb3eb')+
  geom_point(shape=23,fill='black')+
  # geom_path()+
  annotate(geom='text',x=0.6,y=0.38,label='n=')+
  annotate(geom='text',x=hospD$hosp+1,y=0.38,label=hospD$n,size=4,color='grey50')+
  scale_x_discrete(labels=paste0(seq(0,100,by=10),'%'))+
  labs(y=expression(paste("Population growth rate (",lambda,")")),
       x="House sparrow occupation")+
  theme_bw(16)+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
g.hosp

ggsave('~/Downloads/hosp.png',bg='transparent',w=10,h=5.5)


g.hosp <- ggplot(data=lbds,aes(x = as.factor(hosp),y=s3))+
  geom_violin(fill='#50938a')+
  # (aes(ymin=lbd.cil,ymax=lbd.cih),width=.25)+
  stat_summary(fun.min = function(x) quantile(x,0.25),
               fun.max = function(x) quantile(x,0.75),
               geom = 'linerange',size=2,color='#7eb3eb')+
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = 'errorbar',width=.25)+
  stat_summary(fun = function(x) median(x),geom = 'point',shape=23,fill='black')+
  # annotate(geom='text',x=0.6,y=0.11,label='n=')+
  # annotate(geom='text',x=hospD$hosp+1,y=0.11,label=hospD$n,size=4,color='grey50')+
  scale_x_discrete(labels=paste0(seq(0,100,by=10)))+
  labs(y='Survival of ASY',
       x="House sparrow occupation (%)")+
  theme_bw(12)+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
g.hosp


g.hosp2 <- ggplot(data=lbds,aes(x = as.factor(hosp),y=r3))+
  geom_violin(fill='#50938a')+
  # (aes(ymin=lbd.cil,ymax=lbd.cih),width=.25)+
  stat_summary(fun.min = function(x) quantile(x,0.25),
               fun.max = function(x) quantile(x,0.75),
               geom = 'linerange',size=2,color='#7eb3eb')+
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = 'errorbar',width=.25)+
  stat_summary(fun = function(x) median(x),geom = 'point',shape=23,fill='black')+
  # annotate(geom='text',x=0.6,y=0,label='n=')+
  # annotate(geom='text',x=hospD$hosp+1,y=0,label=hospD$n,size=4,color='grey50')+
  scale_x_discrete(labels=paste0(seq(0,100,by=10)))+
  labs(y='Reproductive success of ASY',
       x="House sparrow occupation (%)")+
  theme_bw(12)+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

g.hosp3 <- ggplot(data=lbds,aes(x = as.factor(hosp),y=f3))+
  geom_violin(fill='#50938a')+
  # (aes(ymin=lbd.cil,ymax=lbd.cih),width=.25)+
  stat_summary(fun.min = function(x) quantile(x,0.25),
               fun.max = function(x) quantile(x,0.75),
               geom = 'linerange',size=2,color='#7eb3eb')+
  stat_summary(fun.min = function(x) quantile(x,0.025),
               fun.max = function(x) quantile(x,0.975),
               geom = 'errorbar',width=.25)+
  stat_summary(fun = function(x) median(x),geom = 'point',shape=23,fill='black')+
  annotate(geom='text',x=0.6,y=2,label='n=')+
  annotate(geom='text',x=hospD$hosp+1,y=2,label=hospD$n,size=4,color='grey50')+
  scale_x_discrete(labels=paste0(seq(0,100,by=10)))+
  scale_y_continuous(limits = c(2,6))+
  labs(y='Number of Fledgling',
       x="House sparrow occupation (%)")+
  theme_bw(12)+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
library(cowplot)
plot_grid(g.hosp,g.hosp2,g.hosp3,ncol=1)
ggsave('~/Downloads/SurvHosp.png',bg='transparent',w=11,h=5.5)



# plot yr estimates -------------------------------------------------------

allYrMat <- lapply(1:13,function(t){
  g.Pred.LMat(newd = yrMeanEnvs[t,],betaz = betazz,t = t,pim=t )
})

allYrMat_noImmi <- lapply(1:13,function(t){
  g.Pred.LMat(newd = yrMeanEnvs[t,],betaz = betazz,t = t,pim=0 )
})

meanlambds_yr <- allYrMat%>% map_dfr(function(t){
  lbdz <- apply(t$LMatrix, 3, function(it) eigen.analysis(it)$lambda1)
  tibble(lambda=mean(lbdz),
         lbd.cil=quantile(lbdz,0.025),
         lbd.cih=quantile(lbdz,0.975))
}) %>% cbind(.,yr=c(1:13)+2005)

ggplot(meanlambds_yr,aes(x=yr,y=lambda,ymin=lbd.cil,ymax=lbd.cih))+
  geom_pointrange()+
  geom_hline(yintercept = 1)+  theme_bw(18)+
  labs(x='Year',y= expression(paste("Predicted population growth rate (",lambda,")")))+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

ggsave('~/Downloads/Local_trswDens.png',bg='transparent',w=10,h=5.5)


meanlambds <- allYrMat[c(1:13)] %>% map_dfr(function(t){
  lbdz <- apply(t$LMatrix, 3, function(it) eigen.analysis(it)$lambda1)
  tibble(lambda=median(lbdz),
         lbd.cil=quantile(lbdz,0.025),
         lbd.cih=quantile(lbdz,0.975))
})%>% cbind(.,yr=c(1:13)+2005)


meanMatz <- lapply(allYrMat[c(1:13)],function(x) {apply(x$LMatrix,1:2,median)})
stlambds=stoch.growth.rate(meanMatz)
stlambds=data.frame(yr=2019,
                    stoch=exp(stlambds$sim),
                    cil=exp(stlambds$sim.CI[1]),
                    cih=exp(stlambds$sim.CI[2]))



pgrowth=data.frame(
  yr=2006:2018,
  dndt=exp(diff(log(nb1+nb2+nb3))),
  # lambdaNoimmi=map_dbl(meanMatzNoImmi,function(it) eigen.analysis(it)$lambda1),
  lambda_t=NA,
  N_tm1=NA,
  N=NA
)
pgrowth <- full_join(pgrowth,meanlambds)
for(t in 1:nrow(pgrowth)){
  pgrowth$N_tm1[t] <- sum(c(nb1[t],nb2[t],nb3[t]))
  # pgrowth$N[t] <- sum(c(nb1[t],nb2[t],nb3[t]) %*%  meanMatz[[t]])
  pgrowth$lambda_t[t] <- pgrowth$N[t]/ pgrowth$N_tm1[t]
}


ggplot(pgrowth[,],aes(x=yr))+
  geom_hline(yintercept = 1)+
  geom_path(aes(y=lambda,color='lambda'))+
  geom_linerange(aes(ymin=lbd.cil,ymax=lbd.cih,color='lambda'))+
  geom_path(aes(y=dndt,color='dndt'))+
  geom_path(aes(y=lambda_t,color='lambda_t'))+
  geom_pointrange(data=stlambds,aes(x=2020,y=stoch,ymin=cil,ymax=cih,color='Stoch.'))+
  scale_y_continuous(breaks = seq(0,2,by=0.2),limits = c(0,2))+
  scale_x_continuous(breaks = seq(2006,2020,by=2))

cor(pgrowth$lambda[-13],pgrowth$dndt[-13])
plot(pgrowth$lambda[-13],pgrowth$dndt[-13]);abline(0,1)

pgrowth$dndt[-13] %>% var()
pgrowth$lambda[-13] %>% var()
# ellasticity -------------------------------------------------------------



elest_VR <- allYrMat %>% map(function(.y){
  map_df(1:dim(.y$LMatrix)[3],function(i){
    .x=.y$LMatrix[,,i]
    out=as.data.frame(matrix(as.numeric(eigen.analysis(.x)$elasticities),nrow=1))
    colnames(out) <- paste(rep(1:3,each=3),rep(1:3,3),sep='->')
    out$lambda <- eigen.analysis(.x)$lambda1
    out$e.s1 <-  out$`1->2`+ out$`1->1`
    out$e.s2 <-  out$`2->3`+ out$`2->1`
    out$e.s3 <-  out$`3->3`+ out$`3->1`
    out$e.r1 <-  out$`1->1`
    out$e.r2 <-   out$`2->1`
    out$e.r3 <- out$`3->1`
    out$e.f1 <-   out$`1->1`
    out$e.f2 <-   out$`2->1`
    out$e.f3 <-   out$`3->1`
    out$e.i1 <- out$`1->2`
    out$e.i2 <- out$`2->3`
    out$e.i3 <- out$`3->3`

    Sensi <- eigen.analysis(.x)$sensitivities
    out$s.s1 <- Sensi[2,1]*.y$VR$i1[i] + Sensi[1,1]*.y$VR$r1[i]*.y$VR$f1[i]*0.5
    out$s.s2 <- Sensi[3,2]*.y$VR$i2[i] + Sensi[1,2]*.y$VR$r2[i]*.y$VR$f2[i]*0.5
    out$s.s3 <- Sensi[3,3]*.y$VR$i3[i] + Sensi[1,3]*.y$VR$r3[i]*.y$VR$f3[i]*0.5
    out$s.r1 <- Sensi[1,1]* .y$VR$s1[i] * .y$VR$f1[i] * 0.5
    out$s.r2 <- Sensi[1,2]* .y$VR$s2[i] * .y$VR$f2[i] * 0.5
    out$s.r3 <- Sensi[1,3]* .y$VR$s3[i] * .y$VR$f3[i] * 0.5
    out$s.f1 <- Sensi[1,1]* .y$VR$s1[i] * .y$VR$r1[i] * 0.5
    out$s.f2 <- Sensi[1,2]* .y$VR$s2[i] * .y$VR$r2[i] * 0.5
    out$s.f3 <- Sensi[1,3]* .y$VR$s3[i] * .y$VR$r3[i] * 0.5
    out$s.i1 <- Sensi[2,1]* .y$VR$s1[i]
    out$s.i2 <- Sensi[3,2]* .y$VR$s2[i]
    out$s.i3 <- Sensi[3,3]* .y$VR$s3[i]

    return(out)
  })
} )

elest_summ <- map_dfr(1:length(elest_VR), function(i){
  x=elest_VR[[i]]
  tibble(t=i,
         var=colnames(x),
         mean=apply(x,2, mean),
         sd=apply(x,2, sd),
         cil=apply(x,2, quantile, 0.025),
         cih=apply(x,2, quantile, 0.975)
  )
})


g.sens <- elest_summ %>% filter(str_detect(var,'^s')) %>%
  ggplot(.,aes(x=substr(var,3,12),y=mean,ymin=cil, ymax=cih,color=as.factor(t)))+
  geom_point(position = position_dodge(w=0.4),size=0.75)+
  geom_linerange(position = position_dodge(w=0.4),size=0.5)+
  coord_flip()+
  theme_bw(16)+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) + guides(color=F)+
  labs(x='Fitness component', y='sensitivity')



g.el <- elest_summ %>% filter(str_detect(var,'^e')) %>%
  ggplot(.,aes(x=substr(var,3,12),y=mean,ymin=cil, ymax=cih,color=as.factor(t)))+
  geom_point(position = position_dodge(w=0.4),size=0.75)+
  geom_linerange(position = position_dodge(w=0.4),size=0.5)+
  coord_flip()+
  theme_bw(16)+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) + guides(color=F)+
  labs(x='Fitness component', y='elasticity')
ggsave('~/Downloads/Alest.png',bg='transparent',w=11,h=5.5)

library(cowplot)
plot_grid(g.el,g.sens)


elest_summ %>% filter(var=='lambda') %>%
  ggplot(.,aes(x=t+2005,y=mean,ymin=cil, ymax=cih,color=as.factor(t)))+
  geom_point(position = position_dodge(w=0.4),size=0.75)+
  geom_linerange(position = position_dodge(w=0.4),size=0.5)+
  theme_gab()+guides(color=F)+
  labs(x='year', y='lambda')



# varPart -----------------------------------------------------------------
library(parallel)

varHOSP <- mclapply(1:100,function(i1){
    inp <- mydat$r.farmYrEnv[sample(1:40,1),sample(1:13,1),] %>% t() %>% as.tibble() %>%
        rename(prec.cold=`prec:cold`,
               prec.IA=`prec:IA`,
               IA.cold=`IA:cold`)
    ci <- expand_grid(f=1:40,t=1:13)
    t <- sample(1:13,1)
    f <- sample(1:40,1)
    useit <- sample(1:nrow(betazz))
    out1 <- sapply(sample(1:nrow(ci),size=10),function(i2){
        inp=inp %>% mutate(hosp=mydat$r.farmYrEnv[ci$f[i2],ci$t[i2],'hosp'])
        out1 <- g.Pred.LMat(newd =inp ,betaz = betazz[useit,],f=f,t = t ,pim = t)$LMatrix
        apply(out1,3,function(it){lambda(it)})
    })
    return(apply(out1,1,var))
},mc.cores = 2)
varHOSP %>% unlist() %>% hist(30)


# varAI <- mclapply(1:100,function(i1){
#   inp <- mydat$r.farmYrEnv[sample(1:40,1),sample(1:13,1),] %>% t() %>% as.tibble() %>%
#     rename(prec.cold=`prec:cold`,
#            prec.IA=`prec:IA`,
#            IA.cold=`IA:cold`)
#   ci <- expand_grid(f=1:40,t=1:13)
#   t <- sample(1:13,1)
#   useit <- sample(1:nrow(betazz))
#   out1 <- sapply(sample(1:nrow(ci),size=10),function(i2){
#     inp=inp %>% mutate(IA=mydat$r.farmYrEnv[ci$f[i2],ci$t[i2],'IA'])
#     out1 <- g.Pred.LMat(newd =inp ,betaz = betazz[useit,],f=NA,t = NA ,pim = t)$LMatrix
#     apply(out1,3,function(it){lambda(it)})
#   })
#   return(apply(out1,1,var))
# },mc.cores = 2)
# varAI %>% unlist() %>% hist(30)

varFarm <- mclapply(1:100,function(i1){
    inp <- mydat$r.farmYrEnv[sample(1:40,1),sample(1:13,1),] %>% t() %>% as.tibble() %>%
        rename(prec.cold=`prec:cold`,
               prec.IA=`prec:IA`,
               IA.cold=`IA:cold`)
    t <- sample(1:13,1)
    useit <- sample(1:nrow(betazz))
    out1 <- sapply(sample(1:40,size=10),function(i2){
        out1 <- g.Pred.LMat(newd =inp ,betaz = betazz[useit,],f=i2,t = t ,pim = t)$LMatrix
        apply(out1,3,function(it){lambda(it)})
    })
    return(apply(out1,1,var))
},mc.cores = 2)
varFarm %>% unlist() %>% hist(30)


varYr <- mclapply(1:100,function(i1){
    inp <- mydat$r.farmYrEnv[sample(1:40,1),sample(1:13,1),] %>% t() %>% as.tibble() %>%
        rename(prec.cold=`prec:cold`,
               prec.IA=`prec:IA`,
               IA.cold=`IA:cold`)
    t <- sample(1:13,1)
    f <- sample(1:40,1)
    useit <- sample(1:nrow(betazz))
    out1 <- sapply(sample(1:40,size=10),function(i2){
        out1 <- g.Pred.LMat(newd =inp ,betaz = betazz[useit,],f=f,t = i2 ,pim = t)$LMatrix
        apply(out1,3,function(it){lambda(it)})
    })
    return(apply(out1,1,var))
},mc.cores = 2)

varYr %>% unlist() %>% hist(30)


varTib <- tibble(it=1:length(unlist(varHOSP)),
                 # varAI=unlist(varAI),
                 varHOSP=unlist(varHOSP),varFarm=unlist(varFarm),varYr=unlist(varYr))
hist(varTib$varHOSP)
hist(varTib$varFarm)

varTib %>% filter(varYr<0.1) %>%pivot_longer(-it) %>%
    group_by(it) %>% mutate(value=value/sum(value)) %>%
    ggplot(aes(x=value,fill=name))+geom_density(alpha=0.3)
varTib %>% filter(varYr<0.1) %>%pivot_longer(-it) %>%
    group_by(it) %>% mutate(value=value/sum(value)) %>%
    group_by(name) %>%
    summarise(printMeanCI(value))