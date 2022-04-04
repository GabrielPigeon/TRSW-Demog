library(tidyverse)
load('cache/cleanMultiState.Rdata')
l <- ls()
l <- l[-which(l %in% c('myconst','mydat'))]
rm(list = l)
gc()

longDat <- map_dfr(1:13, function(t){
    out <- data.frame(
        t=t,
        i=1:nrow(mydat$obs),
        id=rownames(mydat$obs),
        age=as.numeric(myconst$age[,t]),
        obs=as.numeric(mydat$obs[,t]),
        obs_tp1=as.numeric(mydat$obs[,t+1]),
        nbFledge=as.numeric(mydat$nbFledge[,t]),
        nbFledge_tp1=as.numeric(mydat$nbFledge[,t+1]),
        ismat=((1:nrow(mydat$obs)) %in% myconst$iMat),
        isIn= (t>=myconst$first),
        farm=as.numeric(mydat$farm[,t]),
        farm_tp1=as.numeric(mydat$farm[,t+1])
    )
    env_t <- mydat$x.farmYrEnv[as.numeric(mydat$farm[,t]),t,]
    env_tp1 <- mydat$x.farmYrEnv[as.numeric(mydat$farm[,t+1]),t+1,]
    colnames(env_tp1) <- paste0(colnames(env_tp1),'_tp1')
    out <- cbind(out,env_t,env_tp1) 
    return(out)
})

longDat <- longDat %>%mutate(rs=case_when(
    obs==1  ~  0 ,
    obs==2  ~  1 ,
    TRUE ~ NA_real_ ),
    rs_tp1=case_when(
        obs_tp1==1  ~  0 ,
        obs_tp1==2  ~  1 ,
        TRUE ~ NA_real_ )
    ) %>% 
    mutate(seen=as.numeric(obs<3),
           seen_tp1=as.numeric(obs_tp1<3) )

longDat %>% filter(!is.na(rs)) %>% group_by(rs,age) %>% 
    summarise(returnProp=mean(seen_tp1),
              nextRepro=mean(rs_tp1,na.rm=T),
    ) %>% arrange(age)


library(lme4)
longDat %>% filter(!is.na(rs),age>1) %>% 
    glmer(seen_tp1~as.factor(rs)*as.factor(age)+(1|t),family='binomial',data=.) %>% 
    summary

longDat %>% filter(!is.na(rs),age>1) %>% 
    glmer(rs_tp1~as.factor(age)*as.factor(rs)+(1|t),family='binomial',data=.) %>% 
    summary

tmp <- filter(longDat,!is.na(rs_tp1) &!is.na(farm_tp1))
modrs <- glmer(rs_tp1~as.factor(age)*(
    precip.rearing_tp1 + coldSnap.rearing_tp1+IA_tp1+`prec:cold_tp1`+`prec:IA_tp1`+`IA:cold_tp1`+triple_tp1 
    )+(1|t)+(1|farm_tp1),family='binomial',data=tmp,na.action = 'na.fail')
modrs%>%summary

dredge <- MuMIn::dredge(modrs)

