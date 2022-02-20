source('R/999_MyFunc.R')
library(readxl)
library(lme4)
library(abind)

options(stringsAsFactors=F)

# Load data ----------------------------
# load trsw data
adult_dt <- read.csv2('data/Adultes_2004-2019.csv',na = c('',' ','NA','Na','nA','na','NA\n','\n'),dec=c('.'))
glimpse(adult_dt)
adult_dt <- adult_dt %>%
  mutate_at(vars(jjulien,masse,laile1,laile2,tarse1,tarse2),~ str_replace(.,pattern = ',',replacement = '.')) %>%
  mutate_at(vars(jjulien,masse,laile1,laile2,tarse1,tarse2),as.numeric) %>%
  mutate(age_morpho=ifelse(age_morpho=="AHY" & sexe_gen=="F","ASY",age_morpho))

ois_dt <- read.csv2('data/Oisillons_2004-2019.csv',na = c('',' ','NA','Na','nA','na','NA\n'))
glimpse(ois_dt)
ois_dt <- ois_dt %>%
  mutate_at(vars(jjulien,masse),~ str_replace(.,pattern = ',',replacement = '.')) %>%
  mutate_at(vars(jjulien,masse),as.numeric)

couv_dt <- read.csv2('data/Couvee_2004-2019.csv',na = c('',' ','NA','Na','nA','na','NA\n'))
couv_hibi=couv_dt %>% filter(!is.na(idF1) & codesp==1) %>% mutate()

ois_dt$envol %>% table(useNA = 'a')
ois_dt%>% filter(envol %in% 1) %>% pull(sexe_gen) %>%  table(useNA = 'a')

tmp <- ois_dt %>% filter(!is.na(idois) & envol==1) %>% group_by(idois,idcouvee) %>%
  summarise(sexe_gen=unique(sexe_gen),
            envol=max(envol,na.rm=T)) %>%
  group_by(idcouvee) %>%
  summarise(nbfF=sum(sexe_gen %in% 'F' & envol %in% 1),
            nbfM=sum(sexe_gen %in% 'M' & envol %in% 1),
            nbfI=sum((sexe_gen %in% 'I' | is.na(sexe_gen)) & envol %in% 1),
            nbfALL= nbfF+nbfM+nbfI  )
couv_hibi <- couv_hibi %>% left_join(tmp) %>%
  mutate(nbfF=ifelse(noisenvol==0,0,nbfF))



### error management
couv_hibi[couv_hibi$idcouvee==200120171,"noisenvol"] <- 0
adult_dt[adult_dt$idadult==218161414 &adult_dt$annee==2017,"age_morpho"] <- rep('ASY',3)
adult_dt <- adult_dt[-which(adult_dt$idadult==262181219 &adult_dt$annee==2016),]


couv_hibi[which(couv_hibi$noisenvol!=couv_hibi$nbfALL),]
couv_hibi[couv_hibi$idcouvee==380220082,c('nbfF','nbfALL')] <- 4
couv_hibi[couv_hibi$idcouvee==260820131,c('nbfI','nbfALL')] <- c(0,5)

# hist(round(couv_hibi$noisenvol/2))
# hist(couv_hibi$nbfF)

mean(couv_hibi$nbfF[couv_hibi$nbfF>0],na.rm=T)
var(couv_hibi$nbfF[couv_hibi$nbfF>0],na.rm=T)

sex.ratio <- couv_hibi %>%ungroup %>%  group_by(annee) %>%
  summarise(nbfF=sum(nbfF,na.rm = T),
            nbfM=sum(nbfM,na.rm = T),
            propFem= nbfF/(nbfF+nbfM))
sex.ratio$propFem[-c(1:2,16)] %>% mean
# couv_hibi$noisenvol <- couv_hibi$nbfF   # using female nb of
# couv_hibi$yrfarm <- paste(couv_hibi$annee,couv_hibi$ferme)
# glmer(cbind(nbfF,nbfM)~1 +(1|annee)+(1|ferme) +(1|yrfarm) ,data=couv_hibi[couv_hibi$annee %in% 2006:2018,],family = 'binomial') %>% summary()

# env ---------------------------------------------------------------------
# Define season
#
# migration duration=27(sd=18) @gowEffectsSpringMigration2019
# breeding arrival date = 115 (sd=11) @gowEffectsSpringMigration2019
# autumn migration period (01 Nov–01 Dec) + spring migration period (10 Mar–10 Apr) @bradleyTransGulfMexicoLoop2014


# MY definition  ------.
# springMig.JJ = 91:120  # (1-avr : 30 Avril)
# breeding.JJ = 121:243   # 1  may to 31 aug
# fallMig.JJ = 244:334   # 1 sept- 30 nov
# winter =                # 1 dec : 31 March
# rearing =              # 1 june : 15-july
# earlyBreeding| 121:151       | 1may: 31may    |   5      |

# springMig.M=4
# breeding.M=5:8
# fallMig.M=9:11
# yrround=5:12,1:4


# #  global env  ------------------------------
# nao_seas <-  read_table2("data/norm.nao.monthly.b5001.current.ascii", col_names = c('annee','month','noa')) %>%
#   pivot_wider(names_from = month,values_from=noa) %>%
#   mutate(annee=annee,NAO.spring=`4`,NAO.breeding=(`5`+`6`+`7`+`8`)/4,NAO.fall=(`9`+`10`+`11`)/3) %>%
#   mutate(l1=lead(`1`,1),l2=lead(`2`,1),l3=lead(`3`,1),l4=lead(`4`,1)) %>%
#   mutate(NAO.winterAft=(`12`+l1+l2+l3)/4) %>%
#   mutate(NAO.2nextYr=(`5`+`6`+`7`+`8`+`9`+`10`+`11`+`12`+l1+l2+l3+l4)/12) %>%
#   select(annee,starts_with('NAO')) %>%
#   mutate(NAO.springAft=lead(NAO.spring)) %>%
#   select(annee,NAO.spring,NAO.2nextYr,NAO.breeding,NAO.fall,NAO.winterAft,NAO.springAft)
# 
# BEST_seas <- read_table2("data/BEST.txt",  col_names = c('annee',1:12)) %>%
#   mutate(annee=annee,BEST.spring=`4`,BEST.breeding=(`5`+`6`+`7`+`8`)/4,BEST.fall=(`9`+`10`+`11`)/3) %>%
#   mutate(l1=lead(`1`,1),l2=lead(`2`,1),l3=lead(`3`,1),l4=lead(`4`,1)) %>%
#   mutate(BEST.winterAft=(`12`+l1+l2+l3)/4) %>%
#   mutate(BEST.2nextYr=(`5`+`6`+`7`+`8`+`9`+`10`+`11`+`12`+l1+l2+l3+l4)/12) %>%
#   select(annee,starts_with('BEST')) %>%
#   mutate(BEST.springAft=lead(BEST.spring)) %>%
#   select(annee,BEST.spring,BEST.2nextYr,BEST.breeding,BEST.fall,BEST.winterAft,BEST.springAft)
# 
# ONI_seas <- read_table2("data/oni.data", skip = 1, col_names = c('annee',1:12)) %>%
#   mutate(annee=annee,ONI.spring=`4`,ONI.breeding=(`5`+`6`+`7`+`8`)/4,ONI.fall=(`9`+`10`+`11`)/3) %>%
#   mutate(l1=lead(`1`,1),l2=lead(`2`,1),l3=lead(`3`,1),l4=lead(`4`,1)) %>%
#   mutate(ONI.winterAft=(`12`+l1+l2+l3)/4) %>%
#   mutate(ONI.2nextYr=(`5`+`6`+`7`+`8`+`9`+`10`+`11`+`12`+l1+l2+l3+l4)/12) %>%
#   select(annee,starts_with('ONI')) %>%
#   mutate(ONI.springAft=lead(ONI.spring)) %>%
#   select(annee,ONI.spring,ONI.2nextYr,ONI.breeding,ONI.fall,ONI.winterAft,ONI.springAft)
# 
# soi_seas <- read_table2("data/soi.long.data",  col_names = c('annee',1:12), skip = 1)%>%
#   mutate_all(function(x)ifelse(x<(-50),NA,x)) %>%
#   mutate(annee=annee,SOI.spring=`4`,SOI.breeding=(`5`+`6`+`7`+`8`)/4,SOI.fall=(`9`+`10`+`11`)/3) %>%
#   mutate(l1=lead(`1`,1),l2=lead(`2`,1),l3=lead(`3`,1),l4=lead(`4`,1)) %>%
#   mutate(SOI.winterAft=(`12`+l1+l2+l3)/4) %>%
#   mutate(SOI.2nextYr=(`5`+`6`+`7`+`8`+`9`+`10`+`11`+`12`+l1+l2+l3+l4)/12) %>%
#   select(annee,starts_with('SOI')) %>%
#   mutate(SOI.springAft=lead(SOI.spring)) %>%
#   select(annee,SOI.spring,SOI.2nextYr,SOI.breeding,SOI.fall,SOI.winterAft,SOI.springAft)
# # soi_seas[soi_seas$annee==2019,'SOI.spring'] <- 0
# 
# cor(soi_seas[139:154,],use = 'p')
# 
# 
# winterEnv <- read_csv('data/winterEnv.csv')
# winds <- read_csv('data/winds.csv')
# winds <- winds %>% mutate(u_springAft=lead(u_spring),
#                           v_springAft=lead(v_spring)) %>%
#   mutate_at(vars(u_fall:v_springAft),scale) %>%
#   select(-u_spring,-v_spring)
# 
# # write_csv(nao_seas,'data/naoSeason.csv')
# # write_csv(soi_seas,'data/soiSeason.csv')
# # write_csv(BEST_seas,'data/bestSeason.csv')
# 
# 
# TenvDT <- nao_seas %>% full_join(BEST_seas) %>%full_join(ONI_seas) %>% full_join(soi_seas) %>%
#   full_join(winds,by=c('annee'='yr')) %>%
#   full_join(winterEnv[,-1],by=c('annee'='yrstart')) %>%
#   arrange(annee) %>% filter(annee %in% 2004:2019)
# 
# write_csv(TenvDT,'data/allYrEnv.csv')
# 
# # Local Environemnt : famr*yr ---------------------------------------------
# 
# 
# Temperature_local <- read_excel("data/Temperature_2004-2019.xlsx", na = "NA", n_max = 103831) %>%
#   mutate_at(vars(annee,jjulien:tempmin),.funs = ~ as.numeric(.)) %>%
#   mutate(breeding=jjulien %in% 121:243 & inout %in% c('out','OUT'),
#          rearing=jjulien %in% 152:196 & inout %in% c('out','OUT'),
#          earlyBreeding=jjulien %in% 121:151 & inout %in% c('out','OUT')) %>%
#   mutate(ferme=ifelse(nchar(ferme)==1,paste0('0',ferme),ferme)) %>%
#   group_by(annee,ferme) %>%
#   summarise(localMeanTemp.breeding=mean(ifelse(breeding,moy24,NA),na.rm = T),
#             localMaxTemp.breeding=mean(ifelse(breeding,tempmax,NA),na.rm = T),
#             localMeanTemp.rearing=mean(ifelse(rearing,moy24,NA),na.rm = T),
#             localMaxTemp.rearing=mean(ifelse(rearing,tempmax,NA),na.rm = T),
#             localMeanTemp.earlyBreeding=mean(ifelse(earlyBreeding,tempmin,NA),na.rm = T),
#             localMaxTemp.earlyBreeding=max(ifelse(earlyBreeding,tempmin,NA),na.rm = T),
#             coldSnap.rearing=sum((tempmax<18.8 & rearing) ,na.rm = T),
#             coldSnap.earlyBreeding=sum((tempmax<18.8 & earlyBreeding) ,na.rm = T)
#   ) %>%
#   mutate_if(is.numeric,function(x)ifelse(is.nan(x)|is.infinite(x),NA,x))
# 
# 
# 
# Temperature_local=expand_grid(ferme=substr(101:140,2,3),annee=2004:2019) %>% left_join(Temperature_local)
# Temperature_local %>% ggplot(aes(x=annee,y=localMeanTemp.breeding,color=ferme,shape=))+geom_path()
# mm=lmer(localMeanTemp.breeding~as.factor(annee)+(1|ferme),data=Temperature_local)
# mmx=lmer(localMaxTemp.breeding~as.factor(annee)+(1|ferme),data=Temperature_local)
# 
# rmm=lmer(localMeanTemp.rearing~as.factor(annee)+(1|ferme),data=Temperature_local)
# rmmx=lmer(localMaxTemp.rearing~as.factor(annee)+(1|ferme),data=Temperature_local)
# 
# mmin=lmer(localMeanTemp.earlyBreeding~as.factor(annee)+(1|ferme),data=Temperature_local)
# rmmin=lmer(localMaxTemp.earlyBreeding~as.factor(annee)+(1|ferme),data=Temperature_local)
# 
# ndt=expand_grid(ferme=substr(101:140,2,3),annee=2004:2019)
# ndt$p.mean=predict(mm,ndt,re.form=~(1|ferme))
# ndt$p.mx=predict(mmx,ndt,re.form=~(1|ferme))
# Temperature_local$localMeanTemp.breeding <- ifelse(is.na(Temperature_local$localMeanTemp.breeding),ndt$p.mean,Temperature_local$localMeanTemp.breeding)
# Temperature_local$localMaxTemp.breeding <- ifelse(is.na(Temperature_local$localMaxTemp.breeding),ndt$p.mx,Temperature_local$localMaxTemp.breeding)
# ndt$p.mean=predict(rmm,ndt,re.form=~(1|ferme))
# ndt$p.mx=predict(rmmx,ndt,re.form=~(1|ferme))
# Temperature_local$localMeanTemp.rearing <- ifelse(is.na(Temperature_local$localMeanTemp.rearing),ndt$p.mean,Temperature_local$localMeanTemp.rearing)
# Temperature_local$localMaxTemp.rearing <- ifelse(is.na(Temperature_local$localMaxTemp.rearing),ndt$p.mx,Temperature_local$localMaxTemp.rearing)
# ndt$p.min=predict(mmin,ndt,re.form=~(1|ferme))
# ndt$p.min=predict(rmmin,ndt,re.form=~(1|ferme))
# Temperature_local$localMeanTemp.earlyBreeding <- ifelse(is.na(Temperature_local$localMeanTemp.earlyBreeding),ndt$p.min,Temperature_local$localMeanTemp.earlyBreeding)
# Temperature_local$localMaxTemp.earlyBreeding <- ifelse(is.na(Temperature_local$localMaxTemp.earlyBreeding),ndt$p.min,Temperature_local$localMaxTemp.earlyBreeding)
# 
# Precip_local <- read_excel("data/Precipitation_2004-2019.xlsx", na = "NA", n_max = 33423) %>%
#   mutate_at(vars(annee,jjulien:intervalle),.funs = ~ as.numeric(.)) %>%
#   mutate(breeding=jjulien %in% 121:243,
#          rearing=jjulien %in% 152:196 ,
#          earlyBreeding=jjulien %in% 121:151 ) %>%
#   mutate(ferme=ifelse(nchar(ferme)==1,paste0('0',ferme),ferme)) %>%
#   group_by(annee,ferme) %>%
#   summarise(localMeanPrec.breeding=mean(ifelse(breeding,precipitations,NA),na.rm=T) ,
#             localMeanPrec.rearing=mean(ifelse(rearing,precipitations,NA),na.rm=T) ,
#             localMeanPrec.earlyBreeding=mean(ifelse(earlyBreeding,precipitations,NA),na.rm=T))
# Precip_local=expand_grid(ferme=substr(101:140,2,3),annee=2004:2019) %>% left_join(Precip_local)
# 
# qplot(data=Precip_local,x=annee,y=ferme,fill=localMeanPrec.earlyBreeding,geom='tile')
# 
# yrFarmEnv=array(NA,dim=c(40,16,14),
#                 dimnames = list(substr(101:140,2,3),2004:2019,c(
#                   'meanTemp.rearing','meanTemp.breeding','meanTemp.earlyBreeding',
#                   'maxTemp.rearing','maxTemp.breeding','maxTemp.earlyBreeding',
#                   'coldSnap.rearing','coldSnap.earlyBreeding',
#                   'precip.breeding','precip.rearing','precip.earlyBreeding',
#                   'hosp','water','IA')))
# 
# 
# yrFarmEnv[1:40,1:16,'meanTemp.rearing'] <- Temperature_local %>%  select(annee, ferme, localMeanTemp.breeding) %>%
#   pivot_wider(names_from = annee,values_from=localMeanTemp.breeding)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'meanTemp.breeding'] <- Temperature_local %>%select(annee, ferme, localMeanTemp.rearing) %>%
#   pivot_wider(names_from = annee,values_from=localMeanTemp.rearing)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'meanTemp.earlyBreeding'] <- Temperature_local %>%select(annee, ferme, localMeanTemp.earlyBreeding) %>%
#   pivot_wider(names_from = annee,values_from=localMeanTemp.earlyBreeding)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# 
# yrFarmEnv[1:40,1:16,'maxTemp.rearing'] <- Temperature_local %>%select(annee, ferme, localMaxTemp.rearing) %>%
#   pivot_wider(names_from = annee,values_from=localMaxTemp.rearing)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'maxTemp.breeding'] <- Temperature_local %>%select(annee, ferme, localMaxTemp.breeding) %>%
#   pivot_wider(names_from = annee,values_from=localMaxTemp.breeding)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'maxTemp.earlyBreeding'] <- Temperature_local %>%select(annee, ferme, localMaxTemp.earlyBreeding) %>%
#   pivot_wider(names_from = annee,values_from=localMaxTemp.earlyBreeding)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# 
# yrFarmEnv[1:40,1:16,'coldSnap.rearing'] <- Temperature_local %>% select(annee, ferme, coldSnap.rearing) %>%
#   pivot_wider(names_from = annee,values_from=coldSnap.rearing)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'coldSnap.earlyBreeding'] <- Temperature_local %>% select(annee, ferme, coldSnap.earlyBreeding) %>%
#   pivot_wider(names_from = annee,values_from=coldSnap.earlyBreeding)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix() %>% round(2)
# 
# yrFarmEnv[1:40,1:16,'precip.breeding'] <- Precip_local %>% select(annee, ferme, localMeanPrec.breeding) %>%
#   pivot_wider(names_from = annee,values_from=localMeanPrec.breeding) %>%
#   arrange(ferme) %>% select(`2004`:`2019`)%>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'precip.rearing'] <- Precip_local %>% select(annee, ferme, localMeanPrec.rearing) %>%
#   pivot_wider(names_from = annee,values_from=localMeanPrec.rearing) %>%
#   arrange(ferme) %>% select(`2004`:`2019`)%>% as.matrix() %>% round(2)
# yrFarmEnv[1:40,1:16,'precip.earlyBreeding'] <- Precip_local %>% select(annee, ferme, localMeanPrec.earlyBreeding) %>%
#   pivot_wider(names_from = annee,values_from=localMeanPrec.earlyBreeding) %>%
#   arrange(ferme) %>% select(`2004`:`2019`)%>% as.matrix() %>% round(2)
# 
# yrFarmEnv[,,'coldSnap.rearing']=apply(yrFarmEnv[,,'coldSnap.rearing'], 2, function(x) ifelse(is.na(x), round(mean(x,na.rm=T)),x))
# yrFarmEnv[,,'coldSnap.earlyBreeding']=apply(yrFarmEnv[,,'coldSnap.earlyBreeding'], 2, function(x) ifelse(is.na(x), round(mean(x,na.rm=T)),x))
# yrFarmEnv[,,'precip.earlyBreeding']=apply(yrFarmEnv[,,'precip.earlyBreeding'], 2, function(x) ifelse(is.na(x), round(mean(x,na.rm=T)),x))
# 
# 
# 
# # get hosp
# hosp_occup <- couv_dt %>% filter(codesp==3,nnich==1) %>%
#   group_by(ferme,annee) %>%
#   summarise(nb.hosp=n()) %>%
#   mutate(ferme=substr(100+ferme,2,3))
# hosp_occup=expand_grid(ferme=substr(101:140,2,3),annee=2004:2019) %>% left_join(hosp_occup) %>%
#   mutate(nb.hosp=ifelse(is.na(nb.hosp),0,nb.hosp))
# 
# yrFarmEnv[1:40,1:16,'hosp'] <- hosp_occup %>% pivot_wider(names_from = annee,values_from=nb.hosp)%>% arrange(ferme) %>% select(`2004`:`2019`) %>% as.matrix()
# 
# ggplot(hosp_occup,aes(x=annee,y=nb.hosp))+geom_point(position = position_jitter(w=0.5,h=0.5))+geom_smooth()
# 
# 
# # get  water
# library(readxl)
# caract <- read_excel("data/caracterisation1a20km.xlsx")
# caract <- caract %>% filter(rayon==5000) %>% select(id, ferme, nichoir,eau,cultannuelle) %>%
#   group_by(ferme) %>%
#   summarise(eau=mean(eau),IA=mean(cultannuelle))
# 
# yrFarmEnv[1:40,1:16,'water'] <- caract %>% arrange(ferme) %>% pull(eau) %>% matrix(.,nrow = 40,ncol=16,byrow = F) %>% round(2)
# # yrFarmEnv[1:40,1:16,'IA'] <- caract %>% arrange(ferme) %>% pull(IA) %>% matrix(.,nrow = 40,ncol=16,byrow = F) %>% round(2)
# 
# 
# scores <- read_csv('data/Site.Scores.2006_2019.csv') %>% arrange(Farm)
# scores$Farm %>% table()
# for(t in 2006:2019){ yrFarmEnv[,t-2003,'IA'] <- scores$Comp.1[scores$Year==t] }
# for(f in 1:40){yrFarmEnv[f,1:2,'IA'] <- mean(yrFarmEnv[f,3:16,'IA'])}
# 
# 
# write_csv(TenvDT,'data/allYrEnv.csv')
# write_rds(yrFarmEnv,'data/yrFarmEnv.rds')
# 
# Temperature_local %>% full_join(Precip_local) %>% write_csv('data/yrFarmEnvLong.csv')
# 


# make Frames -------------------------------------------------------------
# + get IDs   -------------------
id_ad <- adult_dt[,c("idadult",'prefixe','annee','age_morpho','age_exact',"sexe_gen","sexe_morpho")] %>%
  mutate(sex=ifelse(is.na(sexe_gen),sexe_morpho,sexe_gen),
         age_marked=age_morpho) %>%
  # filter(sex=='F') %>%
  filter(!is.na(prefixe)) %>%
  group_by(idadult) %>% arrange(annee) %>% slice_head(n=1) %>% rename(idMark=idadult)


# do any birds change sex within adult file
l=adult_dt[,c("idadult",'prefixe','annee','age_morpho','age_exact',"sexe_gen","sexe_morpho")] %>%
  mutate(sex=ifelse(is.na(sexe_gen),sexe_morpho,sexe_gen),
         age_marked=age_morpho) %>%
  # filter(sex=='F') %>%
  filter(!is.na(prefixe)) %>%
  group_by(idadult,sex) %>%select(idadult,sex) %>%  unique()
l= l[l$idadult %in% l$idadult[duplicated(l$idadult)],] %>% arrange(idadult) %>% na.omit()
l
id_ad$sex[id_ad$idMark %in% c(188197030,188197122,188197154,188197268,188198556)] <- "F"
id_ad$sex[id_ad$idMark==262181507] <- "M"


duplicated(id_ad$idMark) %>% sum() # no dupplicate id

id_ois <- ois_dt %>%
  filter(!is.na(prefixe) & envol==1) %>% select(idois,annee,sexe_gen) %>%
  unique() %>% mutate(age_morpho='ois',age_marked='ois') %>% rename(idMark=idois,sex=sexe_gen)

id_ad[!is.na(id_ad$age_exact),] %>% nrow()

# are there any sex change between ois and adult ; not to bad since keeping gensex of ois by default
merge(id_ois[,c('idMark','sex')], id_ad[,c('idMark','sex')],all=T,by='idMark') %>%  filter(sex.x!=sex.y)



ID_marked <- bind_rows(cbind(id_ois[id_ois$idMark %in% id_ad$idMark,c('idMark','annee','age_marked','sex')],type='b'),
                       cbind(id_ois[!id_ois$idMark %in% id_ad$idMark,c('idMark','annee','age_marked','sex')],type='o'),
                       cbind(id_ad[!id_ad$idMark %in% id_ois$idMark,c('idMark','annee','age_marked','sex')],type='a')
)

ID_marked$idMark %>% unique() %>% length
ID_marked$idMark %>% duplicated() %>% sum()
ID_marked$sex %>% table
ID_marked$age_marked %>% table
is.na(ID_marked$sex) %>% sum()
is.na(ID_marked$age_marked) %>% sum()

ID_marked <- ID_marked %>% filter(sex=='F')
ID_marked$type %>% table

immigrants <- ID_marked %>%filter(type=='a') %>%  group_by(annee,age_marked) %>%
  summarise(Ni=length(unique(idMark))) %>%
  pivot_wider(names_from=age_marked,values_from=Ni,values_fill=0)

ID_marked <- ID_marked %>% filter(sex=='F') %>% filter(annee<2019) %>% arrange(annee)


counts <- id_ad %>% filter(sex=='F' &!is.na(age_morpho)) %>% group_by(annee,age_morpho) %>% summarise(counts=length(unique(idMark))) %>%
  pivot_wider(names_from=age_morpho,values_from=counts,values_fill=0)

tmp <- id_ois %>% filter(sex=='F') %>% group_by(annee) %>% summarise(counts=length(unique(idMark)))
counts <- counts %>% left_join(tmp)

# fill frames  -----------------------------------

years <- adult_dt$annee %>% unique()

state <- obs <- seen <-age <-  RS <-farm <-  nbFledge <- matrix(NA,nrow = nrow(ID_marked),ncol=length(years),
                                         dimnames = list(ID_marked$idMark,years))
yearMarked=match(ID_marked$annee,colnames(obs))
first <- rep(NA,nrow(state))
seen[,] <- 0

obs[,] <- 3
for(i in 1:nrow(obs)){
  couv_id <- couv_hibi[couv_hibi$idF1==ID_marked$idMark[i],]
  adult_id <- adult_dt[adult_dt$idadult==ID_marked$idMark[i],]
  state[i,1:(yearMarked[i]-1)] <- 0
  for(t in yearMarked[i]:ncol(obs)){
    if(t==yearMarked[i]){# if marking year, was obviously seen
      seen[i,t] <- 1
      if(ID_marked$age_marked[i]=="ois"){ #if ois,rs not-important
        curid=unique(na.omit(ois_dt$id[ois_dt$idois==ID_marked$idMark[i]]))[1]
        curferme=as.numeric(substr(curid,1,2))
        age[i,t] <- 1
        if(t<ncol(obs)) age[i,(t+1):ncol(obs)] <- 2
        if((t+1)<ncol(obs)) age[i,(t+2):ncol(obs)] <- 3
        state[i,t] <- 1
        obs[i,t] <- 1
      }else{
        curid=table(na.omit(adult_id$id[adult_id$annee==(years[t])]))
        curid=names(curid[which.max(curid)])
        curferme=as.numeric(substr(curid,1,2))
        rs <- couv_id$noisenvol[couv_id$annee==(years[t])]
        nbFledge[i,t] <- ifelse(length(rs)==0,0,sum(rs,na.rm=T)) # if not in couv but seen, RS=0
        RS[i,t] <- nbFledge[i,t]>0
        if(ID_marked$age_marked[i]=="SY"){
          age[i,t] <- 2
          if(t<ncol(obs)) age[i,(t+1):ncol(obs)] <- 3
          if(RS[i,t]==T){
            state[i,t] <- 2 ; obs[i,t] <- 2
          } else {
            state[i,t] <- 1 ; obs[i,t] <- 1
          }
        }
        if(ID_marked$age_marked[i]=="ASY"){
          age[i,t:ncol(obs)] <- 3
          if(RS[i,t]==T){
            state[i,t] <- 2 ; obs[i,t] <- 2
          } else {
            state[i,t] <- 1 ; obs[i,t] <- 1
          }
        }
      }
    if(curferme %in% 1:40){ farm[i,t] <- as.numeric(as.character(curferme))}
    } else {
    if(sum(adult_id$annee==(years[t]))>0  ){ # was captured at least ounce
      curid=table(na.omit(adult_id$id[adult_id$annee==(years[t])]))
      curid=names(curid[which.max(curid)])
      curferme=as.numeric(substr(curid,1,2))
      if(curferme %in% 1:40){ farm[i,t] <- as.numeric(as.character(curferme))}
      seen[i,t] <- 1
      rs <- couv_id$noisenvol[couv_id$annee==(years[t])]
      nbFledge[i,t] <- ifelse(length(rs)==0,0,sum(rs,na.rm=T))
      RS[i,t] <- nbFledge[i,t]>0
      if("SY" %in% adult_id$age_morpho[adult_id$annee==(years[t])]){
        # Age[i,t] <- 2
        if(RS[i,t]==T){
          state[i,t] <- 2 ; obs[i,t] <- 2
        } else {
          state[i,t] <- 1 ; obs[i,t] <- 1
        }
      }
      if("ASY" %in% adult_id$age_morpho[adult_id$annee==(years[t])]){
        # Age[i,t] <- 3
        if(RS[i,t]==T){
          state[i,t] <- 2 ; obs[i,t] <- 2
        } else {
          state[i,t] <- 1 ; obs[i,t] <- 1
        }
      }


    }

  }

}

}

age <- ifelse(is.na(age),0,age)



# save output to list for use in nimble  ------------------------------------
mydat <- list()
mydat$obs <- as.matrix(obs)
mydat$state <- as.matrix(state)

mydat$nbFledge=nbFledge
mydat$Dummy=rep(1,500)
mydat$imi=as.matrix(immigrants[,c(3,2)])

myconst <- list()
myconst$yearMarked <- yearMarked
myconst$nb.t=ncol(state)
myconst$nb.id=nrow(state)
myconst$age <- as.matrix(age)
myconst$DF=2
myconst$zero=rep(0,myconst$DF)
nb.id=nrow(state)


# reemove 2004:2005 -------------------------------------------------------
mydat_all <- mydat
myconst_all <- myconst

t=years>2005

mydat$obs <- mydat_all$obs[,t]
lid <- rowSums(mydat$obs<3 ,na.rm = T)>0 # & myconst_all$yearMarked>2
lid2 <- (rowSums(!is.na(farm))==0 ) # also remove id with no farm info
# lid3 <- (rownames(mydat$state) %in% c("188197355", "188197820", "188197145", "188197957"))
lid <- lid & !lid2
mydat$obs <- mydat_all$obs[,t][lid,]
mydat$state <- mydat_all$state[,t][lid,]
mydat$nbFledge <- mydat_all$nbFledge[,t][lid,]
mydat$imi <- mydat_all$imi[t,]
mydat$farm <- farm[,t][lid,]


myconst$first <- ifelse(myconst_all$yearMarked<=2,1,myconst_all$yearMarked-2)[lid] # adjust capture index to account for removed 2004,2005
myconst$yearMarked <-(myconst_all$yearMarked-2)[lid]
myconst$nb.t <- ncol(mydat$state)
myconst$nb.id <- nrow(mydat$state)
myconst$age <- myconst$age[,t][lid,]

myconst$firstMat <- apply(myconst$age*(mydat$obs<3),1,function(i) min(which(i>=2))) # get first observed repro to start estimating nb fledged to save time by skipping nestleing
myconst$iMat <- which(myconst$firstMat<=myconst$nb.t )
myconst$firstMat <- myconst$firstMat[myconst$iMat]
myconst$nb.mat <- length(myconst$iMat)
myconst$zero <- rep(0,20)
myconst$W <- diag(20)

# give starting farm for rare birds not observed at first time step (caugh in 2004 or 2005 and not seen in 2006)
lf=sapply(1:nrow(mydat$obs),function(i){mydat$farm[i,myconst$first[i]]}) %>%
    is.na() %>% which
for(i in lf){
  mydat$farm[i,myconst$first[i]] <- na.omit(farm[which(lid==T)[i],])[1]
}

# gives them an arbitratry 1 state
lf=sapply(1:nrow(mydat$obs),function(i){mydat$state[i,myconst$first[i]]}) %>%
  is.na() %>% which
for(i in lf){
  mydat$state[i,myconst$first[i]] <- 1
}


# add dist and envVar to dat list

farmDist <- read_rds('data/FarmDist.rds')
myconst$nb.farm=40
myconst$Dsq=farmDist^2
myconst$D=farmDist


allYrEnv <- read_csv('data/allYrEnv.csv')
yrFarmEnv <- read_rds('data/yrFarmEnv.rds')



aditinalV <- array(NA,c(40,16,dim(allYrEnv)[2]-1),
                   dimnames = list(dimnames(yrFarmEnv)[[1]],dimnames(yrFarmEnv)[[2]],
                                   colnames(allYrEnv)[-1]) )
for(v in 1:(ncol(allYrEnv)-1)){
  for(t in 1:16){     aditinalV[,t,v] <- as.numeric(allYrEnv[t,v+1])  }
}
yrFarmEnv <- abind(yrFarmEnv,aditinalV)
# yrFarmEnv[,,"precip.earlyBreeding"] %>% hist


envScalingPar <- tibble(var=dimnames(yrFarmEnv)[[3]],
                        mean=NA,
                        sd=NA)

for(v in 1:dim(yrFarmEnv)[3]){
  m=mean(yrFarmEnv[,,v],na.rm=T)
  s=sd(yrFarmEnv[,,v],na.rm=T)
  yrFarmEnv[,,v] <- (yrFarmEnv[,,v]-m)/s
  envScalingPar$mean[v] <- m
  envScalingPar$sd[v] <- s
}


x.farmYrEnv <- array(NA,dim = c(40,14,9))
x.farmYrEnv[,,1:5] <-  round(yrFarmEnv[,3:16,c("meanTemp.rearing",'precip.rearing', 'coldSnap.rearing',"hosp" ,'IA')],3)
x.farmYrEnv[,,6] <-  yrFarmEnv[,3:16,c('precip.rearing')]* yrFarmEnv[,3:16,c('coldSnap.rearing')]
x.farmYrEnv[,,7] <- yrFarmEnv[,3:16,c('precip.rearing')]* yrFarmEnv[,3:16,c('IA')]
x.farmYrEnv[,,8] <- yrFarmEnv[,3:16,c('IA')]* yrFarmEnv[,3:16,c('coldSnap.rearing')]
x.farmYrEnv[,,9] <- yrFarmEnv[,3:16,c('precip.rearing')]* yrFarmEnv[,3:16,c('coldSnap.rearing')]*yrFarmEnv[,3:16,c('IA')]
dimnames(x.farmYrEnv)[[3]] <- c("meanTemp.rearing",'precip.rearing', 'coldSnap.rearing',"hosp" ,'IA','prec:cold','prec:IA','IA:cold','triple')
mydat$x.farmYrEnv=round(x.farmYrEnv,3) ;

sum(is.na(x.farmYrEnv))
x.farmYrEnv %>% apply(3, as.numeric) %>% cor(use = 'p') %>% round(2)


# miniData ----------------------------------------------------------------
Reduction=0.25

l=sort(sample(1:nrow(mydat$state),size=round(nrow(mydat$state)*Reduction)))# ; hist(l)

miniDat <- mydat
miniDat$obs <-mydat$obs[l,]
miniDat$state <- mydat$state[l,]

miniDat$RS=mydat$RS[l,]
miniDat$nbFledge=mydat$nbFledge[l,]
miniDat$Dummy=rep(1,10000)
miniDat$imi=round(as.matrix(mydat$imi)*Reduction)
miniDat$farm=mydat$farm[l,]


miniConst <- myconst
miniConst$first=myconst$first[l]
miniConst$yearMarked=miniConst$yearMarked[l]
miniConst$nb.t=ncol(miniDat$state)
miniConst$nb.id=nrow(miniDat$state)
miniConst$age=as.matrix(myconst$age)[l,]
miniConst$nAge=max(miniConst$age,na.rm=T)


miniConst$firstMat <- apply(miniConst$age*(miniDat$obs<3),1,function(i) min(which(i>=2)))
miniConst$iMat <- which(miniConst$firstMat<=miniConst$nb.t )
miniConst$firstMat <- miniConst$firstMat[miniConst$iMat]

miniConst$nb.mat <- length(miniConst$iMat)
miniConst$zero <- rep(0,20)
miniConst$W <- diag(20)

colSums(miniDat$obs)
apply(miniConst$age,2,table)


##micro ----


which(rowSums(!is.na(mydat$nbFledge))==5) %>% sort()
l <- c(12,730,1322,2874,3360,4176,5321, # long life
       # 785   ,   4846    ,  1943  ,
       1885,4730,5671,4169,1236,6073,5351,3780,3900,5659 ,
       3329,3282,3807,4791,1884,2969,5613,2889,6140,4840 ,
       3797,4751,5270,5705,2958,4206,5370,6629,5371,1277 ,
       5655,4206,5608,2467,6045,2446,6629,2982,3382,3280 
       )  %>% unique()
l <- sort(c(l,sample(1:myconst$nb.id,100-length(l))))
microDat <- mydat
microDat$obs <- mydat$obs[l,]
microDat$state <- mydat$state[l,]
microDat$nbFledge <- mydat$nbFledge[l,]
microDat$farm <- mydat$farm[l,]
microConst <- myconst
microConst$first <- myconst$first[l]
microConst$age <- myconst$age[l,]
microConst$nb.id <- length(l)
microConst$first <- myconst$first[l]

# f <- microDat$farm %>% as.numeric() %>% unique() %>% na.omit()
# microConst$D <- microConst$D[f,f]
# microConst$nb.farm <- length(f)
# microDat$farm <- matrix(as.numeric(as.factor(microDat$farm)),nrow = nrow(microDat$farm))

microConst$firstMat <- apply(microConst$age*(microDat$obs<3),1,function(i) min(which(i>=2)))
microConst$iMat <- which(microConst$firstMat<=microConst$nb.t )
microConst$firstMat <- microConst$firstMat[microConst$iMat]
microConst$nb.mat <- length(microConst$iMat)



ArchiveSave(adult_dt,couv_hibi,ois_dt,
            yrFarmEnv,sex.ratio,
            ID_marked,years,envScalingPar,
            mydat,myconst,
            mydat_all,myconst_all,
            miniDat,miniConst,
            microDat,microConst,
            file='cache/cleanMultiState.Rdata')



# checkes -----------------------------------------------------------------

# any rejuvination?
# apply(age,1, function(i){
#   d=diff(na.omit(i))
#   sum(d<0)
# }) %>% table
#
# RS
# table(first)
