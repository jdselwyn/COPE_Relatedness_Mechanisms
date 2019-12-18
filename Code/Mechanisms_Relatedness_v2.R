set.seed(12345)

#### Libraries ####
library(tidyverse)
library(stringr)
library(magrittr)
library(lubridate)
library(adegenet)
library(pegas)
library(HWxtest)
library(related)
library(purrr)
library(broom)
library(modelr)
library(igraph)
library(tidygraph)
library(ggraph)
library(furrr)
library(Cairo)
library(cowplot)
library(fs)
library(INLA)

#### Functions ####
absmax <- function(x) { x[which.max( abs(x) )]}

#### Microsat Data ####
microsats<- read_csv('../Data/Microsat.csv')  %>%
  
  #Remove Cloud 4 since it only has 3 individuals
  filter(Cloud!="COPE-M-04")

#### Shoal Location Data ####
shoal_locations<-read_csv('../Data/COPE clouds.csv', 
                          col_names = c('date','site','bag','Cloud','node.cloud','dist.node','angle.node','dist.cn',
                                        'angle.cn','near.node','number'),
                          col_types = list('c','c','i','c','c','n','n','n','n','c','n'),skip=1) %>%
  mutate(date=dmy(date)) %>%
  
  #Convert angles and distances to XY coordinates
  #X
  group_by(near.node) %>%
  mutate(X=if_else(node.cloud!='Cloud',
                   if_else(node.cloud!='Central Node',round(sin(angle.cn*(pi/180))*(dist.cn*100)),0),
                   round(sin(angle.node*(pi/180))*dist.node))) %>%
  mutate(tmp=1:n(),tmp=ifelse(tmp>1,0,tmp),tmp=tmp*X,tmp=absmax(tmp)) %>%
  mutate(X=if_else((near.node!='Central Node' & node.cloud=='Cloud'),X+tmp,X)) %>%
  
  #Y
  mutate(Y=if_else(node.cloud!='Cloud',
                   if_else(node.cloud!='Central Node',round(cos(angle.cn*(pi/180))*(dist.cn*100)),0),
                   round(cos(angle.node*(pi/180))*dist.node))) %>%
  group_by(near.node) %>%
  mutate(tmp=1:n(),tmp=ifelse(tmp>1,0,tmp),tmp=tmp*Y,tmp=absmax(tmp)) %>%
  mutate(Y=if_else((near.node!='Central Node' & node.cloud=='Cloud'),Y+tmp,Y)) %>%
  select(-tmp) %>%
  mutate(Cloud=recode(Cloud,`COPE-M-01`='A',`COPE-M-03`='B',`COPE-M-05`='C',`COPE-M-08`='D',`COPE-M-09`='E',`COPE-M-10`='F',
                      `COPE-M-11`='G',`COPE-M-12`='H',`COPE-M-13`='I',`COPE-M-14`='J',`COPE-M-15`='K',`COPE-M-18`='L',`COPE-M-20`='M')) %>%
  ungroup() %T>%
  
  #Check that the conversion to XY worked visually 
  plot(Y~X,data=.,pch=16)

#### Length Data ####
cope<-read_csv('../Data/COPE Lengths.csv') %>%
  
  #Remove Cloud 4 since it only has 3 individuals
  filter(Cloud!="COPE-M-04") %>%
  
  #Rename Length column
  rename(SL=`SL (mm)`) %>%
  
  #Rename shoals
  mutate(Cloud=recode(Cloud,`COPE-M-01`='A',`COPE-M-03`='B',`COPE-M-05`='C',`COPE-M-08`='D',`COPE-M-09`='E',`COPE-M-10`='F',
                      `COPE-M-11`='G',`COPE-M-12`='H',`COPE-M-13`='I',`COPE-M-14`='J',`COPE-M-15`='K',`COPE-M-18`='L',`COPE-M-20`='M')) %>%
  
  ## Merge into dataframe with full individual data ##
  #Get individual coordinates
  inner_join(shoal_locations,by='Cloud') %>%
  select(-date:-number) %>%
  
  #Get Microsat data
  inner_join(microsats,by=c('ID'='Sample')) %>%
  select(-Cloud.y) %>%
  rename(Shoal=Cloud.x)

#### Observed Heterozygosity ####
microsats <- cope %>%
  #Reformat for ADEGENET
  unite(COPE5,COPE5.a:COPE5.b,sep='/') %>%
  unite(COPE9,COPE9.a:COPE9.b,sep='/') %>%
  unite(COPE10,COPE10.a:COPE10.b,sep='/') %>%
  unite(CPER26,CPER26.a:CPER26.b,sep='/') %>%
  unite(CPER52,CPER52.a:CPER52.b,sep='/') %>%
  unite(CPER92,CPER92.a:CPER92.b,sep='/') %>%
  unite(CPER99,CPER99.a:CPER99.b,sep='/') %>%
  unite(CPER119,CPER119.a:CPER119.b,sep='/') %>%
  unite(CPER188,CPER188.a:CPER188.b,sep='/') %>%
  
  mutate(site='BZ14-M') %>%
  select(-SL:-Plate) %$% 
  
  #Read into ADEGENET
  df2genind(.[,c(-1,-2,-12)],sep='/',ind.names=ID,pop=site,NA.char =NA,type='codom')

#### Homozygote Excess ####
summary(microsats)$loc.n.all
summary(microsats)$Hobs
summary(microsats)$Hexp

excess.test<-hwx.test(microsats,B=10000,statName='U',method='monte')
excess.out<-hwdf(excess.test,statName='U',showN=T,showk=T)

#### Missing Data ####
cope %>%
  select(contains('.a')) %>%
  gather(Locus,Allele) %>%
  separate(Locus,c('Locus','Remove')) %>%
  select(-Remove) %>%
  group_by(Locus) %>%
  summarise(n=n(),missing=sum(is.na(Allele)),prop_miss=missing/n)

cope <- cope %>%
  select(-contains('CPER52'),-contains('COPE10'))

microsats <- microsats[loc = -which(locNames(microsats) %in% c('CPER52', 'COPE10'))]

marker_check<-cope %>%
  dplyr::select(ID, COPE5.a:CPER188.b) %>% 
  mutate(COPE5 = !is.na(COPE5.a),
         COPE9 = !is.na(COPE9.a),
         CPER26 = !is.na(CPER26.a),
         CPER92 = !is.na(CPER92.a),
         CPER99 = !is.na(CPER99.a),
         CPER119 = !is.na(CPER119.a),
         CPER188 = !is.na(CPER188.a)) %>%
  dplyr::select(-COPE5.a:-CPER188.b) 

cope <- cope %>%
  inner_join(marker_check %>%
               mutate(number_loci = select(., COPE5:CPER188) %>% apply(1, sum, na.rm=TRUE)) %>%
               dplyr::select(ID, number_loci),
             by = 'ID')
  
#### Select Relatedness Estimator ####
sim_rel<-cope %>%
  dplyr::select(-Shoal:-Plate) %>%
  dplyr::select(-starts_with('COPE10'),-starts_with('CPER52')) %>%
  readgenotypedata(.) %$%
  familysim(freqs,500) %>%
  coancestry(.,error.rates = 0.01,allow.inbreeding=T,
             trioml=1,wang=1,lynchli=1,lynchrd=1,ritland=1,quellergt = 1,dyadml = 1,
             trioml.num.reference=150, rng.seed = 12345) %$%
  cleanuprvals(relatedness,500) %>%
  mutate(group=if_else(group=='POPO','Parent-Offspring',
                       if_else(group=='SBSB','Full-Sib',
                               if_else(group=='HSHS','Half-Sib','Unrelated')))) %>%
  mutate(group=factor(group,levels = c('Parent-Offspring','Full-Sib','Half-Sib','Unrelated'))) %>%
  mutate(true_r=if_else(group=='Parent-Offspring' | group=='Full-Sib',0.5,
                        if_else(group=='Half-Sib',0.25,0))) %>%
  select(-pair.no:-ind2.id) %>%
  gather(metric,relatedness,-group,-true_r) %>%
  mutate(metric=if_else(metric=='dyadml','Milligan 2003',
                        if_else(metric=='lynchli','Li et al 1993',
                                if_else(metric=='lynchrd','Lynch and Ritland 1999',
                                        if_else(metric=='quellergt','Queller and Goodnight 1989',
                                                if_else(metric=='ritland','Ritland 1996',
                                                        if_else(metric=='trioml','Wang 2007',
                                                                if_else(metric=='wang','Wang 2002','Missed one'))))))))

# sim_rel<-read_csv('../../Results_tmp/simulated_relatedness.csv') %>%
#   mutate(group=factor(group,levels=c('Parent-Offspring','Full-Sib','Half-Sib','Unrelated'))) %>%
#   mutate(metric=if_else(metric=='dyadml','Milligan 2003',
#                         if_else(metric=='lynchli','Li et al 1993',
#                                 if_else(metric=='lynchrd','Lynch and Ritland 1999',
#                                         if_else(metric=='quellergt','Queller and Goodnight 1989',
#                                                 if_else(metric=='ritland','Ritland 1996',
#                                                         if_else(metric=='trioml','Wang 2007',
#                                                                 if_else(metric=='wang','Wang 2002','Missed one'))))))))

max_unrel <- sim_rel %>%
  filter(group == 'Unrelated',
         metric == 'Milligan 2003') %$%
  quantile(relatedness, 0.975)

#### Probability of identifying a true unrelated pair as unrelated ####
## Adapted from Wang 2007
N<-3
plan(multicore)
probability_exclustion<- microsats %>%
  genind2genpop %>%
  makefreq %>%
  t %>%
  as_tibble(rownames = 'loci.allele') %>%
  rename(frequency = `BZ14-M`) %>%
  separate(loci.allele, into = c('loci','allele')) %>%
  group_by(loci) %>%
  nest(.key = P) %>%
  
  #filter(loci != 'COPE5') %>%
  
  mutate(P = map(P, function(x) x$frequency),
         K = map_int(P, length)) %>%
  
  mutate(sigma.1 = future_map_dbl(P, ~.x^(2*N) %>% sum),
         sigma.2 = future_map2_dbl(P, K, function(p, k) tibble(u = 1:k, v = 1:k, w = 1:k, x = 1:k) %>%
                              tidyr::expand(u, v) %>%
                              filter(u!=v) %>%
                              mutate(u = p[u],
                                     v = p[v]) %>%
                              mutate(tmp = u+v) %$%
                              tmp %>%
                              raise_to_power(2*N) %>%
                              sum
         ),
         sigma.3 = future_map2_dbl(P, K, function(p, k) tibble(u = 1:k, v = 1:k, w = 1:k, x = 1:k) %>%
                              tidyr::expand(u, v, w) %>%
                              filter(u!=v,
                                     w!=v,
                                     w!=u) %>%
                              mutate(u = p[u],
                                     v = p[v], 
                                     w = p[w]) %>%
                              mutate(tmp = 
                                       
                                       3*(u*(u+2*v)+2*w*(u+v))^N-2*(u*(v+2*w)+v*(u+2*w))^N-3*u^N*((u+2*v)^N+(u+2*w)^N)-
                                       2^N*(w^N*(u+v)^N+v^N*(u+w)^N-2*u^N*(v+w)^N-(v*w)^N-(u*v)^N-(u*w)^N)   
                                     
                              ) %$%
                              tmp %>%
                              sum
         ),
         sigma.4 = future_map2_dbl(P, K, function(p, k) tibble(u = 1:k, v = 1:k, w = 1:k, x = 1:k) %>%
                              tidyr::expand(u, v, w, x) %>%
                              filter(u!=v,
                                     w!=v,
                                     w!=u,
                                     x!=u,
                                     x!=v,
                                     x!=w) %>%
                              mutate(u = p[u],
                                     v = p[v], 
                                     w = p[w],
                                     x = p[x]) %>%
                              mutate(tmp = 
                                       
                                       2*(u*v)^N+3*(u*w)^N-(v*x)^N+
                                       2*(w*x)^N+(x^N-3*u^N+(u+x)^N)*(v+w)^N+(v^N-3*w^N)*(u+x)^N+
                                       ((u*w)+(v*x))^N-((u*v)+(u*w)+(v*x))^N-2*((u*v)+(w*x))^N+
                                       3*((u*v)+(u*w)+(w*x))^N-((u*v)+(v*x)+(w*x))^N-
                                       ((u*w)+(v*x)+(w*x))^N   
                                     
                              ) %$%
                              tmp %>%
                              sum
         )
  ) %>%
  mutate(Se = 1-(1/2)*(K-2)*(K-3)*sigma.1-(1/2)*sigma.2-(1/6)*sigma.3-(2^(N-3))*sigma.4)

probability_exclustion %$%
  Se %>%
  subtract(1,.) %>%
  prod %>%
  subtract(1,.)
#0.9989585

#### Pairwise Relatedness ####
pairwise_related <- cope %>% #Did I forget to remove the two shit loci?
  select(-Shoal:-Plate) %>%
  
  #slice(sample(381,5)) %>% #Shrink for testing code
  
  coancestry(.,error.rates = 0.01,
             dyadml=2,
             allow.inbreeding=T,ci95.num.bootstrap=1000,
             rng.seed = 12345)  %$% #Run takes ~2.5 hours with full dataset
  inner_join(relatedness,relatedness.ci95,by='pair.no') %>%
  select(pair.no,ind1.id.x,ind2.id.x,
         dyadml,dyadml.low,dyadml.high) %>%
  rename(ind1.id=ind1.id.x, ind2.id=ind2.id.x) %>%
  
  # Set to unrelated if confidence interval includes 0
  mutate(dyadml=if_else(dyadml.low==0,0,dyadml)) %>%
  
  #Add in Individual A information (Shoals, Sizes, X, Y) 
  inner_join(cope,by=c('ind1.id'='ID')) %>%
  select(-COPE5.a:-CPER188.b) %>%
  
  #Add in Individual B information (Shoals, Sizes, X, Y) 
  inner_join(cope,by=c('ind2.id'='ID')) %>%
  select(-COPE5.a:-CPER188.b) %>%
  
  #Determine pairwise variables (same shoal, |size difference|, distance apart)
  mutate(same.shoal=Shoal.x==Shoal.y) %>%
  mutate(same.plate=Plate.x==Plate.y) %>%
  mutate(size.difference=abs(SL.x-SL.y)) %>%
  mutate(distance=sqrt((X.x-X.y)^2+(Y.x-Y.y)^2)) %>%
  filter(dyadml!=1) %>%
  inner_join(marker_check, by = c('ind1.id'='ID')) %>%
  inner_join(marker_check, by = c('ind2.id'='ID'), suffix = c('.1','.2')) %>%
  mutate(COPE5 = COPE5.1 == COPE5.2,
         COPE9 = COPE9.1 == COPE9.2,
         CPER26 = CPER26.1 == CPER26.2,
         CPER92 = CPER92.1 == CPER92.2,
         CPER99 = CPER99.1 == CPER99.2,
         CPER119 = CPER119.1 == CPER119.2,
         CPER188 = CPER188.1 == CPER188.2) %>%
  dplyr::select(-COPE5.1:-CPER188.2) %>%
  mutate(number_match = COPE5 + COPE9 + CPER26 + CPER92 + CPER99 + CPER119 + CPER188) %>%
  dplyr::select(-COPE5:-CPER188)

# pairwise_related<-read_csv('../../Results_tmp/pairwise_relatedness.csv') %>%
#   inner_join(marker_check, by = c('ind1.id'='ID')) %>%
#   inner_join(marker_check, by = c('ind2.id'='ID'), suffix = c('.1','.2')) %>%
#   mutate(COPE5 = COPE5.1 == COPE5.2,
#          COPE9 = COPE9.1 == COPE9.2,
#          CPER26 = CPER26.1 == CPER26.2,
#          CPER92 = CPER92.1 == CPER92.2,
#          CPER99 = CPER99.1 == CPER99.2,
#          CPER119 = CPER119.1 == CPER119.2,
#          CPER188 = CPER188.1 == CPER188.2) %>%
#   dplyr::select(-COPE5.1:-CPER188.2) %>%
#   mutate(number_match = COPE5 + COPE9 + CPER26 + CPER92 + CPER99 + CPER119 + CPER188) %>%
#   dplyr::select(-COPE5:-CPER188)

#### Marker power relatedness analysis ####
a_col<-1+seq(1,14,by=2)

make_dir<-function(x){
  dir_create(x$tmp_dir)
}

plan(multicore, workers = 20) #breaks running in parallel 
shared_microsat_relatedness<-tibble(number_shared = 1:7) %>% #Something broken with even numbers - no idea what/why c(seq(1,7,by=2), 6)
  mutate(a_cols = map(number_shared, function(x) combn(a_col, m=x, simplify = FALSE))) %>%
  arrange(-number_shared) %>%
  mutate(subset_data = map(a_cols, function(x) x %>%
                             map(function(y) cope %>%
                                   filter(number_loci == 7) %>%
                                   dplyr::select(-Shoal:-Plate,-number_loci) %>%
                                   dplyr::select(ID, sort(c(y, 1+y)))
                             ))) %>%
  unnest %>%
  group_by(number_shared) %>% 
  mutate(combo = 1:n()) %>%
  ungroup %>%
  mutate(tmp_dir = str_c(tempdir(),'/',number_shared,'_',combo)) %T>%
  make_dir %>%
  mutate(pairwise_relatedness = future_map2(subset_data, tmp_dir, function(x,y) x %>%
                                              coancestry(.,error.rates = 0.01,
                                                         dyadml=2, 
                                                         allow.inbreeding=T, ci95.num.bootstrap=1000,
                                                         rng.seed = 12345, working.directory = y) %$% 
                                              inner_join(relatedness,relatedness.ci95,by='pair.no') %>%
                                              # relatedness %>%
                                              
                                              dplyr::select(pair.no,ind1.id.x,ind2.id.x,
                                                            dyadml,dyadml.low,dyadml.high) %>%
                                              rename(ind1.id=ind1.id.x, ind2.id=ind2.id.x)
  )
  )  %>%
  dplyr::select(number_shared, combo, pairwise_relatedness) %>%
  unnest %T>%
  write_csv('../../Results_tmp/all_combos_prefiltering_microsats_relatedness.csv') 
plan(sequential)

shared_microsat_relatedness<-read_csv('../../Results_tmp/all_combos_prefiltering_microsats_relatedness.csv')

shared_microsat_summary <- shared_microsat_relatedness %>%
  mutate(dyadml.f1=if_else(dyadml.low==0,0,dyadml)) %>% #Filter 1
  mutate(dyadml.f2=if_else(dyadml < max_unrel, 0 ,dyadml.f1)) %>% #Filter 2
  mutate(dyadml.f3=if_else(dyadml.low < max_unrel, 0, dyadml.f1)) %>% #Filter 3
  mutate(pair = str_c(ind1.id, ind2.id, sep = '_')) %>%
  dplyr::select(-combo:-ind2.id, -dyadml.low, -dyadml.high) %>%
  gather(stringency, dyadml, -number_shared, -pair) %>%
  mutate(stringency = case_when(stringency == 'dyadml' ~ 'no_changes',
                                stringency == 'dyadml.f1' ~ 'no_0',
                                stringency == 'dyadml.f2' ~ 'no_UR',
                                stringency == 'dyadml.f3' ~ 'no_UR_overlap')) %>%
  group_by(pair, stringency, number_shared) %>%
  summarise(dyad.ind.mean = mean(dyadml)) %>%
  ungroup %>%
  group_by(stringency, number_shared) %>%
  summarise(related_mean = mean(dyad.ind.mean), sd_mean = sd(dyad.ind.mean), n = n()) %>%
  mutate(lwr95 = related_mean - 1.96 * sd_mean/n, upr95 = related_mean + 1.96 * sd_mean/n) 

shared_microsat_summary %>%
  ggplot(aes(x = number_shared, y = related_mean, ymin = lwr95, ymax = upr95, colour = stringency)) +
  geom_line(lty = 2) +
  geom_point() +
  geom_linerange()

tmp<-shared_microsat_summary %>%
  ungroup %>%
  add_row(number_shared = 1:7, stringency = 'no_0', .before = 0) %>%
  add_row(number_shared = 1:7, stringency = 'no_changes', .before = 0) %>%
  add_row(number_shared = 1:7, stringency = 'no_UR', .before = 0) %>%
  add_row(number_shared = 1:7, stringency = 'no_UR_overlap', .before = 0)


rel_common_res <- inla(related_mean ~ f(number_shared, model = 'rw2',
                                        hyper = list(theta = list(prior="pc.prec", param=c(1,0.1)))) + 
                         f(stringency, model = 'iid', 
                           hyper=list(theta=list(prior="loggamma", param=c(1,1)))),
                       
                       family = c("beta"), #change to beta for dyadml
                       #family = c('gaussian'),
                       data = tmp,
                       control.family = list(link='logit'),
                       control.fixed= list(prec.intercept = 1),
                       control.predictor = list(link = 1, compute = TRUE),
                       control.compute = list(dic = TRUE, 
                                              waic = TRUE, 
                                              cpo = TRUE),
                       verbose = TRUE)
summary(rel_common_res)


rel_common_res %$%
  summary.fitted.values %>%
  as_tibble() %>%
  bind_cols(tmp,.) %>%
  filter(is.na(related_mean)) %>%
  rename(q0.025 = `0.025quant`, q0.975 = `0.975quant`) %>%
  
  ggplot(aes(x = number_shared, y = mean, ymin = q0.025, ymax = q0.975, colour = stringency)) +
  geom_line(lty = 2) +
  geom_point() +
  geom_linerange()

rel_common_res$cpo

-mean(log(rel_common_res$cpo$cpo), na.rm=TRUE)



#### Make Graph ####
pairwise_filter<-function(x, y){
  #x: loci options
  #y: relatedness option
  pairwise_related %>%
    
    select(-X1,-pair.no) %>%
    
    #Remove unrelated edges
    #filter(dyadml!=0) %>%
    
    ## Loci filtering options
    {if(x == 'only_7') filter(., number_match == 7) else(.)} %>%
    
    #Relatedness filtering
    #{if(y == 'no_mean_UR') filter(., dyadml > max_unrel) else(.)} %>%
    #{if(y == 'no_overlap_UR') filter(., dyadml.low > max_unrel) else(.)} 
    {if(y == 'no_mean_UR') mutate(., dyadml = if_else(dyadml > max_unrel, dyadml, 0)) else(.)} %>%
    {if(y == 'no_overlap_UR') mutate(., dyadml = if_else(dyadml.low > max_unrel, dyadml, 0)) else(.)} 
}

make_the_graph<-function(x){
  x %>%
    
    filter(dyadml!=0) %>%
    
    #Set relatedness as edge weights
    rename(weight=dyadml) %>% 
    
    #Make Graph
    graph_from_data_frame(d=., vertices=cope, directed=F) %T>%
    
    #Plot to check everything went through right
    plot(.,vertex.label=NA,vertex.size=4)
    
}

#Loci Options
#All, only_7

#Relatedness options
#all, mean_in_UR, overlap_UR

BOOT<-1000
a<-1;b<-1
plan('multiprocess')

cope_graphs<-tibble(loci_option = c('all', 'only_7', NA), related_option = c('no_0','no_mean_UR','no_overlap_UR')) %>%
  tidyr::expand(related_option,loci_option) %>%
  na.omit %>%
  mutate(pairwise = map2(loci_option, related_option, pairwise_filter),
         graphs = map(pairwise, make_the_graph)) %>%
  
  #Summarize Graphs
  mutate(related_dyads = map_int(graphs, function(x) E(x) %>% length),
         percent_related_all_dyads = map_dbl(graphs, function(x) (E(x) %>% length)/(V(x) %>% length %>% choose(.,2))),
         number_without_relatives = map_int(graphs, function(x) sum(degree(x)==0)),
         number_with_a_relative = map_int(graphs, function(x) sum(degree(x) > 0)),
         rel_sum = map(pairwise, function(x) x %>%
                         filter(dyadml > 0) %>%
                         summarise(n = n(), mean.rel = mean(dyadml), sd.rel = sd(dyadml)))) %>%
  unnest(rel_sum) %>%
  
  #Build appropriate random graphs
  mutate(ER_random = map(graphs, function(x) future_map(1:BOOT, function(y) play_erdos_renyi(n=vcount(x),
                                                                                             p=ecount(x)/choose(vcount(x),2),
                                                                                             directed=F)) %>%
                           map(.,set_vertex_attr,name='Shoal',value=V(x)$Shoal))) %>%
  
  #Compare observed to similar random graphs
  mutate(obs_transitivity = map_dbl(graphs, transitivity),
         obs_assortativity = map_dbl(graphs, function(x) assortativity(x,as.factor(V(x)$Shoal)))) %>%
  mutate(ER_transitivity = map(ER_random, function(y) map_dbl(y, transitivity)),
         ER_assortativity = map2(graphs, ER_random, function(x, y) map_dbl(y,assortativity,types1=as.factor(V(x)$Shoal)))) %>%
  mutate(sim_obs_transitivity = map2(obs_transitivity, ER_transitivity, function(x, y) y %>%
                                       quantile(c(0.025,0.975), na.rm = T) %>%
                                       unname %>%
                                       enframe() %>%
                                       mutate(level = c('lwr_0.025_trans', 'upr_0.975_trans')) %>%
                                       dplyr::select(-name) %>%
                                       spread(level, value) %>%
                                       mutate(p_trans = sum(x < y)/BOOT)),
         sim_obs_assortativity = map2(obs_assortativity, ER_assortativity, function(x, y) y %>%
                                        quantile(c(0.025,0.975), na.rm = T) %>%
                                        unname %>%
                                        enframe() %>%
                                        mutate(level = c('lwr_0.025_assort', 'upr_0.975_assort')) %>%
                                        dplyr::select(-name) %>%
                                        spread(level, value) %>%
                                        mutate(p_assort = sum(x < y)/BOOT))) %>%
  unnest(sim_obs_transitivity, sim_obs_assortativity) %>%
  dplyr::select(-starts_with('ER'))

## Number with a relative
prop_with_rel<-cope_graphs %>%
  dplyr::select(loci_option, related_option, number_without_relatives, number_with_a_relative) %>%
  mutate(n = number_with_a_relative,
         total = number_without_relatives + number_with_a_relative,
         Prop_related=(n+a)/(total+a+b),
         sd_prop_intergen=((a+n)*(b+total-n))/((a+b+total)^2*(a+b+total+1)),
         PI_25=qbeta(0.25,n+a,total-n+b),
         PI_75=qbeta(0.75,n+a,total-n+b),
         PI_2.5=qbeta(0.025,n+a,total-n+b),
         PI_97.5=qbeta(0.975,n+a,total-n+b))

fig0<-prop_with_rel %>%
  
  ggplot(aes(x = related_option, y = Prop_related, colour = loci_option)) +
  geom_linerange(aes(ymin = PI_2.5, ymax = PI_97.5), lty = 2, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = PI_25, ymax = PI_75), position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) 
ggsave(filename = '../Figures/Figure 0.tiff',plot = print(fig0),scale=1.5,height = 6,width = 6,dpi = 300)


## Average relatedness 
mean_rel<-cope_graphs %>%
  dplyr::select(loci_option, related_option, n, mean.rel, sd.rel) %>%
  mutate(se.rel = sd.rel/sqrt(n),
         lwr = mean.rel - 1.96*se.rel,
         upr = mean.rel + 1.96*se.rel) 


fig0b<-mean_rel %>%
  ggplot(aes(x = related_option, y = mean.rel, colour = loci_option)) +
  geom_linerange(aes(ymin = lwr, ymax = upr), position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) 
ggsave(filename = '../Figures/Figure 0b.tiff',plot = print(fig0b),scale=1.5,height = 6,width = 6,dpi = 300)

cope_graphs %>% 
  dplyr::select(-graphs, -pairwise)

## Transitivity 
cope_graphs %>%
  dplyr::select(related_option, loci_option, obs_transitivity, lwr_0.025_trans, upr_0.975_trans, p_trans) %T>%
  print %>%
  
  ggplot(aes(x = related_option, y = obs_transitivity, colour = loci_option)) +
  geom_linerange(aes(ymin = lwr_0.025_trans, ymax = upr_0.975_trans), position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5))

## Assortativity 
cope_graphs %>%
  dplyr::select(related_option, loci_option, obs_assortativity, lwr_0.025_assort, upr_0.975_assort, p_assort)%T>%
  print %>%
  
  ggplot(aes(x = related_option, y = obs_assortativity, colour = loci_option)) +
  geom_linerange(aes(ymin = lwr_0.025_assort, ymax = upr_0.975_assort), position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5))

#### Shoal mean relatedness ####
perm_rel<-function(x,type){
  as_tibble(x) %>%
    select(contains(type),dyadml,ind1.id,ind2.id) %>%
    rename_(P1 = names(.)[1], P2=names(.)[2]) %>%
    mutate(same_perm=P1==P2) %>%
    filter(same_perm==T) %>%
    group_by(P1) %>%
    summarise(individuals=length(unique(c(ind1.id,ind2.id))),total_dyads=n(),related.dyads=sum(dyadml>0),
              mean.rel=mean(dyadml),sd.rel=sd(dyadml)) %>%
    mutate(se.rel=sd.rel/sqrt(total_dyads)) %>%
    rename(!!type := P1)
}

full_permutation_related<-function(x,type,NPERM=1000){
  perms<-x %>%
    modelr::permute(NPERM,dyadml) %>%
    group_by(.id) %>%
    do(map_df(.$perm,perm_rel,type=type)) %>%
    ungroup() %>%
    group_by_(type) %>%
    do(lwr=quantile(.$mean.rel,0.025),upr=quantile(.$mean.rel,0.975),full=.$mean.rel) %>%
    unnest(lwr,upr) %>%
    rowwise() %>%
    mutate(n=length(full))
  
  x %>%
    perm_rel(type=type) %>%
    inner_join(perms,by=type) %>%
    rowwise() %>%
    mutate(p=sum(mean.rel<full)/n) %>%
    mutate(p=ifelse(p==0,1/n,p)) %>%
    ungroup()
}

shoal_relatedness<-cope_graphs %>%
  dplyr::select(loci_option, related_option, pairwise) %>%
  mutate(shoal_related = future_map(pairwise, full_permutation_related, 'Shoal', 1000)) %>%
  dplyr::select(-pairwise) %>%
  unnest

#### Spatial strucuring of shoals  ####


#### Cross generational relations ####
GENERATION<-9 #Generation ~9 mm SL per Cole & Robertson 1988

a<-1;b<-1

intergenerational_analysis<-cope_graphs %>%
  mutate(intergenerational_relatives = map(graphs, function(x) x %>%
                                   as_tbl_graph() %>%
                                   activate(edges) %>%
                                   mutate(gens_between=size.difference %/% GENERATION) %>%
                                   as_tibble() %>%
                                   bind_rows((x %>%
                                                complementer %>%
                                                as_tbl_graph() %>%
                                                activate(edges) %>%
                                                mutate(size.difference=abs(.N()$SL[from]-.N()$SL[to]))%>%
                                                mutate(gens_between=size.difference %/% GENERATION) %>%
                                                as_tibble()),
                                             .id='Related') %>%
                                   mutate(Related=if_else(Related=='1','Relatives','Unrelated')) %>%
                                   group_by(Related,gens_between) %>%
                                   summarise(n=n()) %>%
                                   ungroup %>%
                                   group_by(Related) %>%
                                   mutate(total=sum(n)) %>%
                                   ungroup %>%
                                   filter(gens_between==1)  %>%
                                     mutate(Prop_intergen=(n+a)/(total+a+b),
                                            sd_prop_intergen=((a+n)*(b+total-n))/((a+b+total)^2*(a+b+total+1)),
                                            PI_25=qbeta(0.25,n+a,total-n+b),
                                            PI_75=qbeta(0.75,n+a,total-n+b),
                                            PI_2.5=qbeta(0.025,n+a,total-n+b),
                                            PI_97.5=qbeta(0.975,n+a,total-n+b))),
         intergeneration_rel.v.nonrel = map(intergenerational_relatives, function(x) x %>%
                                              ungroup %>%
                                              summarise(difference=diff(Prop_intergen),
                                                        sd_diff=sqrt(sum(sd_prop_intergen)),
                                                        diff_2.5=qnorm(0.025,difference,sd_diff),
                                                        diff_97.5=qnorm(0.975,difference,sd_diff)))) %>%
  dplyr::select(related_option, loci_option, intergenerational_relatives, intergeneration_rel.v.nonrel)

intergenerational_analysis %>%
  dplyr::select(-intergeneration_rel.v.nonrel) %>%
  unnest %T>%
  print %>%
  ggplot(aes(x=interaction(related_option,loci_option), y = Prop_intergen, colour = Related)) +
  geom_linerange(aes(ymin = PI_2.5, ymax = PI_97.5), lty = 2, position = position_dodge(width = 0.2)) +
  geom_linerange(aes(ymin = PI_25, ymax = PI_75), position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2))

intergenerational_analysis %>%
  dplyr::select(-intergenerational_relatives) %>%
  unnest %T>%
  print %>%
  
  ggplot(aes(x = related_option, y = difference, colour = loci_option, ymin = diff_2.5, ymax = diff_97.5)) +
  geom_linerange(position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2)) +
  geom_hline(yintercept = 0, lty = 2)

#### Figures ####
## Figure 1 ##
#Figure showing each shoal with "flow" between shoals showing the number and strength of relatives
fig1<-cope_graphs %>%
  dplyr::select(related_option, loci_option, graphs) %>%
  mutate(figures = map(graphs, function(X) X %>%
                         as_tbl_graph %>%
                         dplyr::select(-Plate:-CPER188.b) %>%
                         activate(edges) %>%
                         dplyr::select(-dyadml.low:-size.difference) %>%
                         activate(nodes) %>%
                         mutate(Shoal = as.factor(Shoal)) %>%
                         contract(., mapping = V(.)$Shoal,
                                  vertex.attr.comb = list(X = 'mean',
                                                          Y = 'mean',
                                                          Shoal = function(x)x[1], "ignore",
                                                          SL = function(x) length(x))) %>%
                         simplify(., remove.multiple = TRUE, remove.loops = FALSE,
                                  edge.attr.comb = list(weight = "sum", 
                                                        distance = 'mean',
                                                        number_match = function(x) length(x))) %>%
                         as_tbl_graph %>%
                         rename(x=X, y=Y, N = SL) %>%
                         mutate(x = x/max(x), y=y/max(y)) %>%
                         activate(edges) %>%
                         rename(number_dyads = number_match) %>%
                         mutate(weight = weight/number_dyads) %>%
                         activate(nodes) %>%
                         mutate(centrality = centrality_authority(weights = weight)) %>%
                         
                         ggraph(layout='manual', node.positions = as_tibble(.)) +
                         geom_edge_link(aes(colour = weight, alpha = weight), show.legend = FALSE) + #
                         geom_edge_loop(aes(colour = weight, alpha = weight), show.legend = FALSE) + #
                         geom_node_point(aes(colour = centrality, size = N), show.legend = FALSE) +
                         #coord_fixed() +
                         scale_x_continuous('X (cm)') +
                         scale_y_continuous('Y (cm)') +
                         theme_classic()
                         )) %$%
  plot_grid(plotlist = figures, labels = LETTERS[1:6], ncol = 2)
ggsave(filename = '../Figures/Figure 1.tiff',plot = print(fig1),scale=1.5,height = 6,width = 4,dpi = 300)

## Figure 2 & Table 2 ##
sim_rel %>%
  group_by(metric,group) %>% 
  mutate(M=mean(relatedness),S=sd(relatedness),CV=S/M) %>%
  ungroup %>%
  group_by(metric) %>%
  mutate(mean_cv=mean(CV),w=cor(true_r,relatedness)) %>%
  select(metric,mean_cv,w,group,M,S,CV) %>%
  distinct() %>%
  arrange(-w) %>%
  write_csv('../Figures/Table 2.csv')

fig2<-sim_rel %>%
  group_by(metric,group) %>%
  summarise(mean_r=mean(relatedness),median_r=median(relatedness),
            sd_r=sd(relatedness),lwr_r=quantile(relatedness,0.25),upr_r=quantile(relatedness,0.75),
            true_r2=mean(true_r)) %>%
  
  ggplot(aes(x=group,y=median_r,ymin=lwr_r,ymax=upr_r,colour=metric)) +
  geom_errorbar(aes(ymax=true_r2,ymin=true_r2), colour='black', lty=2) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5),width=0) +
  coord_cartesian(ylim = c(0, 1)) +
  xlab('True relationship') +
  ylab('Simulated Relatedness') +
  scale_color_discrete('Relatedness Estimator') +
  theme_classic()
ggsave(filename = '../Figures/Figure 2.eps',plot = print(fig2),scale=1.5,height=4,width = 4,device = cairo_ps,dpi = 300)

## Figure 3 ##
cope.layout<-cope_graphs$graphs[[1]] %>%
  as_tbl_graph() %>%
  create_layout(layout = 'nicely') %>%
  select(name,x,y)


fig3<-cope_graphs %>%
  dplyr::select(related_option, loci_option, graphs) %>%
  mutate(figures = map2(graphs, loci_option, function(X, A) X %>%
                         as_tbl_graph() %>%
                         dplyr::select(-Plate:-CPER188.b) %>%
                          inner_join(cope.layout %>% as_tibble %>% mutate(name = as.character(name)),
                                     by = 'name') %>%
                         activate(edges) %>%
                         mutate(gens_between=size.difference %/% GENERATION) %>%
                         mutate(relationship_type=if_else(gens_between==0,'Within Cohort','Intergenerational')) %>%
                         mutate(relationship_type = factor(relationship_type,levels=c('Within Cohort','Intergenerational'))) %>%
                         activate(nodes) %>%
                         {if(A == 'only_7') filter(., number_loci == 7) else(.)} %>%
                          
                         ggraph(layout='manual', node.positions = as_tibble(.)) +
                         #ggraph(layout='nicely') +
                         geom_edge_link(aes(colour=relationship_type, width = weight), show.legend = FALSE) + #
                         geom_node_point(size = 0.5) +
                         #scale_edge_alpha_continuous(range=c(0.4,1)) +
                         scale_edge_width_continuous(range = c(0.05, 0.5)) +
                         coord_fixed() +
                         theme_classic() + 
                         theme(panel.border = element_blank(), 
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_blank(),
                               axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
                               strip.background = element_blank(),
                               strip.text.x = element_blank()))) %$%
  plot_grid(plotlist = figures, labels = LETTERS[1:6], ncol = 2)
ggsave(filename = '../Figures/Figure 3.tiff',plot = print(fig3),scale=1.5,height=4,width = 4,dpi = 300)

## Figure 4
distributions<-shoal_relatedness %>%
  select(loci_option, related_option, Shoal, full) %>%
  unnest() %>%
  rename(mean.rel=full)


shoal_relatedness %>%
  mutate(sig = if_else(p < 0.05, 'sig', 'not')) %>%
  
  ggplot(aes(x = Shoal, y = mean.rel)) +
  geom_violin(data=distributions) +
  geom_linerange(aes(ymin=mean.rel-se.rel,ymax=mean.rel+se.rel, colour = sig)) +
  geom_point(aes(colour = sig)) +
  facet_grid(related_option~loci_option) +
  theme_classic() +
  scale_y_continuous(expression(paste("Mean Pairwise Relatedness (",bar(italic(" r")),')'))) +
  scale_x_discrete('Shoal or Cluster identity')

