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
library(betapart)
library(spdep)
library(adespatial)
library(vegan)
library(RColorBrewer)
library(randomcoloR)
library(Cairo)

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

#### Pairwise Relatedness ####
pairwise_related <- cope %>%
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
  filter(dyadml!=1)

# pairwise_related<-read_csv('../../Results_tmp/pairwise_relatedness.csv')

#### Make Graph ####
cope.graph <- pairwise_related %>%
  
  select(-X1,-pair.no) %>%
  
  #Remove unrelated edges
  filter(dyadml!=0) %>%

  #Set relatedness as edge weights
  rename(weight=dyadml) %>%
  
  #Make Graph
  graph_from_data_frame(d=., vertices=cope, directed=F) %T>%
  
  #Plot to check everything went through right
  plot(.,vertex.label=NA,vertex.size=4)


## Summarize graph ##
E(cope.graph) %>% length

(E(cope.graph) %>% length)/(V(cope.graph) %>% length %>% choose(.,2))

sum(degree(cope.graph)==0)

sum(degree(cope.graph)==0)/(length(V(cope.graph)))

#### Make Random Graphs ####
BOOT<-1000
ER_graphs<-replicate(BOOT,play_erdos_renyi(n=vcount(cope.graph),p=ecount(cope.graph)/choose(vcount(cope.graph),2),directed=F),simplify = F) %>%
  map(.,set_vertex_attr,name='Shoal',value=V(cope.graph)$Shoal)


#### Transitivity ####
transitivity(cope.graph)
sum(transitivity(cope.graph) < map_dbl(ER_graphs,transitivity))/BOOT

map_dbl(ER_graphs,transitivity)  %>%
  tibble(boot.obs='boot',transitivity=.) %>%
  add_row(boot.obs = 'obs', transitivity = transitivity(cope.graph)) %>%
  ggplot(aes(x=transitivity,..density..,fill=boot.obs,col=boot.obs,group=boot.obs))+
  geom_histogram() +
  theme_classic()


#### Assortativity ####
assortativity(cope.graph,as.factor(V(cope.graph)$Shoal))
sum(assortativity(cope.graph,types1=as.factor(V(cope.graph)$Shoal)) < map_dbl(ER_graphs,assortativity,types1=as.factor(V(cope.graph)$Shoal)))/BOOT

map_dbl(ER_graphs,assortativity,types1=as.factor(V(cope.graph)$Shoal))  %>%
  tibble(boot.obs='boot',assortativity=.) %>%
  add_row(boot.obs = 'obs', assortativity = assortativity(cope.graph,types1=as.factor(V(cope.graph)$Shoal))) %>%
  ggplot(aes(x=assortativity,..density..,fill=boot.obs,col=boot.obs,group=boot.obs))+
  geom_histogram() +
  theme_classic()

#### Community detection ####
cope_community <- cope.graph %>% cluster_louvain()

sizes(cope_community) %>%
  as.tibble %>%
  ggplot(aes(x=n,..density..)) +
    geom_histogram() +
    theme_classic()

sum(sizes(cope_community) > 1)

mean(sizes(cope_community)[sizes(cope_community) > 1])

sd(sizes(cope_community)[sizes(cope_community) > 1])

range(sizes(cope_community)[sizes(cope_community) > 1])

cope_community %>% modularity()
sum((cope_community %>% modularity()) < (map(ER_graphs,cluster_louvain) %>% map_dbl(modularity)))/BOOT

map(ER_graphs,cluster_louvain)  %>%
  map_dbl(modularity) %>%
  tibble(boot.obs='boot',modularity=.) %>%
  add_row(boot.obs = 'obs', modularity = (cope_community %>% modularity())) %>%
  ggplot(aes(x=modularity,..density..,fill=boot.obs,col=boot.obs,group=boot.obs))+
  geom_histogram() +
  theme_classic()

cope %<>%
  mutate(community= paste('C',(cope_community %>% membership),sep=''))

pairwise_related %<>%
  inner_join((cope %>% select('ID','community')),by=c('ind1.id'='ID')) %>%
  inner_join((cope %>% select('ID','community')),by=c('ind2.id'='ID')) %>%
  mutate(same.community=community.x==community.y)

#### Shoal v Community Relatedness ####
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


group_relatedness<-bind_rows(full_permutation_related(pairwise_related,'Shoal',1000),
                             full_permutation_related(pairwise_related,'community',1000)) %>%
  mutate(s.c=ifelse(!is.na(Shoal),'Shoal','Community'))

group_relatedness %>%
  mutate(sig=p<0.05) %>% 
  group_by(s.c, sig) %>%
  summarize(mean_rel=mean(mean.rel),sd_rel=sd(mean.rel),min.p=min(p),n=n())

group_relatedness %>%
  filter(s.c=='Shoal') %>%
  dplyr::select(-full,-n,-community,-s.c) %>%
  filter(p<0.05) %>%
  arrange(Shoal) %>%
  dplyr::select(Shoal,mean.rel,p)

group_relatedness %>% 
  aov(mean.rel ~ s.c, data=.) %>%
  summary()

group_relatedness %>%
  filter(s.c=='Community') %>%
  ggplot(aes(x=individuals,y=mean.rel,ymin=mean.rel-se.rel,ymax=mean.rel+se.rel)) +
  geom_point() +
  geom_linerange() +
  ylab('Mean Relatedness') +
  xlab('Community Size') +
  theme_classic()

group_relatedness %>%
  mutate(sig=p<0.05) %>% 
  filter(s.c=='Community',individuals > 10) %>%
  group_by(sig) %>%
  summarize(mean_rel=mean(mean.rel),sd_rel=sd(mean.rel),min.p=min(p),n=n())


#### Spatial strucuring of shoals  ####
cope %>% 
  group_by(Shoal,community) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  arrange(n) %>%
  ggplot(aes(x=Shoal,y=community,fill=n)) + 
  geom_bin2d()

shoal_composition<-cope %>% 
  group_by(Shoal,community) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  spread(community,n,fill=0L) %>%
  right_join(shoal_locations,by=c('Shoal'='Cloud')) %>%
  select(Shoal,X,Y,starts_with('C')) %>%
  filter(!is.na(C1))

shoal_order<-unique(shoal_composition$Shoal)

shoal_XY <- shoal_composition %>%
  select(X,Y) %>%
  as.matrix()

shoal_composition %<>%
  select(-X,-Y,-Shoal)

## Global Bray-Curtis ##
beta.multi.abund(shoal_composition,index.family='bray')

## Make sphere of influence neighborhood based on shoal proximity ##
shoal_sphere<-shoal_XY %>%
  as.matrix() %>%
  tri2nb(coords=.) %>%
  soi.graph(tri.nb=.,coords=shoal_XY) %>%
  graph2nb(gob=.,sym=T) %T>%
  plot(.,coords=shoal_XY)

## Weight edges based on distance between shoals ##
shoal_weights<-shoal_sphere %>%
  nbdists(., shoal_XY) %>%
  map(., function(x) 1-x/max(dist(shoal_XY))) %>%
  nb2listw(shoal_sphere, glist = ., style = "B")

## Moran's Eigen Mapping ##
shoal_mem<-mem(shoal_weights) 
moranI <- moran.randtest(shoal_mem, shoal_weights, 10000, 'two-sided')

plot(shoal_mem[,which(moranI$pvalue<0.05)], SpORcoords = shoal_XY, nb = shoal_sphere)

## Redundancy analysis using significant Eigen mappings ##
shoal_rda<-rda(shoal_composition~.,data=shoal_mem[,which(moranI$pvalue<0.05)])
anova(shoal_rda,permutations=10000)
plot(shoal_rda)

#### Cross generational relations ####
GENERATION<-9 #Generation ~9 mm SL per Cole & Robertson 1988

a<-1;b<-1
cope.graph %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  mutate(gens_between=size.difference %/% GENERATION) %>%
  as_tibble() %>%
  bind_rows((cope.graph %>%
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
  filter(gens_between==1) %>%
  mutate(Prop_intergen=(n+a)/(total+a+b),
         sd_prop_intergen=((a+n)*(b+total-n))/((a+b+total)^2*(a+b+total+1)),
         PI_25=qbeta(0.25,n+a,total-n+b),
         PI_75=qbeta(0.75,n+a,total-n+b),
         PI_2.5=qbeta(0.025,n+a,total-n+b),
         PI_97.5=qbeta(0.975,n+a,total-n+b)) %T>%
  print %>%
  ungroup %>%
  summarise(difference=diff(Prop_intergen),
            sd_diff=sqrt(sum(sd_prop_intergen)),
            diff_2.5=qnorm(0.025,difference,sd_diff),
            diff_97.5=qnorm(0.975,difference,sd_diff))

#### Figures ####
## Figure 1 ##
#Modified from https://stackoverflow.com/questions/43984614/rggplot2geom-points-how-to-swap-points-with-pie-charts

table(sizes(cope_community))

sum(sizes(cope_community)>4)

comm_fills <- c(distinctColorPalette(sum(sizes(cope_community)>4)),'grey50','black')

plot(1:(sum(sizes(cope_community)>4)+2),1:(sum(sizes(cope_community)>4)+2),col=comm_fills,pch=16,cex=5)

g<-shoal_weights %>% 
  as(., "symmetricMatrix") %>% 
  graph.adjacency(., mode="undirected",weighted = T) 
V(g)$Shoal<-shoal_order
cope %>% dplyr::select(Shoal) %>% distinct %$% Shoal

small_comm<-str_c('C',(1:length(cope_community)),sep='')[sizes(cope_community)<10 & sizes(cope_community)!=1]
singles<-str_c('C',(1:length(cope_community)),sep='')[sizes(cope_community)==1]


plot_data<-cope.graph %>%
  as_tbl_graph %>%
  mutate(community= paste('C',(cope_community %>% membership),sep='')) %>% 
  mutate(community_simp=ifelse(community %in% singles,'Singleton',community)) %>%
  mutate(community_simp=ifelse(community %in% small_comm,'Clusters with 2-4 individuals',community_simp)) %>%
  as_tibble() %>%
  dplyr::select(-Plate:-CPER188.b) %>%
  group_by(Shoal,community_simp) %>%
  summarise(n=n(),x=mean(X),y=mean(Y)) %>%
  ungroup() %>%
  mutate(community_simp=as.factor(community_simp)) %>%
  mutate(community_simp=c(str_c('C',1:13,sep=''),'Clusters with 2-4 individuals','Singleton')[as.integer(community_simp)]) %>%
  mutate(community_simp=factor(community_simp,c(str_c('C',1:13,sep=''),'Clusters with 2-4 individuals','Singleton'))) %>%
  group_by(Shoal) %>%
  mutate(total=sum(n)) %>%
  mutate(x_circle=if_else(Shoal=='A' | Shoal=='H',x-110,if_else(Shoal=='B' | Shoal=='E',x+110,x))) %>%
  group_by(Shoal,x,y,total,x_circle)


names(comm_fills)<-c(str_c('C',1:13,sep=''),'Clusters with 2-4 individuals','Singleton')
pie.grobs <- plot_data %>% 
  do(pies = ggplot(., aes(1, n, fill = community_simp)) + 
       geom_col(position = "fill", alpha = 0.5, colour = NA) + 
       scale_fill_manual(values=comm_fills) + 
       coord_polar(theta = "y") + 
       theme_void()+ guides(fill = F))  %>%
  mutate(pie_grobs = list(annotation_custom(ggplotGrob(pies),
                                           xmin = x_circle-100, ymin = y-100,
                                           xmax = x_circle+100, ymax = y+100)))
back.grobs <- plot_data %>% 
  do(backs = ggplot(., aes(1, n)) + 
       geom_col(fill='white',position = "fill", alpha = 1, colour = NA) + 
       scale_fill_manual(values=comm_fills) + 
       coord_polar(theta = "y") + 
       theme_void()+ guides(fill = F))  %>%
  mutate(back_grobs = list(annotation_custom(ggplotGrob(backs),
                                            xmin = x_circle-100, ymin = y-100,
                                            xmax = x_circle+100, ymax = y+100))) %>%
  dplyr::select(Shoal,back_grobs)

fig1<-g %>% 
  as_tbl_graph %>%
  inner_join(pie.grobs,by=c('Shoal')) %>%
  inner_join(back.grobs,by=c('Shoal')) %>%
  mutate(x_label=if_else(Shoal=='D' | Shoal=='J' | Shoal=='C' | Shoal=='A' | Shoal=='H' | Shoal=='K',
                         x_circle-125,x_circle+125)) %>%
  
  ggraph(layout='manual',node.positions=(as.tibble(.) %>% dplyr::select(.,x,y))) +
  geom_edge_link() +
  scale_x_continuous('X (cm)') +
  scale_y_continuous('Y (cm)') +
  coord_fixed ()+
  theme_classic() +
  geom_node_point(size=3) +
  back.grobs$back_grobs +
  pie.grobs$pie_grobs +
  geom_node_text(aes(x=x_label,y=y,label = Shoal)) + 
  geom_col(data = plot_data,
           aes(0,0, fill = community_simp), 
           colour = "white", alpha = 0.5) +
  scale_fill_manual('Cluster',values=comm_fills) +
  theme(legend.position = 'bottom')
ggsave(filename = '../Figures/Figure 1.eps',plot = print(fig1),scale=1.5,height = 6,width = 4,device = cairo_ps,dpi = 300)

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
cope.layout<-cope.graph %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  mutate(cross=crossing(cope_community,cope.graph)) %>%
  filter(!cross) %>%
  create_layout(layout = 'nicely') %>%
  select(name,x,y)


plot_data<-cope.graph %>%
  as_tbl_graph() %>%
  dplyr::select(-Plate:-CPER188.b) %>%
  mutate(community= paste('C',(cope_community %>% membership),sep='')) %>%
  mutate(community=ifelse(community %in% singles,'Singleton',community)) %>%
  mutate(community=ifelse(community %in% small_comm,'Little family',community)) %>%
  mutate(community=as.factor(community)) %>%
  mutate(community=c(str_c('C',1:13,sep=''),'Small Family','Singleton')[as.integer(community)]) %>%
  mutate(community=factor(community,c(str_c('C',1:13,sep=''),'Small Family','Singleton'))) %>%
  activate(edges) %>%
  mutate(cross=crossing(cope_community,cope.graph)) %>%
  mutate(gens_between=size.difference %/% GENERATION) %>%
  mutate(relationship_type=if_else(gens_between==0,'Within Cohort','Intergenerational')) %>%
  mutate(relationship_type = factor(relationship_type,levels=c('Within Cohort','Intergenerational'))) %>%
  filter(!cross) %>%
  activate(nodes) %>%
  inner_join(cope.layout,by='name')


community_data <- plot_data %>%
  as_tibble() %>% 
  filter(community!='Small Family', community!='Singleton') %>%
  group_by(community) %>%
  mutate(hull=1:n(),hull=as.numeric(factor(hull, chull(x, y)))) %>%
  arrange(hull) %>%
  filter(!is.na(hull)) %>%
  mutate(polygon_type=n_distinct(hull)) %>%
  filter(polygon_type>1)
  

dat_text <- data.frame(
  label = c("A", "B"),
  relationship_type   = c('Within Cohort', 'Intergenerational')
)

fig3<-plot_data %>% 
  ggraph(layout='manual',node.positions=cope.layout) +
  geom_polygon(data = community_data,aes(x=x,y=y,fill=community,group=community), alpha = 0.5) +
  #geom_edge_link(aes(colour=relationship_type,alpha=weight)) + 
  geom_edge_link(aes(colour=weight)) + 
  geom_node_text(aes(label=Shoal),size=3) +
  geom_text(data=dat_text, aes(x = -Inf, y = Inf, label = label), size=10, hjust = -0.2, vjust=1.2) +
  facet_wrap(~relationship_type,nrow=1,dir='h') +
  scale_fill_manual(values=comm_fills) + 
  #scale_edge_alpha_continuous(range=c(0.4,1)) +
  scale_edge_color_continuous('Relatedness',low='gray80',high='gray20') +
  guides(shape=F,fill=F)+
  coord_fixed() +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
ggsave(filename = '../Figures/Figure 3.eps',plot = print(fig3),scale=1.5,height=4,width = 8,device = cairo_ps,dpi = 300)

## Figure 4 ##
dat_text <- data.frame(
  name = c("A", "B"),
  s.c   = c('Shoal', 'Community'),
  mean.rel = Inf,
  se.rel=NA,
  label=-Inf
)

plot_data<-group_relatedness %>%
  group_by(s.c) %>%
  mutate(index=1:n()) %>%
  ungroup() %>%
  mutate(s.c=factor(s.c, levels=c('Shoal','Community'))) %>% 
  filter((individuals > 10 & s.c=='Community') | s.c=='Shoal') %>%
  mutate(community=as.factor(community)) %>%
  mutate(community=str_c('C',1:length(unique(.$community)),sep='')[as.integer(community)]) %>%
  mutate(community=factor(community,str_c('C',1:length(unique(.$community)),sep=''))) %>%
  mutate(label=ifelse(is.na(Shoal),as.character(community),Shoal)) %>%
  mutate(label=factor(label,
                      levels=c(LETTERS[1:13],str_c('C',1:length(unique(.$community)),sep=''))))
  

distributions<-plot_data %>%
  select(label,s.c,full) %>%
  unnest() %>%
  rename(mean.rel=full)


fig4<-ggplot(plot_data,aes(x=label,y=mean.rel))+
  #geom_blank()+
  #geom_ribbon(data=plot_data,aes(x=index,ymin=0,ymax=upr),alpha=0.2)+
  geom_violin(data=distributions)+#,draw_quantiles = c(0.025,0.975)) +
  geom_point() +
  geom_linerange(aes(ymin=mean.rel-se.rel,ymax=mean.rel+se.rel)) +
  geom_text(data=dat_text, aes(x = -Inf, y = Inf, label = name), size=10, hjust = -0.5, vjust=1.2) +
  facet_wrap(~s.c,nrow=2,scale='free_x') +
  theme_classic() +
  scale_y_continuous(expression(paste("Mean Pairwise Relatedness (",bar(italic(" r")),')'))) +
  scale_x_discrete('Shoal or Cluster identity')+
  theme(strip.text = element_blank(),
        strip.background = element_blank())
ggsave(filename = '../Figures/Figure 4.eps',plot = print(fig4),scale=1.5,height=4,width = 4,device = cairo_ps,dpi = 300)