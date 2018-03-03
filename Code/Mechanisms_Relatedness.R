rm(list=ls())
setwd('C:/Users/jdsel/Documents/Coryphopterus/Mechanisms to Relatedness (Paper 2 - MS Ch. 2)/Paper/Submission Files/Relatedness_Mechanisms/Data/')

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

#### Functions ####
absmax <- function(x) { x[which.max( abs(x) )]}

#### Microsat Data ####
microsats<- read_csv('../Data/Microsat.csv')  %>% #is this right
  
  #Remove Cloud 4 since it only has 3 individuals
  filter(Cloud!="COPE-M-04") %>% 
  
  #Remove COPE10 and CPER 52 due to too many missing values (106/381) and (149/381)
  select(-starts_with('COPE10')) %>%
  select(-starts_with('CPER52'))

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
  unite(CPER26,CPER26.a:CPER26.b,sep='/') %>%
  unite(CPER92,CPER92.a:CPER92.b,sep='/') %>%
  unite(CPER99,CPER99.a:CPER99.b,sep='/') %>%
  unite(CPER119,CPER119.a:CPER119.b,sep='/') %>%
  unite(CPER188,CPER188.a:CPER188.b,sep='/') %>%

  select(-SL:-Plate) %$% 
  
  #Read into ADEGENET
  df2genind(.[,c(-1,-2)],sep='/',ind.names=ID,pop=Shoal,NA.char =NA,type='codom')

cope.Ho<-lapply(seppop(microsats), function(x) summary(x)$Hobs)
cope.Ho<-matrix(unlist(cope.Ho), ncol = 7, byrow = TRUE)
colnames(cope.Ho)<-levels(microsats@loc.fac)
rownames(cope.Ho)<-levels(microsats@pop)

#### HWE ####
cope.hw<-lapply(seppop(microsats), hw.test,B=10000)
cope.hw<-matrix(p.adjust(unlist(lapply(cope.hw, function(x) x[,'Pr.exact'])),'holm'),ncol=7,byrow=T)
colnames(cope.hw)<-levels(microsats@loc.fac)
rownames(cope.hw)<-levels(microsats@pop)

#### Homozygote Excess ####
excess.test<-hwx.test(microsats,B=10000,statName='U',method='monte')
excess.out<-hwdf(excess.test,statName='U',showN=T,showk=T)
excess.out$`P-val(U)`<-p.adjust(excess.out$`P-val(U)`,'holm')

#### Pairwise Relatedness ####
pairwise_related <- cope %>%
  select(-Shoal:-Plate) %>%
  
  #slice(sample(381,10)) %>% #Shrink for testing code
  
  coancestry(.,error.rates = 0.01,
             trioml=2,allow.inbreeding=T,ci95.num.bootstrap=1000,trioml.num.reference=100) %$% #Run takes ~2.5 hours with full dataset
  inner_join(relatedness,relatedness.ci95,by='pair.no') %>%
  select(pair.no,ind1.id.x,ind2.id.x,trioml,trioml.low,trioml.high) %>%
  rename(ind1.id=ind1.id.x, ind2.id=ind2.id.x) %>%
  
  # Set to unrelated if confidence interval includes 0
  mutate(trioml=if_else(trioml.low==0,0,trioml)) %>%
  filter(trioml!=1) %>%
  
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
  mutate(distance=sqrt((X.x-X.y)^2+(Y.x-Y.y)^2))

write_csv(pairwise_related,'C:/Users/jdsel/Documents/Coryphopterus/Mechanisms to Relatedness (Paper 2 - MS Ch. 2)/Paper/Submission Files/Relatedness_Mechanisms/Data/Relatedness_results.csv')

############################ READ IN FOR TESTING NOT TO GO IN FINAL SCRIPT ##########################################################
pairwise_related<-read_csv('C:/Users/jdsel/Documents/Coryphopterus/Mechanisms to Relatedness (Paper 2 - MS Ch. 2)/Paper/Submission Files/Relatedness_Mechanisms/Data/Relatedness_results.csv')[,-1] %>%
  select(ind1.id,ind2.id,trioml) %>% 
  filter(trioml!=1) %>%
  
  #Add in Individual A information (Shoals, Sizes, X, Y) 
  inner_join(cope,by=c('ind1.id'='ID')) %>%
  select(-COPE5.a:-CPER188.b) %>%
  
  #Add in Individual B information (Shoals, Sizes, X, Y) 
  inner_join(cope,by=c('ind2.id'='ID')) %>%
  select(-COPE5.a:-CPER188.b) %>%
  
  #Determine pairwise variables (same shoal, |size difference|, distance apart)
  mutate(same.shoal=if_else(Shoal.x==Shoal.y,T,F)) %>%
  mutate(same.plate=if_else(Plate.x==Plate.y,T,F)) %>%
  mutate(size.difference=abs(SL.x-SL.y)) %>%
  mutate(distance=sqrt((X.x-X.y)^2+(Y.x-Y.y)^2))

######################################################################################################

#### Make Graph ####
cope.graph <- pairwise_related %>%
  
  #Remove unrelated edges
  filter(trioml!=0) %>%

  #Set relatedness as edge weights
  rename(weight=trioml) %>%
  
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

sizes(cope_community)

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
    select(contains(type),trioml,ind1.id,ind2.id) %>%
    rename_(P1 = names(.)[1], P2=names(.)[2]) %>%
    mutate(same_perm=P1==P2) %>%
    filter(same_perm==T) %>%
    group_by(P1) %>%
    summarise(individuals=length(unique(c(ind1.id,ind2.id))),total_dyads=n(),related.dyads=sum(trioml>0),
              mean.rel=mean(trioml),sd.rel=sd(trioml)) %>%
    mutate(se.rel=sd.rel/sqrt(total_dyads)) %>%
    rename(!!type := P1)
}

full_permutation_related<-function(x,type,NPERM=1000){
  perms<-x %>%
    modelr::permute(NPERM,trioml) %>%
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
  aov(mean.rel ~ s.c, data=.) %>%
  summary()

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
GENERATION<-9 #Generation ~10 mm SL per Cole & Robertson 1988

cope.graph %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  mutate(gens_between=size.difference %/% GENERATION) %>%
  as_tibble() %>%
  group_by(gens_between) %>%
  summarise(n=n(),mean.rel=mean(weight),sd.rel=sd(weight)) %>%
  ungroup() %>%
  mutate(proportion_intergenerational=n/sum(n))


cope.graph %>%
  complementer %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  mutate(size.difference=abs(.N()$SL[from]-.N()$SL[to]))%>%
  mutate(gens_between=size.difference %/% GENERATION) %>%
  as_tibble() %>%
  group_by(gens_between) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(proportion_intergenerational=n/sum(n))
  

#### Figures ####
## Figure 1 ##
#Modified from https://stackoverflow.com/questions/43984614/rggplot2geom-points-how-to-swap-points-with-pie-charts

comm_fills<-c(brewer.pal(12,'Paired'),'black')

comm_fills <- c(distinctColorPalette(sum(sizes(cope_community)>1)),'black')

plot(1:(sum(sizes(cope_community)>1)+1),1:(sum(sizes(cope_community)>1)+1),col=comm_fills,pch=16,cex=5)

g<-shoal_weights %>% 
  as(., "symmetricMatrix") %>% 
  graph.adjacency(., mode="undirected",weighted = T) 
V(g)$Shoal<-LETTERS[1:length(unique(cope$Shoal))]


singles<-paste('C',1:length(cope_community),sep='')[sizes(cope_community)==1]

plot_data<-cope.graph %>%
  as_tbl_graph %>%
  mutate(community= paste('C',(cope_community %>% membership),sep='')) %>% 
  as_tibble() %>%
  mutate(community_simp=ifelse(community %in% singles,'Singleton',community)) %>%
  group_by(Shoal,community_simp) %>%
  summarise(n=n(),x=mean(X),y=mean(Y)) %>%
  ungroup() %>%
  #mutate(community=ifelse(n==1,'Singletons',community)) %>%
  group_by(Shoal) %>%
  mutate(total=sum(n)) %>%
  group_by(Shoal,x,y,total)


names(comm_fills)<-sort(unique(plot_data$community_simp))
df.grobs <- plot_data %>% 
  do(subplots = ggplot(., aes(1, n, fill = community_simp)) + 
       geom_col(position = "fill", alpha = 0.5, colour = NA) + 
       scale_fill_manual(values=comm_fills) + 
       coord_polar(theta = "y") + 
       theme_void() + guides(fill = F))  %>%
  mutate(subgrobs = list(annotation_custom(ggplotGrob(subplots),
                                           xmin = x-100, ymin = y-100,
                                           xmax = x+100, ymax = y+100)))

df.grobs$subplots[[3]]

g %>% 
  as_tbl_graph %>%
  inner_join(df.grobs,by='Shoal') %>%
  ggraph(layout='manual',node.positions=(shoal_XY %>% as_tibble %>% rename(x=X,y=Y))) +
  geom_edge_link() +
  scale_x_continuous('X (cm)',expand = c(0.25, 0)) +
  scale_y_continuous('Y (cm)',expand = c(0.25, 0)) +
  coord_fixed ()+
  theme_classic() +
  df.grobs$subgrobs +
  geom_node_text(aes(label = Shoal)) + 
  geom_col(data = plot_data,
           aes(0,0, fill = community_simp), 
           colour = "white", alpha = 0.5) +
  scale_fill_manual(values=comm_fills)
ggsave(filename = '../Figures/Figure 1.pdf',scale=2)

## Figure 2 ##
cope.layout<-cope.graph %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  mutate(cross=crossing(cope_community,cope.graph)) %>%
  filter(!cross) %>%
  create_layout(layout = 'nicely') %>%
  select(name,x,y)


plot_data<-cope.graph %>%
  as_tbl_graph() %>%
  mutate(community= paste('C',(cope_community %>% membership),sep='')) %>%
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
  group_by(community) %>%
  mutate(hull=1:n(),hull=as.numeric(factor(hull, chull(x, y)))) %>%
  arrange(hull) %>%
  filter(!is.na(hull)) %>%
  mutate(polygon_type=n_distinct(hull)) %>%
  filter(polygon_type>1)
  

plot_data %>% 
  ggraph(layout='manual',node.positions=cope.layout) +
  geom_polygon(data = community_data,aes(x=x,y=y,fill=community,group=community), alpha = 0.5) +
  geom_edge_link(aes(colour=weight)) + 
  geom_node_text(aes(label=Shoal),size=3) +
  facet_wrap(~relationship_type,nrow=1,dir='h') +
  #scale_fill_discrete('Community',na.value=NA) +
  scale_fill_manual(values=comm_fills) + 
  scale_edge_color_continuous('Relatedness',low='gray80',high='gray20') +
  guides(shape=F,fill=F)+
  coord_fixed() +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave(filename = '../Figures/Figure 2.pdf',scale=2)


## Figure 3 ##
plot_data<-group_relatedness %>%
  mutate(label=ifelse(is.na(Shoal),community,Shoal)) %>%
  group_by(s.c) %>%
  mutate(index=1:n()) %>%
  ungroup() %>%
  mutate(s.c=factor(s.c, levels=c('Shoal','Community')))

ggplot(plot_data,aes(x=label,y=mean.rel,ymin=mean.rel-se.rel,ymax=mean.rel+se.rel))+
  geom_blank()+
  geom_ribbon(data=plot_data,aes(x=index,ymin=0,ymax=upr),alpha=0.2)+
  geom_point()+
  geom_linerange() +
  facet_wrap(~s.c,nrow=2,scale='free_x') +
  theme_classic() +
  scale_y_continuous(expression(paste("Mean Pairwise Relatedness (",bar(italic(" r")),')'))) +
  scale_x_discrete('Shoal or Community identity')
ggsave(filename = '../Figures/Figure 3.pdf',scale=2)

## Figure zz ##
cope %>%
  group_by(Shoal,community) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(community=if_else(n==1,'singleton',community)) %>%
  ggplot(aes(x=Shoal,fill=community)) + 
  geom_bar(position='fill')

cope %>% 
  group_by(Shoal,community) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(community=if_else(n==1,'singleton',community)) %>%
  group_by(Shoal,community) %>%
  summarise(n=sum(n)) %>%
  ggplot(aes(x=Shoal,y=community,fill=n)) + 
  geom_bin2d()


## Figure Supplement A ##
plot_data<-cope.graph %>%
  as_tbl_graph() %>%
  mutate(community= paste('C',(cope_community %>% membership),sep='')) %>%
  activate(edges) %>%
  mutate(cross=crossing(cope_community,cope.graph)) %>%
  mutate(gens_between=size.difference %/% GENERATION) %>%
  mutate(relationship_type=if_else(gens_between==0,'Within Cohort','Intergenerational')) %>%
  activate(nodes) %>%
  inner_join(cope.layout,by='name')


community_data <- plot_data %>%
  as_tibble() %>%
  group_by(community) %>%
  mutate(hull=1:n(),hull=as.numeric(factor(hull, chull(x, y)))) %>%
  arrange(hull)


plot_data %>% 
  ggraph(layout='manual',node.positions=cope.layout) +
  #geom_polygon(data = filter(community_data, !is.na(hull)),aes(x=x,y=y,fill=as.factor(community),group=community), alpha = 0.3) +
  geom_edge_link(aes(colour=cross,alpha=weight)) + 
  geom_node_point(aes(shape=Shoal),size=2.5) +
  facet_wrap(~relationship_type,nrow=1,dir='h') +
  scale_shape_manual(values=LETTERS[1:13]) +
  scale_fill_discrete('Community',na.value=NA) +
  #scale_edge_color_discrete('Relatedness',low='gray80',high='gray20') +
  scale_color_brewer('Generation',palette = 'Set1') +
  guides(shape=F,fill=F)+
  coord_fixed() +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
ggsave(filename = '../Figures/Figure S1.pdf',scale=2)
