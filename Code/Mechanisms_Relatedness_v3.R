#### Libraries ####
rm(list=ls())

library(tidyverse)
library(readxl)
library(janitor)
library(magrittr)
library(furrr)
library(rentrez)
library(taxize)
library(msa)
library(seqinr)
library(ape)
library(pegas)
library(igraph)
library(tidygraph)
library(ggraph)
library(RColorBrewer)
library(adegenet)
library(HWxtest)
library(broom)
library(lubridate)
library(scatterpie)
library(related)
library(cowplot)

#### Functions ####
plan('multiprocess')
source('~/R/R Functions/haplotype_network.R')

desaturate <- function(cols, ds=0.4, dv=0.7) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}

absmax <- function(x) { x[which.max( abs(x) )]}

CovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov))) %>% as_tibble()
}

#### Read in Shoal Locations ####
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
  dplyr::rename(cloud = Cloud) %>%
  mutate(site = str_replace(site, '-',' ')) %T>%
  
  #Check that the conversion to XY worked visually 
  plot(Y~X,data=.,pch=16, col = as.factor(site))


#### Read in COPE data ####
cope_data <- unlist(str_split('cdcccccccnnnciiiccclcccccicccc','')) %>%
  tibble(col_type = .) %>%
  mutate(col_type = case_when(col_type == 'c' ~ 'text',
                              col_type == 'd' ~ 'numeric',
                              col_type == 'n' ~ 'numeric',
                              col_type == 'i' ~ 'numeric',
                              col_type == 'l' ~ 'text')) %$%
  # read_excel('~/Coryphopterus/Master COPE collection database.xlsx',sheet = 1, 
  #            col_types = col_type, na = c('N/A','NA')) %>%
  read_excel('~/../Documents/Coryphopterus/Master COPE collection database.xlsx',sheet = 1, 
             col_types = col_type, na = c('N/A','NA')) %>%
  clean_names %>%
  dplyr::select(tube_label, site, cloud, sl_mm) %>%
  dplyr::rename(ID = tube_label) %>%
  mutate(cloud = str_replace(cloud, '-0','-'))

#### Read in post-Sequencher fasta files ####
get_consensus <- function(x){
  if(nrow(x) == 1){
    x$sequence
  } else {
    x$sequence %>%
      str_replace_all("-", ' ') %>%
      str_trim() %>%
      str_replace_all(" ",'-') %>%
      msa(type='dna', verbose = TRUE, order = 'input') %>%
      consensusString()
  }
}

post_sequencher_fasta <- read.fasta('../Data/Coryphopterus_07.12.19.TXT', as.string = TRUE) %>%
  unlist %>%
  enframe() %>%
  dplyr::rename(sequence = value) %>%
  filter(name != 'COHY', name != 'COPE', name != 'Cpersonatus_COI', name != 'Chylalinus_COI') %>% 
  mutate(name = str_replace(name, 'Cope_93_@6/14/2019,_5_00_PM', 'COPE_0093')) %>%
  mutate(name = str_replace(name, 'COPE_', 'COPE-')) %>%
  separate(name, into = c('ID'), sep = '_', extra = 'drop') %>% 
  mutate(ID = str_replace(ID, 'COPE-', 'COPE_')) %>%
  separate(ID, into = c('extra', 'ID'), sep = '-', fill = 'left') %>%
  dplyr::select(-extra) %>%
  mutate(ID = str_replace(ID, 'COPE_', 'COPE-')) %>%
  mutate(sequence = str_to_upper(sequence)) %>%
  group_by(ID, sequence) %>%
  summarise()

post_sequencher_fasta_v2 <- read.fasta('../Data/08.05.19_newseq.TXT', as.string = TRUE) %>%
  unlist %>%
  enframe() %>%
  dplyr::rename(sequence = value) %>%
  filter(name != 'COHY', name != 'COPE', name != 'Cpersonatus_COI', name != 'Chylalinus_COI') %>% 
  mutate(name = str_replace(name, 'Cope_93_@6/14/2019,_5_00_PM', 'COPE_0093')) %>%
  mutate(name = str_replace(name, 'COPE_', 'COPE-'),
         name = str_replace(name, '-Fish', '_Fish')) %>%
  separate(name, into = c('ID'), sep = '_', extra = 'drop') %>% 
  mutate(ID = str_replace(ID, 'COPE-', 'COPE_')) %>%
  separate(ID, into = c('extra', 'ID'), sep = '-', fill = 'left') %>%
  dplyr::select(-extra) %>%
  mutate(ID = str_replace(ID, 'COPE_', 'COPE-')) %>%
  mutate(sequence = str_to_upper(sequence)) %>%
  group_by(ID, sequence) %>%
  summarise()

post_sequencher_fasta_v3 <- read.fasta('../Data/Gobies_09.11.19.TXT', as.string = TRUE) %>%
  unlist %>%
  enframe() %>%
  dplyr::rename(sequence = value) %>%
  filter(name != 'COHY', name != 'COPE', name != 'Cpersonatus_COI', name != 'Chylalinus_COI') %>%
  mutate(name = str_replace(name, 'COPE_', 'COPE-'),
         name = str_replace(name, '-Fish', '_Fish')) %>% 
  separate(name, into = c('ID'), sep = '_', extra = 'drop') %>% 
  mutate(ID = str_replace(ID, 'COPE-', 'COPE_')) %>%
  separate(ID, into = c('extra', 'ID'), sep = '-', fill = 'left') %>%
  dplyr::select(-extra) %>%
  mutate(ID = str_replace(ID, 'COPE_', 'COPE-')) %>%
  mutate(sequence = str_to_upper(sequence)) %>%
  group_by(ID, sequence) %>%
  summarise()

post_sequencher_fasta_v4 <- read.fasta('../Data/11.06.19.TXT', as.string = TRUE) %>%
  unlist %>%
  enframe() %>%
  dplyr::rename(sequence = value) %>%
  filter(name != 'COHY', name != 'COPE', name != 'Cpersonatus_COI', name != 'Chylalinus_COI') %>%
  mutate(name = str_replace(name, 'COPE_', 'COPE-'),
         name = str_replace(name, '-Fish', '_Fish')) %>% 
  separate(name, into = c('ID'), sep = '_', extra = 'drop') %>% 
  mutate(ID = str_replace(ID, 'COPE-', 'COPE_')) %>%
  separate(ID, into = c('extra', 'ID'), sep = '-', fill = 'left') %>%
  dplyr::select(-extra) %>%
  mutate(ID = str_replace(ID, 'COPE_', 'COPE-')) %>%
  mutate(sequence = str_to_upper(sequence)) %>%
  group_by(ID, sequence) %>%
  summarise()


post_sequencher_fasta <- bind_rows(post_sequencher_fasta, post_sequencher_fasta_v2, 
                                   post_sequencher_fasta_v3, post_sequencher_fasta_v4) %>%
  mutate(length = str_length(sequence),
         number_dash = str_count(sequence, '-')) %>%
  distinct %>%
  nest(data = c(sequence, length, number_dash)) %>%
  mutate(n_seq = map_int(data, nrow)) %>%
  mutate(sequence = map_chr(data, get_consensus)) %>%
  select(ID, sequence)


#### Read in microsat data ####
microsats_raw <- read_csv('../Data/Microsat.csv') %>% 
  #dplyr::select(-COPE10.a, -COPE10.b, -CPER52.a, -CPER52.b, -Cloud) %>%
  dplyr::select(-Cloud) %>%
    
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
  dplyr::rename(ID = Sample) %>%
  mutate(ID = str_sub(ID, 3))

#### Join full dataset ####
full_cope_data <- microsats_raw %>%
  left_join(post_sequencher_fasta) %>%
  left_join(cope_data) %>%
  filter(cloud!="COPE-M-4") %>% 
  
  left_join(shoal_locations) %>%
  mutate(shoal = LETTERS[as.integer(as.factor(cloud))]) %>%
  dplyr::rename(shoal_size = number) %>%
  dplyr::select(ID, site, shoal, shoal_size, X, Y, sl_mm, sequence, COPE5:CPER188) %>%
  
  #Remove if I don't have COI sequence data
  filter(!is.na(sequence))
  

full_cope_data %>%
  group_by(ID) %>%
  filter(n() > 1)

#### Read in NCBI COPE/COHY ####
all_cory_search <- '(Coryphopterus personatus[Organism] OR Coryphopterus hyalinus[Organism] OR Coryphopterus lipernes) AND COI[Gene]'
just_baldwin_paper_search <- '(Coryphopterus personatus[Organism] OR Coryphopterus hyalinus[Organism]) AND COI[Gene] AND Reconciling Genetic Lineages with Species in Western Atlantic Coryphopterus (Teleostei: Gobiidae)'

NCBI_sequence_data<-entrez_search(db="nucleotide",term = just_baldwin_paper_search,retmax=9999)$ids %>%
  entrez_fetch(db = 'nucleotide', id = ., rettype = 'fasta') %>%
  str_split(pattern = '>') %>% 
  unlist %>%
  str_subset('mitochondrial') %>%
  str_remove_all('\n') %>%
  str_remove('mitochondrial') %>%
  tibble(sequence = .) %>%
  separate(sequence, into = c('info', 'sequence'), sep = ';') %>%
  mutate(sequence = str_trim(sequence)) %>%
  separate(info, into = c('genbankID_nuc','info'), sep = ' ', extra = 'merge') %>%
  mutate(species = case_when(
    str_detect(info, 'personatus') ~ 'Cpers',
    str_detect(info, 'lipernes') ~ 'Clip',
    str_detect(info, 'hyalinus') ~ 'Chya'
  )) %>%
  filter(str_detect(species, 'Clip', negate = TRUE)) %>%
  filter(str_detect(info, 'UNVERIFIED', negate = TRUE)) %>%
  group_by(species) %>%
  ungroup() %>%
  mutate(node_label = str_c(species,1:n(), sep='_')) %>%
  dplyr::select(node_label, sequence) %>%
  dplyr::rename(ID = 'node_label')

#### Align my fasta's with NCBI ####
aligned_sequence <- NCBI_sequence_data %>%
  bind_rows(full_cope_data, .id = 'source') %>%
  mutate(source = case_when(source == 1 ~ 'NCBI',
                            TRUE ~ 'Collection')) %>%
  filter(!is.na(sequence)) %>%
  dplyr::select(-site:-CPER188) %>%
  nest %>%
  mutate(aligned = map(data, ~.x %$%
                         sequence %>%
                         msa(type='dna', verbose = TRUE, order = 'input') %>%
                         as.character)) %>%
  unnest %>%
  dplyr::select(-sequence)

#### Find diagnostic SNPs ####
## Based on only the NCBI fastas
per_position_base <- aligned_sequence %>%
  mutate(bases = str_split(aligned, '')) %>%
  unnest(bases) %>%
  dplyr::select(-aligned) %>%
  group_by(ID) %>%
  mutate(position = 1:n()) %>%
  ungroup 
  
diagnostic_snps <- per_position_base %>%
  filter(source == 'NCBI') %>%
  dplyr::select(ID, bases, position) %>%
  separate(ID, into = c('species'), sep = '_', remove = FALSE, extra = 'drop') %>%
  group_by(species, position) %>%
  mutate(number_bases = n_distinct(bases)) %>%
  filter(number_bases == 1) %>%
  ungroup %>%
  group_by(position) %>%
  mutate(number_bases = n_distinct(bases)) %>%
  filter(number_bases > 1) %>%
  group_by(species, position, bases) %>%
  summarise %>%
  
  ungroup %>%
  filter(bases != '-', bases != 'N') %>%
  group_by(position) %>%
  filter(n() == 2) %>%
  ungroup %>%
  spread(species, bases)

#### Use diagnostic SNPs to ID fish ####
## Prior is entirely weighted to either 0 (COHY) or 1 (COPE)
a<-0.5; b<-0.5 #Jeffreys Prior - expect to be either one or the other

identification_probability <- per_position_base %>%
  filter(source != 'NCBI') %>%
  inner_join(diagnostic_snps, .) %>%
  filter(bases != '-') %>%
  dplyr::select(-source) %>%
  
  mutate(match = case_when(bases == Chya ~ 'chya',
                           bases == Cpers ~ 'cpers',
                           TRUE ~ 'neither')) %>%
  
  group_by(ID) %>%
  summarise(n_chya = sum(match == 'chya'),
         n_cpers = sum(match == 'cpers'),
         total = n()) %>%
  gather(co1_species, number, -ID, -total) %>%
  mutate(co1_species = str_remove(co1_species, 'n_')) %>%
  
  mutate(prob_id = (number+a)/(total+a+b),
         id_lwr.95 = qbeta(0.025, number+a, total-number+b),
         id_lwr.50 = qbeta(0.25, number+a, total-number+b),
         id_upr.50 = qbeta(0.75, number+a, total-number+b),
         id_upr.95 = qbeta(0.975, number+a, total-number+b)) %>%
  dplyr::select(-number, -total) %>%
  group_by(ID) %>%
  filter(prob_id == max(prob_id)) %>%
  mutate(co1_species = case_when(id_lwr.95 > 0.5 ~ co1_species, 
                                 TRUE ~ 'unknown'))

#### Join CO1 ID to full data ####
full_cope_data <- full_cope_data %>%
  left_join(identification_probability %>%
              dplyr::select(ID, co1_species)) %>%
  mutate(co1_species = if_else(is.na(co1_species), 'unknown', co1_species)) %>%
  filter(co1_species == 'chya')

#### Observed Heterozygosity ####
microsats_species_split <- full_cope_data %>%
  select(-site:-sequence) %$% 
  
  #Read into ADEGENET
  df2genind(.[,c(-1,-11)],sep='/',ind.names=ID,pop=co1_species,NA.char =NA,type='codom')

#### Homozygote Excess ####
summary(microsats_species_split)$loc.n.all
summary(microsats_species_split)$Hobs
summary(microsats_species_split)$Hexp

excess.test<-hwx.test(microsats_species_split,B=10000,statName='U',method='monte')
excess.out<-hwdf(excess.test,statName='U',showN=T,showk=T)

#### Missing Data ####
full_cope_data %>%
  summarise_at(vars(COPE5:CPER188), ~sum(str_detect(., 'NA/NA'))/length(.))


full_cope_data <- full_cope_data %>%
  select(-contains('CPER52'),-contains('COPE10'))

#### Simulate relatives to determine best relatedness measure to use ####
sim_rel <- full_cope_data %>%
  dplyr::select(ID, COPE5:CPER188) %>%
  separate(COPE5, into = c('COPE5.a', 'COPE5.b')) %>%
  separate(COPE9, into = c('COPE9.a', 'COPE9.b')) %>%
  separate(CPER26, into = c('CPER26.a', 'CPER26.b')) %>%
  separate(CPER92, into = c('CPER92.a', 'CPER92.b')) %>%
  separate(CPER99, into = c('CPER99.a', 'CPER99.b')) %>%
  separate(CPER119, into = c('CPER119.a', 'CPER119.b')) %>%
  separate(CPER188, into = c('CPER188.a', 'CPER188.b')) %>%
  readgenotypedata(.) %$%
  familysim(freqs,500) %>%
  coancestry(.,error.rates = 0.01, allow.inbreeding=T,
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
write_csv(sim_rel, '../../Results_tmp/simulated_relatedness.csv')
#sim_rel <- read_csv('../../Results_tmp/simulated_relatedness.csv')

sim_rel %>%
  group_by(metric,group) %>% 
  mutate(M=mean(relatedness),S=sd(relatedness),CV=S/M) %>%
  ungroup %>%
  group_by(metric) %>%
  mutate(mean_cv=mean(CV),w=cor(true_r,relatedness)) %>%
  select(metric,mean_cv,w,group,M,S,CV) %>%
  distinct() %>%
  arrange(-w) %T>%
  write_csv('../Figures/Table 2.csv')

## Need to put back the analysis here to figure out best choice ##

max_unrel <- sim_rel %>%
  filter(group == 'Unrelated',
         metric == 'Milligan 2003') %$%
  quantile(relatedness, 0.975)

#### Probability of identifying a true unrelated pair as unrelated ####
## Adapted from Wang 2007
N<-3
plan(multicore)

probability_exclustion <- full_cope_data %>%
  dplyr::select(ID, COPE5:CPER188) %$%
  df2genind(.[,c(-1)],sep='/',ind.names=ID, NA.char = NA,type='codom') %>%
  
  genind2genpop %>%
  makefreq %>%
  t %>%
  as_tibble(rownames = 'loci.allele') %>%
  dplyr::rename(frequency = `1`) %>%
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

probability_exclustion %$% #%>% filter(loci != 'COPE5')
  Se %>%
  subtract(1,.) %>%
  prod %>%
  subtract(1,.)
#0.9989585

#### KinInfor data prep ####
full_cope_data %>% 
  dplyr::select(ID, COPE5:CPER188) %$%
  df2genind(.[,c(-1)],sep='/',ind.names=ID, NA.char = NA,type='codom') %>%
  
  genind2genpop %>%
  makefreq %>% 
  t %>%
  as_tibble(rownames = 'loci.allele') %>%
  dplyr::rename(frequency = `1`) %>%
  separate(loci.allele, into = c('loci','allele')) %>%
  mutate(allele = as.integer(allele)) %>%
  group_by(loci) %>%
  #summarise(n = n())
  arrange(loci, allele) %>%
  filter(loci == 'CPER188') %$% #Step through loci and copy them into the 'KinInfor.dat' file
  frequency %>%
  str_c(collapse = ' ')

#### Pairwise Relatedness ####
pairwise_related <- full_cope_data %>% 
  dplyr::select(ID, COPE5:CPER188) %>%
  
  #slice(sample(nrow(.),5)) %>% #Shrink for testing code
  
  separate(COPE5, into = c('COPE5.a', 'COPE5.b')) %>%
  separate(COPE9, into = c('COPE9.a', 'COPE9.b')) %>%
  separate(CPER26, into = c('CPER26.a', 'CPER26.b')) %>%
  separate(CPER92, into = c('CPER92.a', 'CPER92.b')) %>%
  separate(CPER99, into = c('CPER99.a', 'CPER99.b')) %>%
  separate(CPER119, into = c('CPER119.a', 'CPER119.b')) %>%
  separate(CPER188, into = c('CPER188.a', 'CPER188.b')) %>%
  
  coancestry(.,error.rates = 0.01,
             dyadml=2,
             allow.inbreeding=T, ci95.num.bootstrap=1000,
             rng.seed = 12345)  %$% #Run takes ~2.5 hours with full dataset
  inner_join(relatedness,relatedness.ci95,by='pair.no') %>%
  select(pair.no,ind1.id.x,ind2.id.x,
         dyadml,dyadml.low,dyadml.high) %>%
  rename(ind1.id=ind1.id.x, ind2.id=ind2.id.x) %>%
  
  # Set to unrelated if confidence interval includes 0
  mutate(dyadml=if_else(dyadml.low==0,0,dyadml)) %>%
  
  #Add in Individual A information (Shoals, Sizes, X, Y) 
  inner_join(full_cope_data,by=c('ind1.id'='ID')) %>%
  
  #Add in Individual B information (Shoals, Sizes, X, Y) 
  inner_join(full_cope_data,by=c('ind2.id'='ID')) %>%
  
  #Determine pairwise variables (same shoal, |size difference|, distance apart)
  mutate(same.shoal=shoal.x==shoal.y) %>%
  mutate(size.difference=abs(sl_mm.x-sl_mm.y)) %>%
  mutate(distance=sqrt((X.x-X.y)^2+(Y.x-Y.y)^2)) %>%
  filter(dyadml!=1) %>%
  
  mutate(COPE5 = !(str_detect(COPE5.x, 'NA') | str_detect(COPE5.y, 'NA')),
         COPE9 = !(str_detect(COPE9.x, 'NA') | str_detect(COPE9.y, 'NA')),
         CPER26 = !(str_detect(CPER26.x, 'NA') | str_detect(CPER26.y, 'NA')),
         CPER92 = !(str_detect(CPER92.x, 'NA') | str_detect(CPER92.y, 'NA')),
         CPER99 = !(str_detect(CPER99.x, 'NA') | str_detect(CPER99.y, 'NA')),
         CPER119 = !(str_detect(CPER119.x, 'NA') | str_detect(CPER119.y, 'NA')),
         CPER188 = !(str_detect(CPER188.x, 'NA') | str_detect(CPER188.y, 'NA'))) %>%
  mutate(number_match = COPE5 + COPE9  + CPER26 + CPER92 + CPER99 + CPER119 + CPER188)
write_csv(pairwise_related, '../../Results_tmp/pairwise_relatedness.csv')
#pairwise_related <- read_csv('../../Results_tmp/pairwise_relatedness.csv')
  
pairwise_related %>%
  filter(dyadml.low > 0)


#### Marker power relatedness analysis ####


#### Make Graph ####
pairwise_filter<-function(x, y){
  #x: loci options
  #y: relatedness option
  pairwise_related %>%
    
    select(-pair.no) %>%
    
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
    graph_from_data_frame(d=., vertices=full_cope_data %>% filter(microsat_species == 'chya'), directed=F) %T>%
    
    #Plot to check everything went through right
    plot(.,vertex.label=NA,vertex.size=4)
  
}

#Loci Options
#All, only_7

#Relatedness options
#all, mean_in_UR, overlap_UR

BOOT<-10
a<-1;b<-1
plan('multiprocess')

cope_graphs<-tibble(loci_option = c('all', 'only_7', NA), related_option = c('no_0','no_mean_UR','no_overlap_UR')) %>%
  tidyr::expand(related_option, loci_option) %>%
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
                           map(.,set_vertex_attr,name='shoal',value=V(x)$shoal))) %>%
  
  #Compare observed to similar random graphs
  mutate(obs_transitivity = map_dbl(graphs, transitivity),
         obs_assortativity = map_dbl(graphs, function(x) assortativity(x,as.factor(V(x)$shoal)))) %>%
  mutate(ER_transitivity = map(ER_random, function(y) map_dbl(y, transitivity)),
         ER_assortativity = map2(graphs, ER_random, function(x, y) map_dbl(y,assortativity,types1=as.factor(V(x)$shoal)))) %>%
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
    rename_(P1 = names(.)[1], P2=names(.)[3]) %>%
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
  mutate(shoal_related = future_map(pairwise, full_permutation_related, 'shoal', NPERM = 1000)) %>%
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
                                                          mutate(size.difference=abs(.N()$sl_mm[from]-.N()$sl_mm[to])) %>%
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
                         dplyr::select(shoal, X, Y, sl_mm) %>%
                         activate(edges) %>%
                         dplyr::select(-dyadml.low:-CPER188) %>%
                         activate(nodes) %>%
                         mutate(shoal = as.factor(shoal)) %>%
                         contract(., mapping = V(.)$shoal,
                                  vertex.attr.comb = list(X = 'mean',
                                                          Y = 'mean',
                                                          shoal = function(x)x[1], "ignore",
                                                          sl_mm = function(x) length(x))) %>%
                         simplify(., remove.multiple = TRUE, remove.loops = FALSE,
                                  edge.attr.comb = list(weight = "sum", 
                                                        distance = 'mean',
                                                        number_match = function(x) length(x))) %>% 
                         as_tbl_graph %>%
                         rename(x=X, y=Y, N = sl_mm) %>%
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
                          dplyr::select(name, shoal, X, Y, sl_mm) %>%
                          inner_join(cope.layout %>% as_tibble %>% mutate(name = as.character(name)),
                                     by = 'name') %>%
                          activate(edges) %>%
                          mutate(gens_between=size.difference %/% GENERATION) %>%
                          mutate(relationship_type=if_else(gens_between==0,'Within Cohort','Intergenerational')) %>%
                          mutate(relationship_type = factor(relationship_type,levels=c('Within Cohort','Intergenerational'))) %>%
                          
                          {if(A == 'only_7') filter(., number_match == 7) else(.)} %>%
                          activate(nodes) %>%
                          
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
  select(loci_option, related_option, shoal, full) %>%
  unnest() %>%
  rename(mean.rel=full)


shoal_relatedness %>%
  mutate(sig = if_else(p < 0.05, 'sig', 'not')) %>%
  
  ggplot(aes(x = shoal, y = mean.rel)) +
  geom_violin(data=distributions) +
  geom_linerange(aes(ymin=mean.rel-se.rel,ymax=mean.rel+se.rel, colour = sig)) +
  geom_point(aes(colour = sig)) +
  facet_grid(related_option~loci_option) +
  theme_classic() +
  scale_y_continuous(expression(paste("Mean Pairwise Relatedness (",bar(italic(" r")),')'))) +
  scale_x_discrete('Shoal or Cluster identity')
