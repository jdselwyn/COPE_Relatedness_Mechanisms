setwd('C:/Users/jdsel/Documents/Coryphopterus/Mechanisms to Relatedness (Paper 2 - MS Ch. 2)/Paper/Submission Files/Relatedness_Mechanisms/Code')

pairwise_related<-read_csv('../../Results_tmp/pairwise_relatedness.csv')

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

matchup_data<-pairwise_related %>%
  inner_join(marker_check, by = c('ind1.id'='ID')) %>%
  inner_join(marker_check, by = c('ind2.id'='ID'), suffix = c('.1','.2')) %>%
  mutate(COPE5 = COPE5.1 == COPE5.2,
         COPE9 = COPE9.1 == COPE9.2,
         CPER26 = CPER26.1 == CPER26.2,
         CPER92 = CPER92.1 == CPER92.2,
         CPER99 = CPER99.1 == CPER99.2,
         CPER119 = CPER119.1 == CPER119.2,
         CPER188 = CPER188.1 == CPER188.2) %>%
  dplyr::select(-X1,-pair.no,-Shoal.x:-CPER188.2) %>%
  mutate(number_match = COPE5 + COPE9 + CPER26 + CPER92 + CPER99 + CPER119 + CPER188) %>%
  dplyr::select(-COPE5:-CPER188) %>%
  mutate(dyadml2 = if_else(dyadml.low <= max_unrel, 0, dyadml)) %>%
  dplyr::select(-dyadml.low, -dyadml.high) 


matchup_data %>%
  gather(type, dyad, -ind1.id, -ind2.id, -number_match) %>%
  mutate(type = if_else(type=='dyadml2', 'Filtered to remove any overlap with UR CI', 'All potentially related pairs')) %>%
  filter(dyad > 0) %>%
  
  ggplot(aes(x=number_match, y=dyad, group=interaction(number_match, type), fill = type)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = 'Number of successsfully amplified loci in common', 
       y = 'Pairwise relatedness of pairs identified as related',
       fill = NULL) +
  theme_classic()

a<-1;b<-1
matchup_data2<-matchup_data %>%
  gather(type, dyad, -ind1.id, -ind2.id, -number_match) %>%
  mutate(type = if_else(type=='dyadml2', 'Filtered to remove any overlap with UR CI', 'All potentially related pairs')) %>%
  filter(number_match > 1) %>%
  group_by(type, number_match) %>%
  summarise(total = n(), n = sum(dyad>0),
            prop = n/total,
            PI_25=qbeta(0.25,n+a,total-n+b),
            PI_75=qbeta(0.75,n+a,total-n+b),
            PI_2.5=qbeta(0.025,n+a,total-n+b),
            PI_97.5=qbeta(0.975,n+a,total-n+b)) 

matchup_data2 %>%
  ggplot(aes(x=number_match, y=prop, group=interaction(number_match, type), colour = type)) +
  geom_linerange(aes(ymin = PI_2.5, ymax = PI_97.5), lty=2, position = position_dodge(0.2)) +
  geom_linerange(aes(ymin = PI_25, ymax = PI_75), position = position_dodge(0.2)) +
  geom_point(position = position_dodge(0.2)) +
  labs(x = 'Number of successsfully amplified loci in common', 
       y = 'Proportion of pairs identified as related',
       colour = NULL) +
  theme_classic()


library(INLA)

tmp_data<-matchup_data %>%
  gather(type, dyad, -ind1.id, -ind2.id, -number_match) %>%
  filter(dyad > 0) %>%
  add_row(.before = 1, number_match = 1:7, type = 'dyadml') %>%
  add_row(.before = 1, number_match = 1:7, type = 'dyadml2')
  

rel_common_res <- inla(dyad ~ number_match + 
                          
                          f(type, model="iid") +
                         f(ind1.id, model="iid") +
                         f(ind2.id, model="iid"), 
                
                family = c("gamma"),
                data = tmp_data,
                control.predictor = list(link = 1),
                control.compute = list(dic = TRUE, 
                                       waic = TRUE, 
                                       cpo = TRUE,
                                       config = TRUE),
                control.results = list(return.marginals.predictor = FALSE, 
                                       return.marginals.random = FALSE),
                verbose = TRUE)
summary(rel_common_res)


rel_common_res %$%
  summary.fitted.values %>%
  as_tibble() %>% 
  slice(1:14) %>%
  mutate(number_match = c(1:7, 1:7),
         type = rep(c('dyadml2','dyadml'), each = 7)) %>%
  rename(q_2.5 = `0.025quant`, 
         q_5 = `0.5quant`, 
         q_975 = `0.975quant`) %>%
  
  ggplot(aes(x = number_match, y = mean, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = q_2.5, ymax = q_975), alpha = 0.2) +
  geom_line() +
  geom_boxplot(data = tmp_data, aes(y=dyad, group=interaction(number_match, type)), alpha = 0.4) +
  scale_color_discrete(labels= c('All potentially related pairs','Filtered to remove any overlap with UR CI')) +
  scale_fill_discrete(labels= c('All potentially related pairs','Filtered to remove any overlap with UR CI')) +
  labs(x = 'Number of successsfully amplified loci in common', 
       y = 'Pairwise relatedness of pairs identified as related',
       fill = NULL,
       colour = NULL,
       title = 'Relatedness values when pair is related') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')


#### Binomial ####
tmp_data2<-matchup_data %>%
  gather(type, dyad, -ind1.id, -ind2.id, -number_match) %>%
  mutate(p.a = as.integer(dyad > 0)) %>%
  filter(number_match != 1) %>%
  group_by(number_match, type) %>%
  summarise(success = sum(p.a), total = n()) %>%
  ungroup %>%
  add_row(.before = 1, number_match = 1:7, type = 'dyadml') %>%
  add_row(.before = 1, number_match = 1:7, type = 'dyadml2') %>%
  mutate(type.numb = str_c(type, number_match, sep='_'))

  
  
  
prop_common_res <- inla(success ~ number_match + 
                          
                          f(type, model="iid") +
                          f(type.numb, model='iid'), 
                        
                        family = c("binomial"),
                        Ntrials = total,
                        data = tmp_data2,
                        control.predictor = list(link = 1),
                        control.compute = list(dic = TRUE, 
                                               waic = TRUE, 
                                               cpo = TRUE,
                                               config = TRUE),
                        control.results = list(return.marginals.predictor = FALSE, 
                                               return.marginals.random = FALSE),
                        verbose = TRUE)
summary(prop_common_res)


prop_common_res %$%
  summary.fitted.values %>%
  as_tibble() %>% 
  slice(1:14) %>%
  mutate(number_match = c(1:7, 1:7),
         type = rep(c('dyadml2','dyadml'), each = 7)) %>%
  #filter(number_match <= 7) %>%
  mutate(type = if_else(type=='dyadml2', 'Filtered to remove any overlap with UR CI', 'All potentially related pairs')) %>%
  rename(q_2.5 = `0.025quant`, 
         q_5 = `0.5quant`, 
         q_975 = `0.975quant`) %>%
  
  ggplot(aes(x = number_match, y = mean, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = q_2.5, ymax = q_975), alpha = 0.2) +
  geom_line() +
  geom_point(data = matchup_data2, aes(y=prop, group=interaction(number_match, type)), position = position_dodge(0.1)) +
  geom_linerange(data = matchup_data2, aes(y=prop, group=interaction(number_match, type),
                                           ymin = PI_2.5, ymax = PI_97.5), lty=2, position = position_dodge(0.1)) +
  geom_linerange(data = matchup_data2, aes(y=prop, group=interaction(number_match, type),
                                           ymin = PI_25, ymax = PI_75), position = position_dodge(0.1)) +
  labs(x = 'Number of successsfully amplified loci in common', 
       y = 'Proportion of pairs identified as related',
       colour = NULL,
       fill = NULL,
       title = 'Proportion relatives') +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal')


