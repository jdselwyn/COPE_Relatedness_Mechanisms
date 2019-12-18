setwd('C:/Users/jdsel/Documents/Coryphopterus/Mechanisms to Relatedness (Paper 2 - MS Ch. 2)/Paper/Submission Files/Relatedness_Mechanisms/Code')

#### Code to run relatedness with all combinations of loci when individuals both have all relevant loci
a_col<-1+seq(1,14,by=2)

library("furrr")
library(fs)

make_dir<-function(x){
  dir_create(x$tmp_dir)
}




cope %>%
  inner_join(marker_check %>%
               mutate(number_loci = select(., starts_with('COPE'), starts_with('CPER')) %>% apply(1, sum, na.rm=TRUE)) %>%
               dplyr::select(ID, number_loci),
             by = 'ID') %>%
  filter(number_loci == 7) %>%
  dplyr::select(-Shoal:-Plate,-number_loci) %>%
  dplyr::select(ID, sort(c(y, 1+y)))

#seq(1,7,by=2)
#seq(4,6,by=2)
plan(multicore, workers = 20) #breaks running in parallel 
test<-tibble(number_shared = 1:7) %>% #Something broken with even numbers - no idea what/why c(seq(1,7,by=2), 6)
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

test <- bind_rows(test_1.3, test_4.5, test_6.7) %T>%
  write_csv('../../Results_tmp/all_combos_prefiltering_microsats_relatedness.csv')

setwd('/home/jason/Documents/Coryphopterus/Mechanisms to Relatedness (Paper 2 - MS Ch. 2)/Paper/Submission Files/Relatedness_Mechanisms/Code')
test<-read_csv('../../Results_tmp/all_combos_prefiltering_microsats_relatedness.csv')





shared_microsat_relatedness %>%
  mutate(dyadml=if_else(dyadml.low==0,0,dyadml)) %>%
  filter(dyadml != 1) %>%
  filter(ind1.id == '01COPE-0001') %>%
  filter(ind2.id == '01COPE-0002' | ind2.id == '01COPE-0003' | ind2.id == '01COPE-0004') %>%
  mutate(pair = str_c(ind1.id, ind2.id, sep = '_')) %>%
  
  group_by(pair, number_shared) %>%
  summarise(n = n(), relatedness = mean(dyadml), sd.relatedness = sd(dyadml)) %>%
  mutate(se.relatedness = sd.relatedness/sqrt(n), lwr = relatedness - se.relatedness,
         upr = relatedness + se.relatedness) %>%
  
  ggplot(aes(x = number_shared, y = relatedness, colour = pair, ymin = lwr, ymax = upr)) +
  geom_point() +
  geom_linerange()



## Model predicting relatedness based on number of amplified loci
library(INLA)

tmp <- shared_microsat_relatedness %>%
  filter(dyadml != 1) %>%
  mutate(dyadml=if_else(dyadml.low==0,0,dyadml)) %>% #Filter 1
  #mutate(dyadml=if_else(dyadml < max_unrel, 0 ,dyadml)) %>% #Filter 2
  #mutate(dyadml=if_else(dyadml.low < max_unrel, 0, dyadml)) %>% #Filter 3
  mutate(pair = str_c(ind1.id, ind2.id, sep = '_')) %>%
  mutate(added_bit = rgamma(n = nrow(.), shape = 1, rate = 1000),
         testing_relatedness = dyadml + added_bit) %>%
  group_by(pair) %>%
  slice(sample(nrow(.), 0.1*nrow(.))) %>%
  unnest %>%
  add_row(number_shared = 1:7, .before = 0)


#https://groups.google.com/forum/#!topic/r-inla-discussion-group/iQELaQF8M9Q
rel_common_res <- inla(testing_relatedness ~ f(number_shared, model = 'rw2')
                       + f(pair, model = 'iid'),
                       
                       family = c("beta"), #change to beta for dyadml
                       data = tmp,
                       control.family = list(link='logit'),
                       control.predictor = list(link = 1, compute = TRUE),
                       control.compute = list(dic = TRUE, 
                                              waic = TRUE, 
                                              cpo = TRUE),
                       control.inla = list(strategy = 'adaptive', 
                                           int.strategy = 'eb',
                                           cmin = 0),
                       verbose = TRUE)
summary(rel_common_res)

rel_common_res %$%
  summary.fitted.values %>%
  as_tibble() %>%
  bind_cols(tmp,.) %>%
  filter(is.na(pair)) %>%
  rename(q0.025 = `0.025quant`, q0.975 = `0.975quant`) %>%
  
  ggplot(aes(x = number_shared, y = mean, ymin = q0.025, ymax = q0.975)) +
  geom_line(lty = 2) +
  geom_point() +
  geom_linerange()


rel_common_res %$%
  summary.fitted.values %>%
  as_tibble() %>%
  bind_cols(tmp,.) %>%
  filter(!is.na(pair)) %>%
  rename(q0.025 = `0.025quant`, q0.975 = `0.975quant`) %>%
  mutate(error = mean - dyadml) %>%
  mutate(related = dyadml > 0) %>%
  #slice(sample(nrow(.),1000)) %>%
  
  ggplot(aes(x = dyadml, y = error, colour = related)) +
  geom_point()


rel_common_res %$%
  summary.fitted.values %>%
  as_tibble() %>%
  bind_cols(tmp,.) %>%
  filter(!is.na(pair)) %>%
  rename(q0.025 = `0.025quant`, q0.975 = `0.975quant`) %>%
  mutate(error = mean - dyadml) %>%
  mutate(related = dyadml > 0) %>%
  #slice(sample(nrow(.),1000)) %>%
  
  ggplot(aes(x = number_shared, y = error, colour = related)) +
  geom_point()
