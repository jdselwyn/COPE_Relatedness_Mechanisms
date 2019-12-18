a_col<-1+seq(1,14,by=2)

library("furrr")
library(fs)

make_dir<-function(x){
  dir_create(x$tmp_dir)
}


tmp<-list.files(path = '../../Results_tmp/', pattern = 'all_combos_prefiltering_microsats_relatedness[0-9]*.csv',full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows() %>%
  group_by(number_shared, combo) %>%
  nest

tmp %>%
  arrange(-number_shared)

the_start <- tibble(number_shared = 1:7) %>% #Something broken with even numbers - no idea what/why c(seq(1,7,by=2), 6)
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
  anti_join(tmp)
  mutate(tmp_dir = str_c(tempdir(),'/',number_shared,'_',combo)) %T>%
  make_dir 

rm('excess.out', 'excess.test', 'marker_check', microsats, shoal_locations, cope)

number_workers<-6
for(i in 1:(nrow(the_start) %/% number_workers+1)){
  plan(multicore, workers = number_workers)
  the_start %>%
    slice( 1:number_workers + number_workers * (i-1)) %>%
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
    unnest %>%
    write_csv(str_c('../../Results_tmp/all_combos_prefiltering_microsats_relatedness',z,'.csv')) 
  z<-z+number_workers
  plan(sequential)
}






the_start <- tibble(number_shared = 1:7) %>% #Something broken with even numbers - no idea what/why c(seq(1,7,by=2), 6)
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
  ungroup



list.files(path = '../../Results_tmp/', pattern = 'all_combos_prefiltering_microsats_relatedness[0-9]*.csv',full.names = TRUE) %>%
  map(read_csv) %>%
  bind_rows() %>%
  write_csv('all_combos_prefiltering_microsats_relatedness.csv')
