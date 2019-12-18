pairwise_related<-read_csv('../../Results_tmp/pairwise_relatedness.csv')

cope.graph <- pairwise_related %>%
  
  select(-X1,-pair.no) %>%
  
  #Remove unrelated edges
  #filter(dyadml!=0) %>%
  
  #Set relatedness as edge weights
  rename(weight=dyadml) %>%
  
  #Make Graph
  graph_from_data_frame(d=., vertices=cope, directed=F) %T>%
  
  #Plot to check everything went through right
  plot(.,vertex.label=NA,vertex.size=4)


E(cope.graph)$weight


start<-cope.graph
n.e<-length(E(cope.graph))
prob.e<-E(cope.graph)$weight


E(cope.graph)$keep<-rbinom(length(E(cope.graph)), 1, prob=E(cope.graph)$weight)
test2<-delete_edges(cope.graph, which(test==0))

plot(test2,vertex.label=NA,vertex.size=4)


graph_update<-function(max_edges, edge_probabilities){
  edge_flips<-rbinom(max_edges, 1, prob=edge_probabilities)
  sum(edge_flips)
}

library(furrr)

plan(multiprocess)
tmp<-tibble(sim = 1:100000) %>%
  mutate(edge_count = future_map_int(sim, function(x) graph_update(n.e, prob.e))) 


tmp %>%
  ggplot(aes(x=edge_count)) +
  geom_histogram(bins=30)

tmp %>%
  ggplot(aes(x=sim, y=edge_count)) +
  geom_line()
