library(dplyr)
library(purrr)
library(tidyr)
library(rcompanion)
library(ggplot2)
library(magrittr)
library(dendextend)
summarize_vcf <- function(df){
  df%>%
    separate(col = V5, sep = '%', into = c('freq','thresh'))%>%
    separate(col = freq,sep=':', into = c(map_chr(c(1:6),~paste0(.,'_thresh')),'mistake_percent'))%>%
    select(-contains('thresh'))%>%
    mutate(mistake_percent=as.numeric(mistake_percent))%>%
    summarize(mean_freq=mean(mistake_percent),sd_freq=sd(mistake_percent), low_interval=mean_freq-3*sd_freq,high_interval=mean_freq+3*sd_freq)
} # provide mean, sd, mean-3*sd, mean+3*sd 
control_1 <- read.table('bioinf_mini_project/control_1.tsv',header = F, sep = '\t') # load data
control_2 <- read.table('bioinf_mini_project/control_2.tsv',header = F, sep = '\t') # load data
control_3 <- read.table('bioinf_mini_project/control_3.tsv',header = F, sep = '\t') # load data
roommate <- read.table('bioinf_mini_project/roommate.tsv',header = F, sep = '\t') # load data
list_of_dataframes <- list(control_1,control_2,control_3)
statistic_for_each_control <- map_df(list_of_dataframes,~summarize_vcf(.))
roommate_results <- roommate%>%
  separate(col = V5, sep = '%', into = c('freq','thresh'))%>%
  separate(col = freq,sep=':', into = c(map_chr(c(1:6),~paste0(.,'_thresh')),'mistake_percent'))%>%
  select(-contains('thresh'))%>%
  mutate(mistake_percent=as.numeric(mistake_percent))%>%
  mutate(results= if_else(mistake_percent>max(statistic_for_each_control$high_interval),'Interesting result','Mistake'))%>%
  filter(results!='Mistake') # returns result for 3 sigma approac
names(roommate_results) <- c('Region', 'position', 'ref_allele','new_allele','mistake_percent','results')
# But this approach seems to be too obvious. Making assumptions that the results of the controls were obtained from one sample,
# using one batch of reagents and on the same equipment, we can combine these experiments into one.
full_control <- rbind(control_1,control_2, control_3)
full_control_results <-  full_control%>%
  separate(col = V5, sep = '%', into = c('freq','thresh'))%>%
  separate(col = freq,sep=':', into = c(map_chr(c(1:6),~paste0(.,'_thresh')),'mistake_percent'))%>%
  select(-contains('thresh'))%>%
  mutate(mistake_percent=as.numeric(mistake_percent))
names(full_control_results) <- c('Region', 'position', 'ref_allele','new_allele','mistake_percent')
# We plot the distribution density function. As we see the frequency distribution is abnormal and has a clear bias.
# Which is actually logical, because the probability of an error comes from a beta distribution, not a normal one, and the 3 sigma rule does not work there.
ggplot(full_control_results,aes(x=mistake_percent))+geom_density()+geom_rug()+theme_bw()
# Using bootstrap BCA
real_interval <- groupwisePercentile(mistake_percent~1,data = full_control_results,conf=0.99, tau=0.99, R=5000, bca = T)
roommate_results_method_2 <- roommate%>%
  separate(col = V5, sep = '%', into = c('freq','thresh'))%>%
  separate(col = freq,sep=':', into = c(map_chr(c(1:6),~paste0(.,'_thresh')),'mistake_percent'))%>%
  select(-contains('thresh'))%>%
  mutate(mistake_percent=as.numeric(mistake_percent))%>%
  mutate(results= if_else(mistake_percent>real_interval[1,7],'Interesting result','Mistake'))%>%
  filter(results!='Mistake') # results with new thresh hold
# Cluster analysis
control_cluster<- as.data.frame(scale(full_control_results$mistake_percent))
dist_mat <- dist(control_cluster, method = 'euclidean') # distance matrices
hclust_avg <- hclust(dist_mat, method = 'average') # actual cluster
cut_avg <- cutree(hclust_avg, k = 2) # cluster tree
avg_dend_obj <- as.dendrogram(hclust_avg) # beautify
avg_col_dend <- color_branches(avg_dend_obj,h=2) # beautify
plot(avg_col_dend) # plot tree
cluster_results = mutate(control_cluster,cluster=cut_avg)
full_control_results$mistake_cluster <- cluster_results$cluster
full_control_results %<>%
  mutate(mistake_cluster = factor(mistake_cluster,levels = c(1,2), labels=c('PCR','Sequence')))# add cluster column
ggplot(full_control_results, aes(x=position,y=mistake_percent, color=mistake_cluster))+geom_point()+theme_bw()+scale_y_continuous('Percent of the mistake')+ggtitle('Visualizing mistakes by position')