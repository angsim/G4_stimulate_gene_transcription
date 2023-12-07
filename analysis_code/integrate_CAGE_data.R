# import all those data. 


# all_cage_prom
# k562_aggreg
# hek293_cage_clusters.csv
# k562_individ.csv
# k562_individE1.csv
# k562_individE2.csv
# k562_individE3.csv
# 
# all_cage_prom.hek293.centers.hg19.bed
# all_cage_prom.K562.centers.hg19.bed
# 
# relative_dist.ENCFF736NYC_ctcf.all_cage_prom.bed
# relative_dist.k562_G4.all_cage_prom.bed
# relative_dist.hek293_G4.all_cage_prom.bed
# 
# CCAATbox_200_fimo.temp
# GCbox_200_fimo.temp
# TATAbox_50_fimo.temp
# initiator_200_fimo.temp
# 
# slop50L.CpG.bed
# slop200L.CpG.bed

## definition of custom function to import data
import_cage_clusters_individual <- function(file, label_to_use){
  cage_clusters <- read.delim(file,sep = ",", stringsAsFactors = F, header=T)
  colnames(cage_clusters)<- paste0(label_to_use,'_',colnames(cage_clusters))
  cage_clusters$cage_id <- paste0(cage_clusters[,2],'_',cage_clusters[,3],'_',cage_clusters[,4])
  head(cage_clusters)
  
  cage_clusters$cell_exp = rep(label_to_use,length(cage_clusters$cage_id))
  names(cage_clusters)[names(cage_clusters) == "cell_exp"] <- paste0("cell_exp_",label_to_use)
  return(cage_clusters)
}


library(tidyverse)
setwd("/Users/simeon01/Documents/Yuqi/Analasis_distance_Cage_G4")
all_cage_prom <- read.delim('all_cage_prom.bed',sep = "\t", stringsAsFactors = F, header=F)
head(all_cage_prom)
colnames(all_cage_prom) <- c('chr','start','end','feat','cage_id',"strand")


#import hek293
hek293_cage_clusters_fin <- import_cage_clusters_individual('hek293_cage_clusters.csv','hek293')
#import k562 -0
k562_cage_clusters_rep1_fin <- import_cage_clusters_individual('k562_individE1.csv','k562_rep1')
#import k562 -1
k562_cage_clusters_rep2_fin <- import_cage_clusters_individual('k562_individE2.csv','k562_rep2')
#import k562 -2
k562_cage_clusters_rep3_fin <- import_cage_clusters_individual('k562_individE3.csv','k562_rep3')
#import k562 -3
k562_cage_clusters_rep0_fin<- import_cage_clusters_individual('hek293_cage_clusters.csv','k562_rep0')
#import k562  - combined
k562_cage_clusters_comb_fin <- import_cage_clusters_individual('k562_aggreg_consensus.csv','k562_comb')

head(k562_cage_clusters_comb_fin)

temp_join_M <- left_join(all_cage_prom,hek293_cage_clusters_fin, by = "cage_id") %>% left_join(.,k562_cage_clusters_rep1_fin,by = "cage_id") %>% left_join(.,k562_cage_clusters_rep2_fin,by = "cage_id") %>% left_join(.,k562_cage_clusters_rep3_fin,by = "cage_id")%>% left_join(.,k562_cage_clusters_rep0_fin,by = "cage_id") %>% left_join(.,k562_cage_clusters_comb_fin,by = "cage_id")


temp_join_M %>% select(k562_rep1_width) %>% drop_na() %>% ggplot(aes(k562_rep1_width)) +geom_histogram(color='red', fill = 'grey')+ xlim(c(0,100))
ggsave('K562_TC_width.pdf')
temp_join_M %>% select(hek293_width) %>% drop_na() %>% ggplot(aes(hek293_width)) +geom_histogram(color='red', fill = 'grey') + xlim(c(0,100))
ggsave('hek293_TC_width.pdf')
# import BG4- peak centers
# all_cage_prom.hek293.centers.hg19.bed
# all_cage_prom.K562.centers.hg19.bed
hek293_centers <- read.delim("all_cage_prom.hek293.centers.hg19.bed",sep = "\t", stringsAsFactors = F, header=F)
head(hek293_centers)
hek293_centers_sel <- hek293_centers[which(hek293_centers$V10 > -5000 & hek293_centers$V10 <  5000),]
hek293_centers_sel_df <- data.frame(bg4_ID_hek293 = paste0(hek293_centers_sel$V1,"_",hek293_centers_sel$V2,"_",hek293_centers_sel$V3),dist_hek293G4_to_cage=hek293_centers_sel$V10,cage_id=hek293_centers_sel$V8)

k562_centers <- read.delim("all_cage_prom.K562.centers.hg19.bed",sep = "\t", stringsAsFactors = F, header=F)
head(k562_centers)
k562_centers_sel <- k562_centers[which(k562_centers$V10 > -5000 & k562_centers$V10 <  5000),]
head(k562_centers_sel)
k562_centers_sel_df <- data.frame(bg4_ID_k562 = paste0(k562_centers_sel$V1,"_",k562_centers_sel$V2,"_",k562_centers_sel$V3),dist_k562_to_cage=k562_centers_sel$V10,cage_id=k562_centers_sel$V8)

rm(k562_centers,hek293_centers)

# temp_join_M_distances <- left_join(temp_join_M,hek293_centers_sel_df,by="cage_id") %>%left_join(.,k562_centers_sel_df,by="cage_id")

plot_density_distances <- function(dataframe, var_to_use, title){
  plot_density <- ggplot(temp_join_M_distances,aes(dist_k562_to_cage)) + geom_density() +xlim(c(-200,200))
  density_distr_details_log <- density(temp_join_M_distances$dist_k562_to_cage,na.rm = T)
  ymax <- which.max(density(temp_join_M_distances$dist_k562_to_cage,na.rm=T,n=1000)$y)
  x_to_print <- round(density(temp_join_M_distances$dist_k562_to_cage,na.rm=T,n=1000)$x[ymax],2)
  
  plot_density_with_maxima <- ggplot(temp_join_M_distances,aes(dist_k562_to_cage)) + geom_density() +xlim(c(-200,200)) + geom_vline(xintercept = density(temp_join_M_distances$dist_k562_to_cage,na.rm=T,n=1000)$x[ymax],linetype="dashed",col = 'red') + geom_text(aes(x=x_to_print,y=0,label=x_to_print,col = 'red'), vjust = -0.4,hjust = -0.1)
  
  ggsave(filename="plot_density_with_maxima.pdf",plot=plot_density_with_maxima)
  
  
}
library(ggplot2)

# make a function that produce all plots

#import other data to stratify
ccaatbox <- read.delim("CCAATbox_200_fimo.temp",stringsAsFactors = F)
tmp <- ccaatbox$sequence_name[1:(length(ccaatbox$sequence_name)-3)]
ccaatbox_df <- data.frame(ccaat=rep(1,length(tmp)),cage_id=gsub(":|-","_",gsub("\\(.*\\)","",tmp)))
rm(tmp)

GCbox <- read.delim("GCbox_200_fimo.temp",stringsAsFactors = F)
tmp <- GCbox$sequence_name[1:(length(GCbox$sequence_name)-3)]
GCbox_df <- data.frame(GCbox=rep(1,length(tmp)),cage_id=gsub(":|-","_",gsub("\\(.*\\)","",tmp)))
rm(tmp)

tatabox <- read.delim("TATAbox_50_fimo.temp",stringsAsFactors = F)
tmp <- tatabox$sequence_name[1:(length(tatabox$sequence_name)-3)]
tatabox_df <- data.frame(tatabox=rep(1,length(tmp)),cage_id=gsub(":|-","_",gsub("\\(.*\\)","",tmp)))
rm(tmp)

init <- read.delim("initiator_200_fimo.temp",stringsAsFactors = F)
tmp <- init$sequence_name[1:(length(init$sequence_name)-3)]
init_df <- data.frame(init=rep(1,length(tmp)),cage_id=gsub(":|-","_",gsub("\\(.*\\)","",tmp)))
rm(tmp)

CpG_50 <- read.delim("slop50L.CpG.bed",stringsAsFactors = F,header=F)
CpG_50_df <- data.frame(CpG50=rep(1,length(CpG_50$V1)),cage_id=CpG_50$V5)
CpG_200 <- read.delim("slop200L.CpG.bed",stringsAsFactors = F,header=F)
CpG_200_df <- data.frame(CpG200=rep(1,length(CpG_200$V1)),cage_id=CpG_200$V5)
# slop50L.CpG.bed
# slop200L.CpG.bed

# make a function what estimate density, max and also perform prpoportion test (over representation vs under)
temp_join_M_distances_elements <- left_join(temp_join_M,ccaatbox_df,by='cage_id') %>% left_join(.,GCbox_df,by='cage_id') %>% left_join(.,tatabox_df,by='cage_id') %>% left_join(.,init_df,by='cage_id') %>% left_join(.,CpG_50_df,by='cage_id') %>% left_join(.,CpG_200_df,by='cage_id')

#features
#ccaatbox,GCbox_df,tatabox_df,init_df,CpG_50_df,CpG_200_df

plot_density_distances <- function(dataframe, title_label){

  # dataframe <- temp_to_use
  # title_label <- "k562_repE2"

  plot_density <- ggplot(dataframe,aes(dist)) + geom_density() +xlim(c(-200,200)) + labs(title=title_label)
  ggsave(filename=paste0(title_label,"_plot_density.pdf"))
  density_distr_details_log <- density(dataframe$dist,na.rm = T)
  ymax <- which.max(density(dataframe$dist,na.rm=T,n=1000,bw=0.5)$y)
  x_to_print <- round(density(dataframe$dist,na.rm=T,n=1000,bw=0.5)$x[ymax],2)
  
  plot_density_with_maxima <- ggplot(dataframe,aes(dist)) + geom_histogram(aes(y=..density..), colour="black", fill="grey",bins=20) + geom_density() +xlim(c(-200,200)) + geom_vline(xintercept = density(dataframe$dist,na.rm=T,n=1000)$x[ymax],linetype="dashed",col = 'red') + geom_text(aes(x=x_to_print,y=0,label=x_to_print,col = 'red'), vjust = -0.4,hjust = -0.1)
  
  ggsave(filename=paste0(title_label,"_plot_density_with_maxima.pdf"),plot=plot_density_with_maxima)
  
  #filter dataframe for broad and narrow
  names(dataframe)[grep('width',names(dataframe))] <- "width"
  names(dataframe)[grep('tpm',names(dataframe))] <- "tpm"
  dataframe_updated <- dataframe %>% mutate(type_cage=case_when(width<= 4 ~ "narrow",width> 4 ~ "broad"))
  ref_tpm_quantile <- quantile(dataframe$tpm,c(0.25,0.75))
  dataframe_updated <- dataframe_updated %>% mutate(tpm_class=case_when(tpm<=ref_tpm_quantile[1] ~ "tpm_lower_quart",tpm> ref_tpm_quantile[1] & tpm< ref_tpm_quantile[2] ~"tpm_iqr",tpm > ref_tpm_quantile[2] ~'tpm_upper_quart'))
  dataframe_updated <- dataframe_updated  %>% mutate(dist_to_cage=case_when(dist<= 0 & dist> -200 ~ "upstream200", dist> 0 & dist<200 ~ "downstram200",dist <= -200 | dist > 200 ~ "distal"))
  
  
  
  contingency_table <- dataframe_updated %>% group_by(dist_to_cage,type_cage,ccaat ,GCbox ,tatabox ,init ,CpG200,tpm_class) %>% tally()
  
  write.table(contingency_table, file =paste0(title_label,"_contingency_table.csv"), sep = ",", quote = F, col.names = NA)
  
  dataframe_updated_narrow <- dataframe_updated %>%  filter(type_cage=="narrow")
  ymax_narrow <- which.max(density(dataframe_updated_narrow$dist,na.rm=T,n=1000,bw=0.5)$y)
  x_to_print_narrow <- round(density(dataframe_updated_narrow$dist,na.rm=T,n=1000,bw=0.5)$x[ymax],2)
  ggplot(dataframe_updated_narrow,aes(dist)) + geom_histogram(aes(y=..density..), colour="black", fill="grey",bins=20) + geom_density(adjust = 0.7) +xlim(c(-200,200)) + geom_vline(xintercept = density(dataframe_updated_narrow$dist,na.rm=T,n=1000,bw=0.5)$x[ymax_narrow],linetype="dashed",col = 'red') + geom_text(aes(x=x_to_print_narrow,y=0,label=x_to_print,col = 'red'), vjust = -0.4,hjust = -0.1)
  ggsave(filename=paste0(title_label,"_plot_density_with_maxima_narrow.pdf"))
  
  dataframe_updated_broad <- dataframe_updated %>%  filter(type_cage=="broad")
  ymax_broad <- which.max(density(dataframe_updated_broad$dist,na.rm=T,n=1000,bw=0.5)$y)
  x_to_print_broad <- round(density(dataframe_updated_broad$dist,na.rm=T,n=1000,bw=0.5)$x[ymax],2)
  ggplot(dataframe_updated_broad,aes(dist)) + geom_histogram(aes(y=..density..), colour="black", fill="grey",bins=20) + geom_density(adjust = 0.7) +xlim(c(-200,200)) + geom_vline(xintercept =x_to_print_broad,linetype="dashed",col = 'red') + geom_text(aes(x=x_to_print_broad,y=0,label=x_to_print_broad,col = 'red'), vjust = -0.4,hjust = -0.1)
  ggsave(filename=paste0(title_label,"_plot_density_with_maxima_broad.pdf"))
  
  # include an additional stratification by promoter elements
  
}

# there are 5 additional promoter elements features
temp_df <-temp_join_M_distances_elements %>%  drop_na(starts_with('k562_rep1')) %>% select('cage_id','k562_rep1_tpm','k562_rep1_width',"ccaat","GCbox","tatabox","init","CpG50","CpG200")
head(temp_df)
temp_df$ccaat[is.na(temp_df$ccaat)] <- 0
temp_df$GCbox [is.na(temp_df$GCbox )] <- 0
temp_df$tatabox [is.na(temp_df$tatabox )] <- 0
temp_df$init [is.na(temp_df$init )] <- 0
temp_df$CpG50[is.na(temp_df$CpG50)] <- 0
temp_df$CpG200[is.na(temp_df$CpG200)] <- 0

temp_to_use <- left_join(k562_centers_sel_df,temp_df,by='cage_id') %>% distinct() %>% drop_na()
names(temp_to_use)[grep('dist',names(temp_to_use))] <- paste0("dist")
head(temp_to_use)
plot_density_distances(temp_to_use,"k562_repE1")

temp_df <-temp_join_M_distances_elements %>%  drop_na(starts_with('k562_rep2')) %>% select('cage_id','k562_rep2_tpm','k562_rep2_width',"ccaat","GCbox","tatabox","init","CpG50","CpG200")
head(temp_df)
temp_df$ccaat[is.na(temp_df$ccaat)] <- 0
temp_df$GCbox [is.na(temp_df$GCbox )] <- 0
temp_df$tatabox [is.na(temp_df$tatabox )] <- 0
temp_df$init [is.na(temp_df$init )] <- 0
temp_df$CpG50[is.na(temp_df$CpG50)] <- 0
temp_df$CpG200[is.na(temp_df$CpG200)] <- 0
temp_to_use <- left_join(k562_centers_sel_df,temp_df,by='cage_id') %>% distinct() %>% drop_na()
names(temp_to_use)[grep('dist',names(temp_to_use))] <- paste0("dist")
plot_density_distances(temp_to_use,"k562_repE2")



temp_df <-temp_join_M_distances_elements %>%  drop_na(starts_with('k562_rep3')) %>% select('cage_id','k562_rep3_tpm','k562_rep3_width',"ccaat","GCbox","tatabox","init","CpG50","CpG200")
temp_df$ccaat[is.na(temp_df$ccaat)] <- 0
temp_df$GCbox [is.na(temp_df$GCbox )] <- 0
temp_df$tatabox [is.na(temp_df$tatabox )] <- 0
temp_df$init [is.na(temp_df$init )] <- 0
temp_df$CpG50[is.na(temp_df$CpG50)] <- 0
temp_df$CpG200[is.na(temp_df$CpG200)] <- 0
temp_to_use <- left_join(k562_centers_sel_df,temp_df,by='cage_id') %>% distinct() %>% drop_na()
names(temp_to_use)[grep('dist',names(temp_to_use))] <- paste0("dist")
plot_density_distances(temp_to_use,"k562_repE3")


temp_df <-temp_join_M_distances_elements %>%  drop_na(starts_with('k562_rep0')) %>% select('cage_id','k562_rep0_tpm','k562_rep0_width',"ccaat","GCbox","tatabox","init","CpG50","CpG200")
temp_df$ccaat[is.na(temp_df$ccaat)] <- 0
temp_df$GCbox [is.na(temp_df$GCbox )] <- 0
temp_df$tatabox [is.na(temp_df$tatabox )] <- 0
temp_df$init [is.na(temp_df$init )] <- 0
temp_df$CpG50[is.na(temp_df$CpG50)] <- 0
temp_df$CpG200[is.na(temp_df$CpG200)] <- 0
temp_to_use <- left_join(k562_centers_sel_df,temp_df,by='cage_id') %>% distinct() %>% drop_na()
names(temp_to_use)[grep('dist',names(temp_to_use))] <- paste0("dist")
plot_density_distances(temp_to_use,"k562_rep0")


temp_df <-temp_join_M_distances_elements %>%  drop_na(starts_with('hek293')) %>% select('cage_id','hek293_tpm','hek293_width',"ccaat","GCbox","tatabox","init","CpG50","CpG200")
temp_df$ccaat[is.na(temp_df$ccaat)] <- 0
temp_df$GCbox [is.na(temp_df$GCbox )] <- 0
temp_df$tatabox [is.na(temp_df$tatabox )] <- 0
temp_df$init [is.na(temp_df$init )] <- 0
temp_df$CpG50[is.na(temp_df$CpG50)] <- 0
temp_df$CpG200[is.na(temp_df$CpG200)] <- 0
temp_to_use <- left_join(hek293_centers_sel_df,temp_df,by='cage_id') %>% distinct() %>% drop_na()
names(temp_to_use)[grep('dist',names(temp_to_use))] <- paste0("dist")
plot_density_distances(temp_to_use,"hek293")


