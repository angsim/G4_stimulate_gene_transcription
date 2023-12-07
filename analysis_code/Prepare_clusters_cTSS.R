library(CAGEr)
data(FANTOM5humanSamples)
setwd('/Users/simeon01/Documents/Yuqi/Analasis_distance_Cage_G4')

#look for K562
grep('K562',FANTOM5humanSamples$data_url, value = T)

#[1] "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/chronic%2520myelogenous%2520leukemia%2520cell%2520line%253aK562.CNhs11250.10454-106G4.hg19.ctss.bed.gz"                              
#[2] "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/chronic%2520myelogenous%2520leukemia%2520cell%2520line%253aK562%2520ENCODE%252c%2520biol_rep1.CNhs12334.10824-111C5.hg19.ctss.bed.gz"
#[3] "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/chronic%2520myelogenous%2520leukemia%2520cell%2520line%253aK562%2520ENCODE%252c%2520biol_rep2.CNhs12335.10825-111C6.hg19.ctss.bed.gz"
#[4] "http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.hCAGE/chronic%2520myelogenous%2520leukemia%2520cell%2520line%253aK562%2520ENCODE%252c%2520biol_rep3.CNhs12336.10826-111C7.hg19.ctss.bed.gz"


#select K562
K562Samples <- FANTOM5humanSamples[grep("K562", 
                                        FANTOM5humanSamples[,"description"]),]

#import TSS data using importPublicData() function and set the argument source = "FANTOM5"
# The dataset argument can be set to either "human" or "mouse", and the sample argument is provided by a vector of sample lables/names. 

K562Samples[,"sample"]
# [1] "chronic_myelogenous_leukemia_cell_line_K562"                  
# [2] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep1"
# [3] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep2"
# [4] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep3"

K562CAGEset <- importPublicData(source = "FANTOM5", dataset = "human", 
                                sample = K562Samples[1:4,"sample"])

# K562CAGEset contains the output of fetching the samples from fantom5
corr.m <- plotCorrelation2( K562CAGEset, samples = "all"
                            , tagCountThreshold = 1, applyThresholdBoth = FALSE
                            , method = "pearson")


cum_dis <- cumulativeCTSSdistribution(K562CAGEset, clusters = "tagClusters", useMulticore = T)

quantilePositions(K562CAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
test_aggregate <- aggregateTagClusters(K562CAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)


normalizeTagCount(K562CAGEset, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)
# one clustering for each experiment
K562CAGEset2 <- K562CAGEset
clusterCTSS( object = K562CAGEset
             , threshold = 1
             , thresholdIsTpm = TRUE
             , nrPassThreshold = 1
             , method = "distclu"
             , maxDist = 20
             , removeSingletons = TRUE
             , keepSingletonsAbove = 5)

cumulativeCTSSdistribution(K562CAGEset, clusters = "tagClusters", useMulticore = T)

quantilePositions(K562CAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

#ranges of clusters for each case
# [1] "chronic_myelogenous_leukemia_cell_line_K562"                  
# [2] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep1"
# [3] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep2"
# [4] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep3"
k562_individ <- tagClustersGR( K562CAGEset, "chronic_myelogenous_leukemia_cell_line_K562"
               , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
k562_individE1 <- tagClustersGR( K562CAGEset, "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep1"
                               , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
k562_individE2 <- tagClustersGR( K562CAGEset, "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep2"
                               , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
k562_individE3 <- tagClustersGR( K562CAGEset, "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep3"
                               , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

write.table(k562_individ,'k562_individ.csv',sep = ",",quote = F, col.names = NA)
write.table(k562_individE1,'k562_individE1.csv',sep = ",",quote = F, col.names = NA)
write.table(k562_individE2,'k562_individE2.csv',sep = ",",quote = F, col.names = NA)
write.table(k562_individE3,'k562_individE3.csv',sep = ",",quote = F, col.names = NA)


# # this generates the ensembl consensus
aggregateTagClusters(K562CAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
k562_aggreg <- consensusClustersGR(K562CAGEset)
write.table(k562_aggreg,'k562_aggreg_consensus.csv',sep = ",",quote = F, col.names = NA)

# this generates the consensus
cumulativeCTSSdistribution(K562CAGEset, clusters = "consensusClusters", useMulticore = TRUE)
quantilePositions(K562CAGEset, clusters = "consensusClusters", qLow = 0.1, qUp = 0.9, useMulticore = TRUE)
k562_aggreg_all  <- consensusClustersGR( K562CAGEset, 
                     returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)




# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
HEKSamples <- FANTOM5humanSamples[grep("HEK293%252fSLAM%2520untreated", 
                                        FANTOM5humanSamples[,"data_url"]),]
HEKCAGEset <- importPublicData(source = "FANTOM5", dataset = "human", 
                                sample = HEKSamples[1,"sample"])
normalizeTagCount(HEKCAGEset, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)
clusterCTSS( object = HEKCAGEset
             , threshold = 1
             , thresholdIsTpm = TRUE
             , nrPassThreshold = 1
             , method = "distclu"
             , maxDist = 20
             , removeSingletons = TRUE
             , keepSingletonsAbove = 5)

cumulativeCTSSdistribution(HEKCAGEset, clusters = "tagClusters", useMulticore = T)

quantilePositions(HEKCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
HEK293_individ <- tagClustersGR( HEKCAGEset, "embryonic_kidney_cell_line__HEK293_SLAM_untreated"
                              , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)

write.table(HEK293_individ,'hek293_cage_clusters.csv',sep = ",",quote = F, col.names = NA)

par(mfrow=c(2,1))
hist(k562_individE3$interquantile_width, 100,freq=F, ylim=c(0,0.2),xlim=c(0,50), density = F)
hist(HEK293_individ$interquantile_width, 100,freq=F, ylim=c(0,0.2),xlim=c(0,50))

d_k562_individE3 <- density(k562_individE3$interquantile_width)
d_HEK293_individ <- density(HEK293_individ$interquantile_width)

plot(d_k562_individE3, y=c(0,0.10),xlim=c(0,50))
plot(d_HEK293_individ,y=c(0,0.10),xlim=c(0,50))

library('fitdistrplus')
plotdist(k562_individE3$interquantile_width,histo = TRUE, demp = TRUE, ylim=c(0,0.10),)
plotdist(HEK293_individ$interquantile_width,histo = TRUE, demp = TRUE)

par(mfrow=c(1,1))
descdist(k562_individE3$interquantile_width, boot = 1000)
fw <- fitdist(k562_individE3$interquantile_width, "weibull")
fg <- fitdist(k562_individE3$interquantile_width, "gamma")
fln <- fitdist(k562_individE3$interquantile_width, "lnorm")

summary(fw)

par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)


gofstat(list(fw, fg, fln), fitnames = c("weibull", "gamma", "lnorm"))


par(mfrow=c(1,1))
descdist(HEK293_individ$interquantile_width, boot = 1000)
fw <- fitdist(HEK293_individ$interquantile_width, "weibull")
fg <- fitdist(HEK293_individ$interquantile_width, "gamma")
fln <- fitdist(HEK293_individ$interquantile_width, "lnorm")

summary(fw)

par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)


gofstat(list(fw, fg, fln), fitnames = c("weibull", "gamma", "lnorm"))

library('mixtools')
# https://hal.archives-ouvertes.fr/hal-00717545/file/RR_Chauveau.pdf
# https://tinyheero.github.io/2015/10/13/mixture-model.html

set.seed(12)
a_K562 <- npEM(k562_individE1$interquantile_width, mu0=2, bw=1.5,maxiter=150)

par(mfrow = c(1,1))
plot(a_K562)

set.seed(12)
a_hek293 <- npEM(HEK293_individ$interquantile_width, mu0=2, bw=1.5,maxiter=150)

par(mfrow = c(1,1))
plot(a_hek293)
          
#periodicitiy
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4207034/


## try to read hg38
clusterCTSS( object = K562CAGEset2
             , threshold = 1
             , thresholdIsTpm = TRUE
             , nrPassThreshold = 1
             , method = "distclu"
             , maxDist = 10
             , removeSingletons = TRUE
             , keepSingletonsAbove = 2)

cumulativeCTSSdistribution(K562CAGEset2, clusters = "tagClusters", useMulticore = T)

quantilePositions(K562CAGEset2, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

#ranges of clusters for each case
# [1] "chronic_myelogenous_leukemia_cell_line_K562"                  
# [2] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep1"
# [3] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep2"
# [4] "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep3"
k562_individ_sec <- tagClustersGR( K562CAGEset2, "chronic_myelogenous_leukemia_cell_line_K562"
                               , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
k562_individE1_sec <- tagClustersGR( K562CAGEset2, "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep1"
                                 , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
k562_individE2_sec <- tagClustersGR( K562CAGEset2, "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep2"
                                 , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
k562_individE3_sec <- tagClustersGR( K562CAGEset2, "chronic_myelogenous_leukemia_cell_line_K562_ENCODE__biol_rep3"
                                 , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
HEKCAGEset2 <- HEKCAGEset

clusterCTSS( object = HEKCAGEset2
             , threshold = 1
             , thresholdIsTpm = TRUE
             , nrPassThreshold = 1
             , method = "distclu"
             , maxDist = 10
             , removeSingletons = TRUE
             , keepSingletonsAbove = 1)
cumulativeCTSSdistribution(HEKCAGEset2, clusters = "tagClusters", useMulticore = T)

quantilePositions(HEKCAGEset2, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
HEK293_individ_sec <- tagClustersGR( HEKCAGEset2, "embryonic_kidney_cell_line__HEK293_SLAM_untreated"
                                 , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)


HEKCAGEset3 <- HEKCAGEset

clusterCTSS( object = HEKCAGEset3
             , threshold = 1
             , thresholdIsTpm = FALSE
             , nrPassThreshold = 1
             , method = "distclu"
             , maxDist = 10
             , removeSingletons = TRUE
             , keepSingletonsAbove = 1)

cumulativeCTSSdistribution(HEKCAGEset3, clusters = "tagClusters", useMulticore = T)

quantilePositions(HEKCAGEset3, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
HEK293_individ_third <- tagClustersGR( HEKCAGEset3, "embryonic_kidney_cell_line__HEK293_SLAM_untreated"
                                     , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)


HEKCAGEset4 <- HEKCAGEset

clusterCTSS( object = HEKCAGEset4
             , threshold = 1
             , thresholdIsTpm = FALSE
             , nrPassThreshold = 1
             , method = "distclu"
             , maxDist = 10
             , removeSingletons = TRUE
             , keepSingletonsAbove = 3)

cumulativeCTSSdistribution(HEKCAGEset4, clusters = "tagClusters", useMulticore = T)

quantilePositions(HEKCAGEset4, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
HEK293_individ_forth <- tagClustersGR( HEKCAGEset4, "embryonic_kidney_cell_line__HEK293_SLAM_untreated"
                                       , returnInterquantileWidth = TRUE,  qLow = 0.1, qUp = 0.9)
