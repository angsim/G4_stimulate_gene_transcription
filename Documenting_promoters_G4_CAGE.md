Distance CAGE G4
================

## Resources and papers

<https://www.sciencedirect.com/science/article/pii/S0022283619302530>

refTSS: A Reference Data Set for Human and Mouse Transcription Start
Sites

<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4581360/> ElemeNT: a
computational tool for detecting core promoter elements

<https://www.sciencedirect.com/science/article/abs/pii/S0378111906006238?via%3Dihub>
Prevalence of the initiator over the TATA box in human and yeast genes
and identification of DNA motifs enriched in human TATA-less core
promoters – pioneering work about TATA box in human and yeast

<https://www.frontiersin.org/articles/10.3389/fgene.2019.00286/full#B46>
DeePromoter: Robust Promoter Predictor Using Deep Learning

–\> eukariotic promoter database
<https://epd.epfl.ch/human/human_database.php?db=human>

–FindM - to find motifs

Promoter selection tool <https://epd.epfl.ch/get_promoters.php>

After identifying regions (cluster) of interest extract only 50 and
200bp upstreat start site (stranded) and get fasta

``` bash
fasta=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa
genome=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome

cd /Users/simeon01/Documents/Yuqi/Analasis_distance_Cage_G4

awk -F "," '{print $2"\t"$3"\t"$4"\t.\t.\t"$6}' k562_individE3.csv  | tail -n +2 | sortBed -i - | slopBed -i - -l 50 -r 0 -g $genome -s | bedtools getfasta -fi $fasta -bed - -s >test.fa


 for file in *.csv ; do awk -F "," '{print $2"\t"$3"\t"$4"\t.\t.\t"$6}' $file  | tail -n +2 | sortBed -i - | uniq| sort > ${file%%.csv}.bed; done
 cat hek293_cage_clusters.bed k562_aggreg_consensus.bed k562_individ.bed k562_individE1.bed k562_individE2.bed  k562_individE3.bed |  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$1"_"$2"_"$3"\t"$6}' > all_cage_prom.bed
 
 

 
fasta=/scratcha/sblab/simeon01/reference_genomes/hg19_gen/hsa.hg19.fa
genome=/scratcha/sblab/simeon01/reference_genomes/hg19_gen/hsa.hg19.genome
cd /scratcha/sblab/simeon01/Data/20201126_yuqi_promoters/CAGE_clusters_cTSS
 # extract both TSS-50 and TSS-200
sortBed -i all_cage_prom.bed | slopBed -i - -l 50 -r 0 -g $genome -s |sort | uniq | bedtools getfasta -fi $fasta -bed - -s > all_cage_prom.slop50L.fa
sortBed -i all_cage_prom.bed | slopBed -i - -l 50 -r 0 -g $genome -s |sort | uniq > all_cage_prom.slop50L.bed
sortBed -i all_cage_prom.bed | slopBed -i - -l 200 -r 0 -g $genome -s |sort | uniq | bedtools getfasta -fi $fasta -bed - -s > all_cage_prom.slop200L.fa
sortBed -i all_cage_prom.bed | slopBed -i - -l 200 -r 0 -g $genome -s |sort | uniq >  all_cage_prom.slop200L.bed


subtractBed -b all_cage_prom.bed -a all_cage_prom.slop50L.bed > slop50L.bed
subtractBed -b all_cage_prom.bed -a all_cage_prom.slop200L.bed > slop200L.bed

cat slop50L.bed | sort | uniq |bedtools getfasta -fi $fasta -bed - -s > slop50L.fa
cat slop200L.bed | sort | uniq |bedtools getfasta -fi $fasta -bed - -s > slop200L.fa

# overlap with promoters and CpG

cut -f 2,3,4 cpgIslandExt.hg19.txt | sortBed -i - | intersectBed -a slop50L.bed -b - -wa > slop50L.CpG.bed
cut -f 2,3,4 cpgIslandExt.hg19.txt | sortBed -i - | intersectBed -a slop200L.bed -b - -wa > slop200L.CpG.bed


# check relative distance of BG4 to promoters
k562_G4=/scratcha/sblab/simeon01/Data/20201126_yuqi_promoters/g4_k562/trimmed/aligned/GSE107690_K562_High_confidence_peaks.bed
hek293_hg19=/scratcha/sblab/simeon01/Data/20201126_yuqi_promoters/KASseq_data/Hek.bio2.sorted.hg19.bed
cd /scratcha/sblab/simeon01/Data/20201126_yuqi_promoters/CAGE_clusters_cTSS
awk '{print $1"\t"$2+int(($3-$2)/2+0.49)-1"\t"$2+int(($3-$2)/2+0.49)}' $k562_G4 | sortBed -i - > K562.centers.hg19.bed
awk '{print $1"\t"$2+int(($3-$2)/2+0.49)-1"\t"$2+int(($3-$2)/2+0.49)}' $hek293_hg19 | sortBed -i - > hek293.centers.hg19.bed

# to find relative distance metric - specific metric
bedtools reldist -a K562.centers.hg19.bed -b all_cage_prom.bed > relative_dist.k562_G4.all_cage_prom.bed
bedtools reldist -a hek293.centers.hg19.bed -b all_cage_prom.bed > relative_dist.hek293_G4.all_cage_prom.bed


#ctctf as a ref
bedtools reldist -a ENCFF736NYC_ctcf_hg19.sorted.bed -b all_cage_prom.bed > relative_dist.ENCFF736NYC_ctcf.all_cage_prom.bed



# to find relative distance metrix between centers of regions and beginning of prom
awk '{if ($6=="+") {print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6} else{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6} }' all_cage_prom.bed |sortBed -i - | bedtools closest -b - -a K562.centers.hg19.bed -D "b" -t "first"| sort | uniq > all_cage_prom.K562.centers.hg19.bed

awk '{if ($6=="+") {print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6} else{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6} }' all_cage_prom.bed |sortBed -i - |  bedtools closest -b - -a hek293.centers.hg19.bed -D "b" -t "first" | sort | uniq > all_cage_prom.hek293.centers.hg19.bed


# process FIMO output
for file in *temp
do
cat $file | sed 's/:/\t/g' | sed 's/(-)/\tminus/g'| sed 's/(+)/\tplus/g' | sed 's/-/\t/g' | sed 's/minus/-/g'| sed 's/plus/+/g'| tail -n +2 | head -n -4  | sortBed -i -   > ${file%%.temp}.fin.bed
done

# annotate promoters cage



1. left join all_cage_prom and all individual cells (merging by cluster id where cluster id=chr_start_end)
2. left join with presence of not of any of the studied motifs (the tmp file contains the coordinate of the promoter, it needs to be parsed 
3. left join the output of closest (*centers.hg19.bed)

1b . independently, plot deistibution of relative distances. 

---> plots all possible cases (stratifying, shyny?)
```

``` r
import:

all_cage_prom
k562_aggreg
hek293_cage_clusters.csv
k562_individ.csv
k562_individE1.csv
k562_individE2.csv
k562_individE3.csv

all_cage_prom.hek293.centers.hg19.bed
all_cage_prom.K562.centers.hg19.bed

relative_dist.ENCFF736NYC_ctcf.all_cage_prom.bed
relative_dist.k562_G4.all_cage_prom.bed
relative_dist.hek293_G4.all_cage_prom.bed

CCAATbox_200_fimo.temp
GCbox_200_fimo.temp
TATAbox_50_fimo.temp
initiator_200_fimo.temp

slop50L.CpG.bed
slop200L.CpG.bed
```
