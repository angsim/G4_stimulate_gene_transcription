Yuqi_feb2022
================
Angela Simeone
16/02/2022

## Download data

``` bash
mkdir /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip

java -jar /Users/simeon01/bin/clarity-tools.jar -l SLX-21315
```

# Processing

## rename files

``` bash
## ==  rename files
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355
cut -d ',' -f 1,2,4 SLX-18355.HVN7MBGXB.s_1.contents.csv | sed 's/"//g' | awk  -F',' '{print $1"."$2".HVN7MBGXB.s_1.r_1\t" $3"_"$1"."$2}' | sed 's/-i/_i/g' >  SLX-18355.HVN7MBGXB.s_1.contents_for_renaming.csv
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355
while read -r old_name new_name; do
  echo $old_name
  ls -lth ${old_name}.fq.gz
  echo $new_name
  #ls -lth ${old_name}.bam
  #mv ${old_name}.flagstat ${new_name}.flagstat
  mv ${old_name}.fq.gz ${new_name}.fq.gz
done < /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/SLX-18355.HVN7MBGXB.s_1.contents_for_renaming.csv

## === cut adapters ===

#cut illumina adapters
out_dir=/scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355
cd $out_dir
for f in *.fq.gz 
do
  sbatch --mem 8G --wrap "cutadapt -m 10 -q 20 -O 3 -a CTGTCTCTTATACACATCT -o ${f%%.fastq.gz}.trimmed.fastq.gz $f"
done


## === alignment ===
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355
ls
path=/scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355
g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'
mkdir aligned
for f in *trimmed.fastq.gz
do
    sbatch --mem 16000 --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > $path/aligned/${f%%.trimmed.fastq.gz}.hg38.sort.bam"
done

## === remove duplicates ===
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned
for file in *.hg38.sort.bam
do
  #/Users/marsic01/applications/picard-2.8.2.jar MarkDuplicates ==> this was the one used before
  # now picard in home
  #picard MarkDuplicates --version
  #2.18.12-SNAPSHOT
  sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$file OUTPUT=${file%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${file%%.bam}.markduplicates_metrics.txt"
  
  echo "==========="
  echo "          "
done

## === generate stats ===

for file in *.hg38.sort.bam
  do
  
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam_hg38
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

# === generate tracks === 

cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned
for bam in *.markduplicates.bam
    do
  sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

for bam in *.markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        echo "=============="
done

# === call peaks === 
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned
out_dir_narrow='macs2_individual_rep'
out_dir_broad='macs2_broad_individual_rep'
mkdir $out_dir_narrow
mkdir $out_dir_broad

rep=(1 3 4)
cases=(CNCC_Rep2_3xBG4)
#StemCell3_NT_BG4 # this one has rep2-3-4 and not 1-2-3
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.hg38.sort.markduplicates.bam`
    c=`ls ${j}_input*.hg38.sort.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done

rep=(2 3 4)
cases=(ESC_Rep2_3xBG4)
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.hg38.sort.markduplicates.bam`
    c=`ls ${j}_input*.hg38.sort.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done

rep=(1 2 4)
cases=(NSC_Rep2_3xBG4)
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.hg38.sort.markduplicates.bam`
    c=`ls ${j}_input*.hg38.sort.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done
```

## Processing results of peak calling - intersecting with OQs map and using

On local machine

``` r
compare_tech_rep_overlap_and_OQs <- function(conditions,list_FOI,OQs_map,search_thresholds) {
  # list_FOI is the list of files of interest
  summary_OQs <- c()
  counter <- 1
  temp_label <- c()
  for (j in 1:length(search_thresholds)) { 
    
    for (i in 1:length(conditions)){
      
      curr_condition <- conditions[i]
      
      if (search_thresholds[j] != 'default') {
        
        matches_temp1 <- unique(grep(curr_condition,list_FOI, value=TRUE))
        matches_temp2 <- unique(grep(search_thresholds[j],matches_temp1, value = TRUE))
      
        }  else   {
        
        matches_temp1 <- unique(grep(curr_condition,list_FOI, value=TRUE))
        matches_temp2 <- unique(matches_temp1[!grepl(paste(search_thresholds[(1+1):length(search_thresholds)], collapse ="|"),matches_temp1)])
        }
      
      out_file_name_multi2 <- paste0(curr_condition,'_',search_thresholds[j],'_multi2.bed')
      out_file_name_multi3 <- paste0(curr_condition,'_',search_thresholds[j],'_multi3.bed')
        
      system(paste("multiIntersectBed -i", paste(matches_temp2, collapse=' ') ,"| awk '{if($4>=2) print $0}'| sortBed -i - | mergeBed -i - >",out_file_name_multi2))
      system(paste("multiIntersectBed -i", paste(matches_temp2, collapse=' ') ,"| awk '{if($4>=3) print $0}'| sortBed -i - | mergeBed -i - >",out_file_name_multi3))
      multi2_overlap_OQs <- as.numeric(system(paste("intersectBed -a ",out_file_name_multi2, "-b ",OQs_map," -wa | sortBed -i - | mergeBed -i - | wc -l"),intern = TRUE))
      multi3_overlap_OQs <- as.numeric(system(paste("intersectBed -a ",out_file_name_multi3, "-b ",OQs_map," -wa | sortBed -i - | mergeBed -i - | wc -l"),intern = TRUE)) 
      len_multi2 <- as.numeric(system(paste("wc -l",out_file_name_multi2,"| awk '{print $1}'"),intern = TRUE))
      len_multi3 <- as.numeric(system(paste("wc -l",out_file_name_multi3,"| awk '{print $1}'"),intern = TRUE))
      summ <- rbind(len_multi2,len_multi3,multi2_overlap_OQs,multi3_overlap_OQs,multi2_overlap_OQs/len_multi2,multi3_overlap_OQs/len_multi3)
      print(j)
      print(i)
      summary_OQs <- cbind(summary_OQs,summ)
      temp_label <- c(temp_label,paste0(curr_condition,"_",search_thresholds[j]))
      counter <- counter+1
      
      }
  }
  rownames(summary_OQs) <- c('len_multi2','len_multi3','multi2_overlap_OQs','multi3_overlap_OQs','perc_multi2_overlap_OQs','perc_multi3_overlap_OQs')
  colnames(summary_OQs) <- temp_label
  return(summary_OQs)
} 

## ************************************************
# narrowPeaks
path_peaks <- '/Users/simeon01/Documents/Katie/Sequencing_Sept19/macs2_individual_rep'
setwd(path_peaks)

list_files_narrow_peaks <- list.files(path = path_peaks, pattern= '_peaks.narrowPeak') 

# for human print peak number and save it into a structure that we can plot/export
N_peaks_hg38 <- vector(mode="numeric", length=length(list_files_narrow_peaks))

for (i in 1:length(list_files_narrow_peaks)){
  N_peaks_hg38[i] <- as.numeric(system(paste("wc -l",list_files_narrow_peaks[i],"| awk '{print $1}' "),intern = TRUE))
}
names(N_peaks_hg38) <- list_files_narrow_peaks

# write output to file
write.table(N_peaks_hg38,'SLX-18496_narrowPeak_summary_number_peaks.csv',col.names=NA,row.names = T,quote = F, sep = ",")

# for each condition check multiintersect (2ou3 and 3out3) and then generate dthe intersection with oqs map
conditions_hg38 <- c("ESC_Rep1_3xBG4_rep", "NSC_Rep1_3xBG4_rep","CNCC_Rep1_3xBG4_rep")
pval_hg38 <- c('default')
OQs_hg38 <- '/Users/simeon01/Documents/Katie/OQ_hits.lifted_hg19_to_hg38.bed'
sumary_OQs_hg38 <- compare_tech_rep_overlap_and_OQs(conditions_hg38,list_files_narrow_peaks,OQs_hg38,pval_hg38)
write.table(t(sumary_OQs_hg38),'SLX-18496_summary_confirmed_peaks_overlap_OQS_hg38.csv',col.names=NA,row.names = T,quote = F, sep = ",")


## ************************************************
# broad peaks
path_peaks <- '/Users/simeon01/Documents/Katie/Sequencing_Sept19/macs2_broad_individual_rep'
setwd(path_peaks)
# get filenames of narrowPeaks - all -  even various conditions tested.
list_files_broad_peaks <- list.files(path = path_peaks, pattern= '_peaks.broadPeak') 
# for human print peak number and save it into a structure that we can plot/export
N_peaks_hg38 <- vector(mode="numeric", length=length(list_files_broad_peaks))
for (i in 1:length(list_files_broad_peaks)){
  N_peaks_hg38[i] <- as.numeric(system(paste("wc -l",list_files_broad_peaks[i],"|  awk '{print $1}' "),intern = TRUE))
}
names(N_peaks_hg38) <- list_files_broad_peaks
# write output to file
write.table(N_peaks_hg38,'SLX-18496_broadPeak_summary_number_peaks.csv',col.names=NA,row.names = T,quote = F, sep = ",")

# for each condition check multiintersect (2ou3 and 3out3) and then generate dthe intersection with oqs map
conditions_hg38 <- c("ESC_Rep1_3xBG4_rep", "NSC_Rep1_3xBG4_rep","CNCC_Rep1_3xBG4_rep")
pval_hg38 <- c('default')
OQs_hg38 <- '/Users/simeon01/Documents/Katie/OQ_hits.lifted_hg19_to_hg38.bed'
sumary_OQs_hg38 <- compare_tech_rep_overlap_and_OQs(conditions_hg38,list_files_broad_peaks,OQs_hg38,pval_hg38)
write.table(t(sumary_OQs_hg38),'SLX-18496_summary_confirmed_broadpeaks_overlap_OQS_hg38.csv',col.names=NA,row.names = T,quote = F, sep = ",")
```

# Second Biological Replicate

## Dowload fastq from basespace and start processing second biological replicate

``` bash
bs list runs
#+--------------------------------+-----------+----------------------------------------------------------------+----------+
#|              Name              |    Id     |                         ExperimentName                         |  Status  |
#+--------------------------------+-----------+----------------------------------------------------------------+----------+
#| 190913_NS500222_0532_HC2M5BGXC | 190873736 | 20190913_SLX18495_StemCell_Rep3_3xBG4                          | Complete |
#| 190903_NS500222_0531_HVN3KBGXB | 190706553 | 20190903_SLX-18351_PDTX                                        | Complete |


bs list projects
#+-------------------------------------------------+-----------+--------------+
#|                      Name                       |    Id     |  TotalSize   |
#+-------------------------------------------------+-----------+--------------+
#| 20190913_SLX18495_StemCell_Rep3                 | 141038899 | 16717067092  |
#| 20190903_SLX-18351_PDTX_NEWchromatin_3xupscale  | 140676543 | 22975479257  |

katie_path=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4
mkdir $katie_path
cd $katie_path

sbatch --mem 18G --wrap "bs download project -i 141038899 -o $katie_path --extension=fastq.gz"

# replace names 

for f in `find . -name '*.gz'`; do f2=${f##./}; f3=${f2/\//_}; f4=${f3/ds.*SLX/SLX}; mv $f $f4; done

# remove all the files that are not fastq and start processing
```

## Processing

``` bash
## === cut adapters ===

#cut illumina adapters
out_dir=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4
cd $out_dir
for f in *.fastq.gz 
do
  sbatch --mem 8G --wrap "cutadapt -m 10 -q 20 -O 3 -a CTGTCTCTTATACACATCT -o ${f%%.fastq.gz}.trimmed.fastq.gz $f"
done

## === alignment ===
cd /scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4
ls
path=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4
g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'
mkdir aligned
for f in *trimmed.fastq.gz
do
    sbatch --mem 16000 --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > $path/aligned/${f%%.trimmed.fastq.gz}.hg38.sort.bam"
done


## === merge lanes of the same library ==
cd /scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned
out_folder=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged
mkdir $out_folder

for f in *_L001*.hg38.sort.bam
do
  base_name=${f%%_SLX*}
  echo $base_name
  
  hg38_f=`ls ${base_name}*hg38.sort.bam`
  
  echo $hg38_f
  ls ${hg38_f%%_L001*}*hg38.sort.bam
  echo ${hg38_f%%_L001*}.all.hg38.merged.bam
  
  #echo sbatch --mem 8G --wrap "samtools merge -@8 $out_folder/${hg38_f%%_L001*}.all.hg38.merged.bam ${hg38_f%%_L001*}*hg38.sort.bam"
  sbatch --mem 8G --wrap "samtools merge -@8 $out_folder/${hg38_f%%_L001*}.all.hg38.merged.bam ${hg38_f%%_L001*}*hg38.sort.bam"
  
  #echo "
  echo " ======================================= "
  echo " "
done

# chip1 from ESC escaped (because wrongly renamed) 
# manual merging
sbatch --mem 8G --wrap "samtools merge -@8 $out_folder/ESC_Rep3_3xBG4_ChIP1.all.hg38.merged.bam ESC_Rep3_3xBG4_ChIP1_L00*"


## === remove duplicates ===
cd /scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged
for file in *.merged.bam
do
  #/Users/marsic01/applications/picard-2.8.2.jar MarkDuplicates ==> this was the one used before
  # now picard in home
  #picard MarkDuplicates --version
  #2.18.12-SNAPSHOT
  sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$file OUTPUT=${file%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${file%%.bam}.markduplicates_metrics.txt"
  
  echo "==========="
  echo "          "
done

## === generate stats ===

for file in *.hg38.merged.bam
  do
  
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam_hg38
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done


touch SLX-18355.out.stat
for bam in *hg38.sort.bam
  do
  base_name=${bam/.bam/}
  files_stat=`ls *.stat* | grep $base_name`
  echo $bam > temp
  cat $files_stat >> temp
  paste -d '\t' SLX-18355.out.stat temp >> temp2
  mv temp2 SLX-18355.out.stat
done

# === generate tracks === 

cd /scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged
for bam in *.markduplicates.bam
    do
  sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

for bam in *.markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        echo "=============="
done



# === call peaks === 
cd /scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged
out_dir_narrow='macs2_individual_rep'
out_dir_broad='macs2_broad_individual_rep'
mkdir $out_dir_narrow
mkdir $out_dir_broad

rep=(2 3 4)
cases=(CNCC_Rep3_3xBG4 NSC_Rep3_3xBG4)

for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.all.hg38.merged.markduplicates.bam`
    c=`ls ${j}_input.all.hg38.merged.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done

rep=(1 2 4)
cases=(ESC_Rep3_3xBG4)
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.all.hg38.merged.markduplicates.bam`
    c=`ls ${j}_input.all.hg38.merged.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done
```

## Processing results of peak calling - intersecting with OQs map and using

On local machine

``` r
compare_tech_rep_overlap_and_OQs <- function(conditions,list_FOI,OQs_map,search_thresholds) {
  # list_FOI is the list of files of interest
  summary_OQs <- c()
  counter <- 1
  temp_label <- c()
  for (j in 1:length(search_thresholds)) { 
    
    for (i in 1:length(conditions)){
      
      curr_condition <- conditions[i]
      
      if (search_thresholds[j] != 'default') {
        
        matches_temp1 <- unique(grep(curr_condition,list_FOI, value=TRUE))
        matches_temp2 <- unique(grep(search_thresholds[j],matches_temp1, value = TRUE))
      
        }  else   {
        
        matches_temp1 <- unique(grep(curr_condition,list_FOI, value=TRUE))
        matches_temp2 <- unique(matches_temp1[!grepl(paste(search_thresholds[(1+1):length(search_thresholds)], collapse ="|"),matches_temp1)])
        }
      
      out_file_name_multi2 <- paste0(curr_condition,'_',search_thresholds[j],'_multi2.bed')
      out_file_name_multi3 <- paste0(curr_condition,'_',search_thresholds[j],'_multi3.bed')
        
      system(paste("multiIntersectBed -i", paste(matches_temp2, collapse=' ') ,"| awk '{if($4>=2) print $0}'| sortBed -i - | mergeBed -i - >",out_file_name_multi2))
      system(paste("multiIntersectBed -i", paste(matches_temp2, collapse=' ') ,"| awk '{if($4>=3) print $0}'| sortBed -i - | mergeBed -i - >",out_file_name_multi3))
      multi2_overlap_OQs <- as.numeric(system(paste("intersectBed -a ",out_file_name_multi2, "-b ",OQs_map," -wa | sortBed -i - | mergeBed -i - | wc -l"),intern = TRUE))
      multi3_overlap_OQs <- as.numeric(system(paste("intersectBed -a ",out_file_name_multi3, "-b ",OQs_map," -wa | sortBed -i - | mergeBed -i - | wc -l"),intern = TRUE)) 
      len_multi2 <- as.numeric(system(paste("wc -l",out_file_name_multi2,"| awk '{print $1}'"),intern = TRUE))
      len_multi3 <- as.numeric(system(paste("wc -l",out_file_name_multi3,"| awk '{print $1}'"),intern = TRUE))
      summ <- rbind(len_multi2,len_multi3,multi2_overlap_OQs,multi3_overlap_OQs,multi2_overlap_OQs/len_multi2,multi3_overlap_OQs/len_multi3)
      print(j)
      print(i)
      summary_OQs <- cbind(summary_OQs,summ)
      temp_label <- c(temp_label,paste0(curr_condition,"_",search_thresholds[j]))
      counter <- counter+1
      
      }
  }
  rownames(summary_OQs) <- c('len_multi2','len_multi3','multi2_overlap_OQs','multi3_overlap_OQs','perc_multi2_overlap_OQs','perc_multi3_overlap_OQs')
  colnames(summary_OQs) <- temp_label
  return(summary_OQs)
} 

## ************************************************
# narrowPeaks
path_peaks <- '/Users/simeon01/Documents/Katie/Sequencing_Sept19_thirdBioRep/macs2_individual_rep'
setwd(path_peaks)

list_files_narrow_peaks <- list.files(path = path_peaks, pattern= '_peaks.narrowPeak') 

# for human print peak number and save it into a structure that we can plot/export
N_peaks_hg38 <- vector(mode="numeric", length=length(list_files_narrow_peaks))

for (i in 1:length(list_files_narrow_peaks)){
  N_peaks_hg38[i] <- as.numeric(system(paste("wc -l",list_files_narrow_peaks[i],"| awk '{print $1}' "),intern = TRUE))
}
names(N_peaks_hg38) <- list_files_narrow_peaks

# write output to file
write.table(N_peaks_hg38,'SLX-18355_narrowPeak_summary_number_peaks.csv',col.names=NA,row.names = T,quote = F, sep = ",")

# for each condition check multiintersect (2ou3 and 3out3) and then generate dthe intersection with oqs map
conditions_hg38 <- c("ESC_Rep2_3xBG4_rep", "NSC_Rep2_3xBG4_rep","CNCC_Rep2_3xBG4_rep")
pval_hg38 <- c('default')
OQs_hg38 <- '/Users/simeon01/Documents/Katie/OQ_hits.lifted_hg19_to_hg38.bed'
sumary_OQs_hg38 <- compare_tech_rep_overlap_and_OQs(conditions_hg38,list_files_narrow_peaks,OQs_hg38,pval_hg38)
write.table(t(sumary_OQs_hg38),'SLX-18355_summary_confirmed_peaks_overlap_OQS_hg38.csv',col.names=NA,row.names = T,quote = F, sep = ",")


## ************************************************
# broad peaks
path_peaks <- '/Users/simeon01/Documents/Katie/Sequencing_Sept19_thirdBioRep/macs2_broad_individual_rep'
setwd(path_peaks)
# get filenames of narrowPeaks - all -  even various conditions tested.
list_files_broad_peaks <- list.files(path = path_peaks, pattern= '_peaks.broadPeak') 
# for human print peak number and save it into a structure that we can plot/export
N_peaks_hg38 <- vector(mode="numeric", length=length(list_files_broad_peaks))
for (i in 1:length(list_files_broad_peaks)){
  N_peaks_hg38[i] <- as.numeric(system(paste("wc -l",list_files_broad_peaks[i],"|  awk '{print $1}' "),intern = TRUE))
}
names(N_peaks_hg38) <- list_files_broad_peaks
# write output to file
write.table(N_peaks_hg38,'SLX-18355_broadPeak_summary_number_peaks.csv',col.names=NA,row.names = T,quote = F, sep = ",")

# for each condition check multiintersect (2ou3 and 3out3) and then generate dthe intersection with oqs map
conditions_hg38 <- c("ESC_Rep2_3xBG4_rep", "NSC_Rep2_3xBG4_rep","CNCC_Rep2_3xBG4_rep")
pval_hg38 <- c('default')
OQs_hg38 <- '/Users/simeon01/Documents/Katie/OQ_hits.lifted_hg19_to_hg38.bed'
sumary_OQs_hg38 <- compare_tech_rep_overlap_and_OQs(conditions_hg38,list_files_broad_peaks,OQs_hg38,pval_hg38)
write.table(t(sumary_OQs_hg38),'SLX-18355_summary_confirmed_broadpeaks_overlap_OQS_hg38.csv',col.names=NA,row.names = T,quote = F, sep = ",")
```

``` bash
oqs=/Users/simeon01/Documents/references/OQs/OQ_hits.lifted_hg19_to_hg38.bed
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19/macs2_individual_rep

for file in *narrowPeak
do
intersectBed -a $file -b $oqs -wa | sort | uniq | wc -l
done
```

## Merge Replicate that has been resequemced

Katie has resequenced a library prepared before that shower low depth.

I 1) merge the aligned bam (prior dup-removal) 2) remove duplicates 3)
call narrow (3A) and broad(3B) peaks. Use them for NSC.

``` bash
NSC_bam_seq1=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged/NSC_Rep3_3xBG4_ChIP2.all.hg38.merged.bam 
NSC_bam_seq2=/scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned/NSC_Rep3b_3xBG4_ChIP2_SLX-18355.i708_i505.hg38.sort.bam
NSC_bam_merged=/scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned/NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam
NSC_bam_merged_nodup=/scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned/NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.markduplicates.bam

## == mergeBam == 
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned
sbatch --mem 16G --wrap="samtools merge NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam $NSC_bam_seq1 $NSC_bam_seq2 && samtools index NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam"

## == remove dup == 
sbatch --mem 16G  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam OUTPUT=NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.markduplicates_metrics.txt"

## == generate stats ==
sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $NSC_bam_merged >  ${NSC_bam_merged%%.bam}.stat2"
sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $NSC_bam_merged_nodup >  ${NSC_bam_merged_nodup%%.bam}.stat5"

## == call peaks only for this rep ==
# note: input has been sequenced before! /scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged/NSC_Rep3_3xBG4_input.all.hg38.merged.markduplicates.bam
cd /scratcha/sblab/simeon01/Data/20220216_SLX-21315_yuqi_BG4Chip/SLX-18355/aligned
out_dir_narrow='macs2_individual_rep'
out_dir_broad='macs2_broad_individual_rep'
NSC_rep3_input=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged/NSC_Rep3_3xBG4_input.all.hg38.merged.markduplicates.bam
sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $NSC_bam_merged_nodup -c $NSC_rep3_input -n $out_dir_narrow/NSC_Rep3b_and_3_3xBG4_rep3" # human is default
sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $NSC_bam_merged_nodup -c $NSC_rep3_input -n $out_dir_broad/NSC_Rep3b_and_3_3xBG4_rep3_hg38_broad" # human is default
```
