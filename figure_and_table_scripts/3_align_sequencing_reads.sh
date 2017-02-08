#!/bin/bash
#PBS -P rnaseq_nod
#PBS -N align-pipe
#PBS -J 1-61
#PBS -j oe
#PBS -q workq
#PBS -o /lustre/scratch/users/michael.schon/log/170202_align/align-rna_^array_index^_perfectmapping.log
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=48gb

# === begin ENVIRONMENT SETUP ===
####set to 0 (false) or 1 (true) to let the repsective code block run
#1. run rsem
run_rsem=1
#2. make plots or not
make_plots=0
#3. delete unecessary files from temp_dir
clean=0
##### specify RSEM parameters
aligner="star"
seq_mode="SE"
file_type="fastq"
threads=8 #set this to the number of available cores
##### specify folders and variables #####
#set script dir
pipe_dir=/lustre/scratch/users/$USER/contamination_pipe
#set ouput base dir
base_dir=/lustre/scratch/users/$USER/rna_seq
#add folder basename as prefix (follows convention from rsem_make_reference)
rsem_ref=$rsem_ref_dir/$(basename $rsem_ref_dir)
#location of the mapping file for the array job
pbs_mapping_file=$pipe_dir/mapping_table.txt
#super folder of the temp dir, script will create subfolders with $sample_name
temp_dir=$base_dir/temp

#####loading of the required modules #####
module load Python/2.7.12-foss-2016b
module load RSEM/1.2.31-foss-2016a
module load BEDTools/v2.17.0-goolf-1.4.10
module load SAMtools/1.3-foss-2015b
# conditional loading of modules based on aligner to be used by RSEM
if [ $aligner == "bowtie" ]; then
  module load Bowtie/1.1.2-foss-2015b
fi
if [ $aligner == "bowtie2" ]; then
  module load Bowtie2/2.2.7-foss-2015b
fi
if [ $aligner == "star" ]; then
  module load rna-star/2.5.2a-foss-2016a
fi
if [ $make_plots -eq 1 ]; then
  module load R/3.2.3-foss-2016a
fi
##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as line number
input_mapper=`sed -n "${PBS_ARRAY_INDEX} p" $pbs_mapping_file` #read mapping file
names_mapped=($input_mapper)
sample_dir=${names_mapped[1]} # get the sample dir
sample_name=`basename $sample_dir` #get the base name of the dir as sample name
#folder for rsem reference
rsem_ref_dir=/lustre/scratch/users/$USER/indices/rsem/$aligner/${names_mapped[2]}
star_folder_name=star_u_cross

#print some output for logging
echo '#########################################################################'
echo 'Starting RSEM RNA-seq pipeline for: '$sample_name
echo 'Rsem reference: ' $rsem_ref
echo 'Aligner to be used: ' $aligner
echo 'Mapping file: ' $pbs_mapping_file
echo 'Selected file type: ' $file_type
echo 'Selected sequencing mode: ' $seq_mode
echo '#########################################################################'

#some paramter checking
if [ $seq_mode != "PE" ] && [ $seq_mode != "SE" ]; then
  echo "Wrong parameters selected for seq_mode! Aborting." 1>&2
  exit 1
fi
if [ $file_type != "bam" ] && [ $file_type != "fastq" ]; then
  echo "Wrong parameters selected for file_type! Aborting." 1>&2
  exit 1
fi

#make output folder
mkdir -p $sample_dir/$star_folder_name/
cd $sample_dir

#folders for temp files
temp_dir_s=$temp_dir/$sample_name
#rm -rf $temp_dir_s

if [ $run_rsem -eq 1 ]; then
  #TODO: think about how to replace the ugly ifs with a case switch
  #initalize variable

  # TEMPORARY FIX: add --no-qualities to rsem_opts. REENABLE AFTER RUNNING NODINE 2012
#######
  rsem_opts=""
#######
if [ $seq_mode = "PE" ]; then #add paired-end flag if data is PE
    rsem_opts=$rsem_opts"--paired-end "
  fi
  if [ $file_type = "bam" ]; then #convert to fastq if input is bam
    f=($(ls  $sample_dir | grep -e ".bam")) # get all bam files in folder
    file_number=${#f[@]} # get length of the array
    if [ "$file_number" = "1" ]; then
      #sort bam file
      samtools sort -n -m 4G -@ $threads -o $sample_dir/${f%.*}.sorted.bam \
      $sample_dir/$f
      if [ $seq_mode = "PE" ]; then
        #convert bam to fastq then add to rsem_opts string
        bedtools bamtofastq -i $sample_dir/${f%.*}.sorted.bam \
          -fq $sample_dir/${f%.*}.1.fq \
          -fq2 $sample_dir/${f%.*}.2.fq
        rsem_opts=$rsem_opts"$sample_dir/${f%.*}.1.fq $sample_dir/${f%.*}.2.fq"
      fi
      if [ $seq_mode = "SE" ]; then
        #convert bam to fastq then add to rsem_opts string
        bedtools bamtofastq -i $sample_dir/${f%.*}.sorted.bam \
          -fq $sample_dir/${f%.*}.fq
        rsem_opts=$rsem_opts"$sample_dir/${f%.*}.fq "
      fi
    else
      echo "Only one bam file per sample folder allowed! Aborting."\
           "Files present: $file_number" 1>&2
      exit 1
    fi
  elif [ $file_type = "fastq" ]; then
    rsem_opts=$rsem_opts
    #check if fastq files are zipped and unzip them if needed
    f=($(ls  $sample_dir | grep -e ".fq.gz\|.fastq.gz"))
    file_number=${#f[@]}
    if [ $file_number -eq 1 ] || [ $file_number -eq 2 ]; then
      cd $sample_dir
      gunzip ${f[@]}
    fi
    #get files with .fq or .fastq extention
    f=($(ls  $sample_dir| grep -e ".fq\|.fastq"))
    file_number=${#f[@]}
    #some error handling. Check if only the expected number of fq files is there
    if [ $file_number -eq 1 ]  && [ "$seq_mode" = "SE" ]; then
      rsem_opts=$rsem_opts"$sample_dir/$f"
    elif [ $file_number -eq 2 ]  && [ "$seq_mode" = "PE" ]; then
      rsem_opts=$rsem_opts"$sample_dir/${f[0]} $sample_dir/${f[1]}"
    else
      echo "Wrong number of fastq files in sample folder! Aborting."\
           "Files present: $file_number" 1>&2
      exit 1
    fi
  else
    echo "Unsupported file type selected! Aborting." 1>&2
    exit 1
  fi

# run rsem to calculate the expression levels
# --estimate-rspd: estimate read start position to check if the data has bias
# --output-genome-bam: output bam file as genomic, not transcript coordinates
# --seed 12345 set seed for reproducibility of rng
# --calc-ci calcutates 95% confidence interval of the expression values
# --ci-memory 30000 set memory
# --solexa-quals only for the 8cell pilot studies
rsem_params="--$aligner \
--num-threads $threads \
--temporary-folder $temp_dir_s \
--append-names \
--estimate-rspd \
--output-genome-bam \
--sort-bam-by-coordinate \
--seed 12345 \
--calc-ci \
--ci-memory 40000 \
--solexa-quals \
$rsem_opts \
$rsem_ref \
$sample_name"

star_params="--runMode alignReads \
--runThreadN $threads \
--runRNGseed 12345 \
--genomeDir $rsem_ref_dir \
--readFilesIn $sample_dir/u_"$sample_name".fastq \
--clip3pNbases 2 \
--clip5pNbases 2 \
--outFileNamePrefix $sample_dir/$star_folder_name/ \
--outTmpDir $temp_dir_s \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMmultNmax 1 \
--outFilterMultimapNmax 1 \
--outFilterMismatchNmax 0 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--quantMode GeneCounts"
#rsem command that should be run

### FASTQ UNIQUE
echo "python $pipe_dir/fastq_unique.py $sample_dir/$(basename $(find $sample_dir/$sample_name*)) $sample_dir/u_"$sample_name".fastq"
eval "python $pipe_dir/fastq_unique.py $sample_dir/$(basename $(find $sample_dir/$sample_name*)) $sample_dir/u_"$sample_name".fastq"

cd $sample_dir/$star_folder_name/
echo "STAR $star_params >& $sample_name.log"
eval "STAR $star_params >& $sample_name.log"
fi

#run the rsem plot function
if [ $make_plots -eq 1 ]; then
  rsem-plot-model $sample_dir/rsem/$sample_name $sample_dir/rsem/$sample_name.pdf
fi

#delete the temp files
if [ $clean -eq 1 ]; then
  gzip $sample_dir/*.fq $sample_dir/*.fastq
  rm $sample_dir/${f%.*}.sorted.bam
  rm $sample_dir/rsem/*.transcript.bam
  rm -rf $temp_dir_s
fi

rm $sample_dir/u_"$sample_name".fastq
echo 'Finished RSEM RNA-seq pipeline for: '$sample_name
