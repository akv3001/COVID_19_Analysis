#!/bin/bash -l
#SBATCH --job-name=RNACov
#SBATCH --output=RNACov_%A_%a.out
#SBATCH --error=RNASeqCov_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=2
#SBATCH --workdir="/athena/elementolab/scratch/akv3001/SLURM_Scripts/SLURM_WorkLog"
#SBATCH --mem=60G

# LOAD SPACK DEPENDENCIES
spack load bzip2@1.0.6+shared%gcc@6.3.0

path=$1 #path to all the samples
gtf_path=$2 #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]
PROJECT_NAME=$(basename "$path")
file=$(ls ${path} | tail -n +${SLURM_ARRAY_TASK_ID}| head -1) # Uses job array for each sample in the folder

PATH_TO_STORE=$(dirname "$path")

cd $TMPDIR  #COPY THE SAMPLES FOLDER



###  Path to Store ###
command_1=${PATH_TO_STORE}/BAMS_${PROJECT_NAME}
command_2=${PATH_TO_STORE}/${PROJECT_NAME}_FPKM
command_3=${PATH_TO_STORE}/${PROJECT_NAME}_Counts

if [[ ! -d "$command_1" ]];then mkdir $command_1; fi
if [[ ! -d "$command_2" ]];then mkdir $command_2; fi
if [[ ! -d "$command_3" ]];then mkdir $command_3 ; fi


echo "$path"

Sample=$(basename "$file")
Sample=${Sample%.*}
echo $Sample
rsync -r -v -a -z --exclude 'Summary' $path/$file/*.gz ./ #copy the samples foldeR
#rsync -r -v -a -z --exclude 'Summary' $path/$file/*.fastq ./
#rsync -r -v -a -z --exclude 'Summary' $path/*.gz ./
#rsync -r -v -a -z --exclude 'Summary' $file ./


echo 'Processing' ${Sample}


#Select gtf based on option
if [ $gtf_path == 1 ]
then
	echo "Aligning and quantifying against Human-Hg19"
	gtf="/athena/elementolab/scratch/akv3001/GenomeReference_hg19/gencode.v19.annotation.gtf"
	STARREF="/athena/elementolab/scratch/akv3001/GenomeReference_hg19/Indexed_genome/hg19_Genomedir"
elif [ $gtf_path == 2 ]
then
	echo "Aligning and quantifying against Mouse-mm10"
	gtf="/athena/elementolab/scratch/akv3001/GenomeReference_mm10/gencode.vM10.annotation.gtf"
	STARREF="/athena/elementolab/scratch/akv3001/GenomeReference_mm10/mm10_Genome_Index"
elif [ $gtf_path == 3 ]
then
        echo "Aligning and quantifying against Human - hg38"
        gtf="/athena/elementolab/scratch/akv3001/GenomeReference_hg38/gencode.v22.annotation.gtf"
        STARREF="/athena/elementolab/scratch/akv3001/GenomeReference_hg38/GRCh38/Sequence/hg38_STAR_index"

else
	gtf=$2
	STARREF=$3
	echo " Aligning to $STARREF "

	if [ $gtf == " "]
	then
		echo "Reference gtf not detected "
		exit 1
	fi

fi



########## Merge #######################3
mkdir $TMPDIR/${Sample}_merged

      echo "Inside $file"

      for R in 1 2; do

              cat $TMPDIR/*_$R*.fastq.gz  > $TMPDIR/${Sample}_merged/"$Sample"_R$R.fastq.gz &

              echo "...$R"

      done

        wait

#########################

# Paired End RNA-Seq
F1=$(ls $TMPDIR/${Sample}_merged/*R1*.gz)
F2=$(ls $TMPDIR/${Sample}_merged/*R2*.gz)



echo "F1 = $F1"
echo "F2 = $F2"

rsync -v $gtf ./


gtf_name=$(basename $gtf)
echo $gtf_name

##-----------------------------------------------------------

	mkdir $TMPDIR/${Sample}_STAR

	cd $TMPDIR/${Sample}_STAR

#------------STAR/tophat Alignment Command--------------------------------#


		echo "[Aligning fastq.gz]"
		STAR --genomeDir ${STARREF} \
		--readFilesIn ${F1}  \
		--outSAMstrandField intronMotif \
		--outFileNamePrefix ${Sample}  \
		--runThreadN 4 --readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate


	echo $(ls -l )

        #--------------Samtools index ------------------------------

        /home/akv3001/Programs/samtools-1.2/samtools index \
	$TMPDIR/${Sample}_STAR/${Sample}*.bam


      java  -Xmx4g -jar /home/akv3001/Programs/picard/dist/picard.jar  AddOrReplaceReadGroups \
      I=$TMPDIR/${Sample}_STAR/*Aligned.sortedByCoord.out.bam \
      O=$TMPDIR/${Sample}_STAR/${Sample}_sorted_RG.bam \
      RGID=1 \
      RGLB=paired \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=${Sample}

      rsync -r -v --exclude="*.sam" $TMPDIR/${Sample}_STAR ${PATH_TO_STORE}/BAMS_${PROJECT_NAME}

 	/home/akv3001/Programs/samtools-1.2/samtools sort \
        -T ${Sample}_name_sorted \
        -o ${Sample}_name_sorted.bam \
        -n $TMPDIR/${Sample}_STAR/${Sample}Aligned.sortedByCoord.out.bam

#fi

echo $(ls -l )

#------------Running HTSeq Count-----------------------------
# a Running htseq-count

samtools index $TMPDIR/${Sample}_STAR/${Sample}_name_sorted.bam
ls -lh

python /home/akv3001/HTSeq-0.6.1/scripts/htseq-count -q -r name \
-s reverse -f bam $TMPDIR/${Sample}_STAR/${Sample}_name_sorted.bam \
$TMPDIR/$gtf_name  > $TMPDIR/${Sample}.bam.count

#-------------------------COPY RESULTS--------------------------------------

rsync -r -v $TMPDIR/${Sample}.bam.count  ${PATH_TO_STORE}/${PROJECT_NAME}_Counts

#-------------------CuffLinks Command---------------------------

cufflinks -q  \
 -p 2 --library-type fr-firststrand -G $TMPDIR/$gtf_name \
 -o $TMPDIR/${Sample}_CuffLinks  $TMPDIR/${Sample}_STAR/${Sample}Aligned.sortedByCoord.out.bam

#------------------ COPY RESULTS ---------------------------------------------
rsync -r -v $TMPDIR/${Sample}_CuffLinks ${PATH_TO_STORE}/${PROJECT_NAME}_FPKM
