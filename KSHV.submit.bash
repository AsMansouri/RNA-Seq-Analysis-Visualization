#!/bin/sh
#SBATCH --job-name=KSHV
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=40:00:00
#SBATCH --mem=20gb
#SBATCH --output=KSHV.%J.out
#SBATCH --error=KSHV.%J.err

module load samtools/1.3 bowtie/2.2 tophat/2.1 cufflinks/2.2 trimmomatic/0.33 fastqc/0.10

mail -s "KSHV ..." amirsalar.mansouri@gmail.com <<< "Starting ... "

Type1=hg19
GTFhg19=GTFs/gencode.v17.annotation.gtf

Type2=KSHV1
GFFKSHV1=GFFs/NC_009333.gff

bowtie2-build -f ../Genomes/$Type2/NC_009333.fna ../Genomes/$Type2/Index$Type2/$Type2

Type3=KSHV2
GFFKSHV2=GFFs/GQ994935.gtf

set a="First"
set lanes=0
for i in ../Samples/Old/*; do
	if [ -d $i ]; then
    # if the folder is there
		Samplename=${i#../Samples/Old/}
		echo Sample: $Samplename
		mkdir -p Results/$Samplename
		mkdir -p Results/$Samplename/FastQC
		mkdir -p Results/$Samplename/Trimmomatic
		mkdir -p Results/$Samplename/Trimmomatic/FastQC
		mkdir -p Results/$Samplename/Bowtie
		mkdir -p Results/$Samplename/TopHat
		mkdir -p Results/$Samplename/SamToFastq
		mkdir -p Results/$Samplename/SamToFastq/TopHat
		mkdir -p Results/$Samplename/Cufflinks_output/$Type1
		mkdir -p Results/$Samplename/Cufflinks_output/$Type2
		for tmpr in ../Samples/Old/$Samplename/*; do
			if [ -f $tmpr ]; then
				tmprr=../Samples/Old/$Samplename/
				r=${tmpr#$tmprr}

				let "lanes += 1"
				fastqc -o Results/$Samplename/FastQC/ --quiet -t $SLURM_NTASKS_PER_NODE -f fastq ../Samples/Old/$Samplename/$r
				read=${r/.fastq.gz/}
				
				echo FastQC $Samplename - $read
				fastqcFol=$read
				fastqcFol+="_fastqc"
				echo $fastqcFol
				python Trims.py -i Results/$Samplename/FastQC/$fastqcFol/fastqc_data.txt -o Results/$Samplename/FastQC/Adapters
				echo Extreacting the Over representated reads of $Samplename - $read 


				eval allreads[$lanes]=$read
					
				a="First"
				read1=$read
				echo read1 : $read1
				readU1=$read1
				readU1+="_U"
				
				echo Trimmimg the reads...
				java -jar $TM_HOME/trimmomatic.jar SE -threads $SLURM_NTASKS_PER_NODE -phred33 -trimlog Results/$Samplename/Trimmomatic/Trimmomatictrimlog \
				../Samples/Old/$Samplename/$read1.fastq.gz Results/$Samplename/Trimmomatic/$read1.fastq \
				ILLUMINACLIP:Results/$Samplename/FastQC/Adapters.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
				adapt=Results/$Samplename/FastQC/Adapters.fa
				rm $adapt
				
				fastqc -o Results/$Samplename/Trimmomatic/FastQC/ --quiet -t $SLURM_NTASKS_PER_NODE -f fastq Results/$Samplename/Trimmomatic/$read1.fastq
				# set r!lanes!=$read

				if [ $((lanes%4)) -eq 0 ]; then
					lanes=0
						
					echo TopHat $Samplename ...
					tophat2 -o Results/$Samplename/TopHat -p $SLURM_NTASKS_PER_NODE ../Genomes/$Type1/Index$Type1/$Type1 \
					Results/$Samplename/Trimmomatic/${allreads[1]}.fastq,Results/$Samplename/Trimmomatic/${allreads[2]}.fastq,Results/$Samplename/Trimmomatic/${allreads[3]}.fastq,Results/$Samplename/Trimmomatic/${allreads[4]}.fastq
					
					samtools view -h -o Results/$Samplename/TopHat/accepted_hits.sam Results/$Samplename/TopHat/accepted_hits.bam
					samtools index -b Results/$Samplename/TopHat/accepted_hits.bam Results/$Samplename/TopHat/accepted_hits.bai
					samtools view -h -o Results/$Samplename/TopHat/unmapped.sam Results/$Samplename/TopHat/unmapped.bam
						
					echo Cufflinks $Type1 - $Samplename ...
					cufflinks --GTF ../Genomes/$Type1/$GTFhg19 -o Results/$Samplename/Cufflinks_output/$Type1 -p $SLURM_NTASKS_PER_NODE Results/$Samplename/TopHat/accepted_hits.sam
	
								
					
					samtools fastq -1 Results/$Samplename/SamToFastq/paired1.fastq Results/$Samplename/TopHat/unmapped.bam
		
					echo TopHat $Samplename $Type2 ...
					tophat2 -o Results/$Samplename/SamToFastq/TopHat -p $SLURM_NTASKS_PER_NODE ../Genomes/$Type2/Index$Type2/$Type2 \
					Results/$Samplename/SamToFastq/paired1.fastq			
		
					samtools view -h -o Results/$Samplename/SamToFastq/TopHat/accepted_hits.sam Results/$Samplename/SamToFastq/TopHat/accepted_hits.bam
					samtools index -b Results/$Samplename/SamToFastq/TopHat/accepted_hits.bam Results/$Samplename/SamToFastq/TopHat/accepted_hits.bai
		
					echo Cufflinks $Type2 - $Samplename ...
					cufflinks --GTF ../Genomes/$Type2/$GFFKSHV1 -o Results/$Samplename/Cufflinks_output/$Type2 -p $SLURM_NTASKS_PER_NODE Results/$Samplename/SamToFastq/TopHat/accepted_hits.sam

				fi
			
			fi
		echo allreads : ${allreads[*]}
		done
	fi
done

mkdir -p Results/Cuffdiff
mkdir -p Results/Cuffdiff/$Type1
mkdir -p Results/Cuffdiff/$Type2

echo Cuffdiff $Type1  ...

cuffdiff -L NC,GFP,RGpos -p $SLURM_NTASKS_PER_NODE -o Results/Cuffdiff/$Type1/ ../Genomes/$Type1/$GTFhg19 \
            Results/NC/TopHat/accepted_hits.sam Results/GFP/TopHat/accepted_hits.sam Results/RGpos/TopHat/accepted_hits.sam

echo Cuffdiff $Type2  ...

cuffdiff -L NEG,GFP,RGpos -p $SLURM_NTASKS_PER_NODE -o Results/Cuffdiff/$Type2/ ../Genomes/$Type2/$GFFKSHV1 \
            Results/NEG/SamToFastq/TopHat/accepted_hits.sam Results/GFP/SamToFastq/TopHat/accepted_hits.sam Results/RGpos/SamToFastq/TopHat/accepted_hits.sam


mail -s "KSHV ..." amirsalar.mansouri@gmail.com <<< "Done..."
