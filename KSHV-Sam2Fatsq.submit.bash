#!/bin/sh
#SBATCH --job-name=KSHV
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=20:00:00
#SBATCH --mem=20gb
#SBATCH --output=KSHV.%J.out
#SBATCH --error=KSHV.%J.err


mail -s "KSHV2 ..." amirsalar.mansouri@gmail.com <<< "Starting ... "

Type1=hg19
GTFhg19=GTFs/gencode.v17.annotation.gtf

Type2=KSHV
GFFKSHV=GFFs/NC_009333.gff

for i in ../Samples/Old/*; do
	if [ -d $i ]; then
    # if the folder is there
		Samplename=${i#../Samples/Old/}
		echo Sample: $Samplename
		mkdir -p Results/$Samplename/SamToFastq
		mkdir -p Results/$Samplename/SamToFastq/TopHat
		
		samtools collate -o Results/$Samplename/TopHat/Sortedunmapped.bam Results/$Samplename/TopHat/unmapped.bam
		samtools fastq -0 Results/$Samplename/SamToFastq/SE.fastq -n Results/$Samplename/TopHat/unmapped.bam
		
		echo TopHat $Samplename $Type2 ...
		tophat2 -o Results/$Samplename/SamToFastq/TopHat -p $SLURM_NTASKS_PER_NODE ../Genomes/$Type2/Index$Type2/$Type2 \
		Results/$Samplename/SamToFastq/SE.fastq				
		
		samtools view -h -o Results/$Samplename/SamToFastq/TopHat/accepted_hits.sam Results/$Samplename/SamToFastq/TopHat/accepted_hits.bam
		samtools index -b Results/$Samplename/SamToFastq/TopHat/accepted_hits.bam Results/$Samplename/SamToFastq/TopHat/accepted_hits.bai
		
		echo Cufflinks $Type2 - $Samplename ...
		cufflinks --GTF ../Genomes/$Type2/$GFFKSHV -o Results/$Samplename/Cufflinks_output/$Type2 -p $SLURM_NTASKS_PER_NODE Results/$Samplename/SamToFastq/TopHat/accepted_hits.sam
	fi
done

echo Cuffdiff $Type1  ...

cuffdiff -L NC,GFP,RGpos -p $SLURM_NTASKS_PER_NODE -o Results/Cuffdiff/$Type1/ ../Genomes/$Type1/$GTFhg19 \
            Results/NC/TopHat/accepted_hits.sam Results/G/TopHat/accepted_hits.sam Results/RG2/TopHat/accepted_hits.sam

echo Cuffdiff $Type2  ...

cuffdiff -L NC,GFP,RGpos -p $SLURM_NTASKS_PER_NODE -o Results/Cuffdiff/$Type2/ ../Genomes/$Type2/$GFFKSHV \
            Results/NC/SamToFastq/TopHat/accepted_hits.sam \
			Results/G/SamToFastq/TopHat/accepted_hits.sam Results/RG2/SamToFastq/TopHat/accepted_hits.sam

echo Cuffdiff $Type2  ...

cuffdiff -L GFP,RGpos -p $SLURM_NTASKS_PER_NODE -o Results/Cuffdiff/$Type2/ ../Genomes/$Type2/$GFFKSHV1 \
            Results/G/SamToFastq/TopHat/accepted_hits.sam Results/RG2/SamToFastq/TopHat/accepted_hits.sam

mail -s "KSHV ..." amirsalar.mansouri@gmail.com <<< "Done..."
