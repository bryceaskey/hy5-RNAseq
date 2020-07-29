#!/bin/bash
#SBATCH --job-name=convert-to-BAM-PRJNA549284        # Job name
#SBATCH --mail-type=END,FAIL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                  # Where to send mail	
#SBATCH --account=jkim6                              # Group providing CPU and memory resources
#SBATCH --qos=jkim6                                  # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                                   # Number of CPU cores to use
#SBATCH --mem=4gb                                    # Job memory request
#SBATCH --time=24:00:00                              # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=convert-to-BAM-PRJNA549284_%j.log   # Standard output and error log

pwd; hostname; date

echo "Finding genes associated with ChIPseq peaks"

module load samtools/1.10

sam="/blue/jkim6/share/braskey/data/PRJNA549284/HISAT2/sam/"
bam="/blue/jkim6/share/braskey/data/PRJNA549284/HISAT2/bam/"

samtools view -b ${sam}SRR9312923.sam > ${bam}hy5_1.bam
samtools view -b ${sam}SRR9312924.sam > ${bam}hy5_2.bam
samtools view -b ${sam}SRR10811458.sam > ${bam}hy5_3.bam
samtools view -b ${sam}SRR9312925.sam > ${bam}hy5-VP16_1.bam
samtools view -b ${sam}SRR9312926.sam > ${bam}hy5-VP16_2.bam
samtools view -b ${sam}SRR10811459.sam > ${bam}hy5-VP16_3.bam
samtools view -b ${sam}SRR9312927.sam > ${bam}hy5-SRDX_1.bam
samtools view -b ${sam}SRR9312928.sam > ${bam}hy5-SRDX_2.bam

for id in hy5_1 hy5_2 hy5_3 hy5-VP16_1 hy5-VP16_2 hy5-VP16_3 hy5-SRDX_1 hy5-SRDX_2
do
  samtools sort ${bam}${id}.bam -o ${bam}${id}.sorted.bam
  samtools index ${bam}${id}.sorted.bam
done