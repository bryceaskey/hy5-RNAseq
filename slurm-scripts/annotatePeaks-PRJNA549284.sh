#!/bin/bash
#SBATCH --job-name=annotatePeaks-PRJNA549284        # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                 # Where to send mail	
#SBATCH --account=jkim6                             # Group providing CPU and memory resources
#SBATCH --qos=jkim6                                 # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                                  # Number of CPU cores to use
#SBATCH --mem=4gb                                   # Job memory request
#SBATCH --time=24:00:00                             # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=annotatePeaks-PRJNA549284_%j.log   # Standard output and error log

pwd; hostname; date

echo "Finding genes associated with ChIPseq peaks"

module load homer/4.10 bedtools/2.29.2

peaks=/blue/jkim6/share/braskey/data/PRJNA549284/peaks/
bed=/blue/jkim6/share/braskey/data/PRJNA549284/bed/
bedAnnotated=/blue/jkim6/share/braskey/data/PRJNA549284/bedAnnotated/
ref=/blue/jkim6/share/braskey/data/TAIR10/

for id in hy5_1 hy5_2 hy5_3 hy5-VP16_1 hy5-VP16_2 hy5-VP16_3 hy5-SRDX_1 hy5-SRDX_2
do
  pos2bed.pl ${peaks}${id}.txt > ${bed}${id}.bed
  bedtools window -a ${bed}${id}.bed -b ${ref}TAIR10_genes-only.gff \
    -l 0 -r 2000 -header > ${bedAnnotated}${id}.txt
done