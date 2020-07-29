#!/bin/bash
#SBATCH --job-name=findPeaks-PRJNA549284        # Job name
#SBATCH --mail-type=END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu             # Where to send mail	
#SBATCH --account=jkim6                         # Group providing CPU and memory resources
#SBATCH --qos=jkim6                             # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                              # Number of CPU cores to use
#SBATCH --mem=4gb                               # Job memory request
#SBATCH --time=24:00:00                         # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=findPeaks-PRJNA549284_%j.log   # Standard output and error log

pwd; hostname; date

module load homer/4.10

echo "Finding peaks from ChIPseq data for PRJNA549284"

aln=/blue/jkim6/share/braskey/data/PRJNA549284/HISAT2/
tags=/blue/jkim6/share/braskey/data/PRJNA549284/tags/
peaks=/blue/jkim6/share/braskey/data/PRJNA549284/peaks/

# Make tag directories for all samples
for id in SRR9312923 SRR9312924 SRR10811458 SRR9312925 SRR9312926 SRR10811459 SRR9312927 SRR9312928 SRR9312922 SRR10811456 SRR10811457
do
  makeTagDirectory ${tags}${id}/ ${aln}${id}.sam
done

# hy5 peak identification
findPeaks ${tags}SRR9312923/ -i ${tags}SRR9312922/ -o ${peaks}hy5_1.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${tags}SRR9312924/ -i ${tags}SRR10811456/ -o ${peaks}hy5_2.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${tags}SRR10811458/ -i ${tags}SRR10811457/ -o ${peaks}hy5_3.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

# hy5-VP16 peak identification
findPeaks ${tags}SRR9312925/ -i ${tags}SRR9312922/ -o ${peaks}hy5-VP16_1.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${tags}SRR9312926/ -i ${tags}SRR9312922/ -o ${peaks}hy5-VP16_2.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${tags}SRR10811459/ -i ${tags}SRR10811457/ -o ${peaks}hy5-VP16_3.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

# hy5-SRDX peak identification
findPeaks ${tags}SRR9312927/ -i ${tags}SRR9312922/ -o ${peaks}hy5-SRDX_1.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${tags}SRR9312928/ -i ${tags}SRR9312922/ -o ${peaks}hy5-SRDX_2.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001