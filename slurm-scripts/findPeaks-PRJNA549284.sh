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
peaks=/blue/jkim6/share/braskey/data/PRJNA549284/peaks/

# hy5
findPeaks ${aln}SRR9312923.sam -i ${aln}SRR9312922.sam -o ${peaks}hy5_1.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${aln}SRR9312924.sam -i ${aln}SRR10811456.sam -o ${peaks}hy5_2.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${aln}SRR10811458.sam -i ${aln}SRR10811457.sam -o ${peaks}hy5_3.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

# hy5-VP16
findPeaks ${aln}SRR9312925.sam -i ${aln}SRR9312922.sam -o ${peaks}hy5-VP16_1.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${aln}SRR9312926.sam -i ${aln}SRR9312922.sam -o ${peaks}hy5-VP16_2.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${aln}SRR10811459.sam -i ${aln}SRR10811457.sam -o ${peaks}hy5-VP16_3.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

# hy5-SRDX
findPeaks ${aln}SRR9312927.sam -i ${aln}SRR9312922.sam -o ${peaks}hy5-SRDX_1.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001

findPeaks ${aln}SRR9312928.sam -i ${aln}SRR9312922.sam -o ${peaks}hy5-SRDX_2.txt \
  -style factor -F 3 -P 0.0001 -L 2 -LP 0.001