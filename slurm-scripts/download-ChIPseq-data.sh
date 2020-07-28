#!/bin/bash
#SBATCH --job-name=download-ChIPseq-data        # Job name
#SBATCH --mail-type=END,FAIL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu             # Where to send mail	
#SBATCH --account=jkim6                         # Group providing CPU and memory resources
#SBATCH --qos=jkim6                             # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                              # Number of CPU cores to use
#SBATCH --mem=1gb                               # Job memory request
#SBATCH --time=24:00:00                         # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=download-ChIPseq-data_%j.log   # Standard output and error log

pwd; hostname; date

module load sra/2.10.4

echo 'Downloading arabidopsis ChIP-seq data'

dest=/blue/jkim6/share/braskey/data/PRJNA549284/

# hy5
fastq-dump SRR9312923 -O ${dest}
fastq-dump SRR9312924 -O ${dest}
fastq-dump SRR10811458 -O ${dest}

# hy5-VP16
fastq-dump SRR9312925 -O ${dest}
fastq-dump SRR9312926 -O ${dest}
fastq-dump SRR10811459 -O ${dest}

# hy5-SRDX
fastq-dump SRR9312927 -O ${dest}
fastq-dump SRR9312928 -O ${dest}

# controls
fastq-dump SRR9312922 -o ${dest}
fastq-dump SRR10811456 -o ${dest}
fastq-dump SRR10811457 -o ${dest}