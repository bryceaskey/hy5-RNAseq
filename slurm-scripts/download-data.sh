#!/bin/bash
#SBATCH --job-name=download-data        # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu     # Where to send mail	
#SBATCH --account=jkim6                 # Group providing CPU and memory resources
#SBATCH --qos=jkim6                     # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                      # Number of CPU cores to use
#SBATCH --mem=1gb                       # Job memory request
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=download-data_%j.log   # Standard output and error log

pwd; hostname; date

module load sra/2.10.4

echo "Downloading arabidopsis RNA seq data"

dest=/ufrc/jkim6/share/braskey/data/PRJNA549285/

# PRJNA549285
fastq-dump SRR9313209 -O ${dest} # Col-0, dark 3 days
fastq-dump SRR9313210 -O ${dest} # Col-0, dark 3 days
fastq-dump SRR9313211 -O ${dest} # Col-0, dark 3 days

fastq-dump SRR9313212 -O ${dest} # hy5, dark 3 days
fastq-dump SRR9313213 -O ${dest} # hy5, dark 3 days
fastq-dump SRR9313214 -O ${dest} # hy5, dark 3 days

fastq-dump SRR9313223 -O ${dest} # Col-0, dark 3 days + light 1.5 hrs
fastq-dump SRR9313224 -O ${dest} # Col-0, dark 3 days + light 1.5 hrs
fastq-dump SRR9313225 -O ${dest} # Col-0, dark 3 days + light 1.5 hrs

fastq-dump SRR9313226 -O ${dest} # hy5, dark 3 days + light 1.5 hrs
fastq-dump SRR9313227 -O ${dest} # hy5, dark 3 days + light 1.5 hrs
fastq-dump SRR9313228 -O ${dest} # hy5, dark 3 days + light 1.5 hrs