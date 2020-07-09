#!/bin/bash
#SBATCH --job-name=align-reads              # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=...@ufl.edu             # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=1gb                           # Job memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=align-reads_%j.log         # Standard output and error log

pwd; hostname; date

module load hisat2/2.2.0

echo 'Aligning reads to TAIR10 reference genome'

index=/ufrc/jkim6/...
reads=/ufrc/jkim6/...
aln=/ufrc/jkim6/...
mkdir -p ${aln}

# Generate index for HISAT2
hisat2-build ${index}TAIR10.fa ${index}TAIR10

for id in ...
do
  # Apply HISAT2 to align trimmed and filtered reads to TAIR10 reference genome
  hisat2 -U ${reads}${id}_trimmed.fastq \
    -x ${index}TAIR10 \
    -S ${aln}${id}.sam
done