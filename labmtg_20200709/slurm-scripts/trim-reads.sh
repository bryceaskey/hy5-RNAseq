#!/bin/bash
#SBATCH --job-name=trim-reads               # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=...@ufl.edu             # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=1gb                           # Job memory request
#SBATCH --time=24:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=trim-reads_%j.log          # Standard output and error log

pwd; hostname; date

module load adapterremoval/2.2.2

echo 'Trimming and filtering reads'

reads=/ufrc/jkim6/...

for id in ...
do
  # Apply AdapterRemoval to remove adapter sequences and filter low quality reads
  AdapterRemoval --file1 ${reads}${id}.fastq \
    --basename ${reads}${id} --output1 ${reads}${id}_trimmed.fastq \
    --trimns --trimqualities --minlength 36
done