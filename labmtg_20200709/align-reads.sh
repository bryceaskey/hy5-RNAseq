#!/bin/bash
#SBATCH --job-name=HISAT2-PRJNA546251       # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=7gb                           # Job memory request
#SBATCH --time=48:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=HISAT2-PRJNA546251_%j.log  # Standard output and error log

pwd; hostname; date

module load adapterremoval/2.2.2 hisat2/2.2.0

echo 'Trimming, filtering, and aligning selected RNAseq data in PRJNA546251 to TAIR10 reference genome'

ref=/ufrc/jkim6/share/braskey/data/TAIR10/
index=/ufrc/jkim6/share/braskey/data/TAIR10/HISAT2-index/
reads=/ufrc/jkim6/share/braskey/data/PRJNA546251/
aln=/ufrc/jkim6/share/braskey/data/PRJNA546251/HISAT2/

# Copy reference fasta into index folder, and generate index
#cp ${ref}TAIR10.fa ${index}
#hisat2-build ${index}TAIR10.fa ${index}TAIR10

for id in SRR9313209 SRR9313210 SRR9313211 SRR9313212 SRR9313213 SRR9313214 SRR9313223 SRR9313224 SRR9313225 SRR9313226 SRR9313227 SRR9313228
do
  # Apply AdapterRemoval to remove adapter sequences and filter low quality reads
  AdapterRemoval --file1 ${reads}${id}.fastq \
    --basename ${reads}${id} --output1 ${reads}${id}_trimmed.fastq \
    --trimns --trimqualities --minlength 80

  # Apply HISAT2 to align trimmed and filtered reads to TAIR10 reference genome
  hisat2 -U ${reads}${id}_trimmed.fastq \
    -x ${index}TAIR10 \
    -S ${aln}${id}.sam
done