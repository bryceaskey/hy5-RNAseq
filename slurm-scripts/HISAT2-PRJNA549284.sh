#!/bin/bash
#SBATCH --job-name=HISAT2-PRJNA549284       # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=48:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=HISAT2-PRJNA549284_%j.log  # Standard output and error log

pwd; hostname; date

module load trimmomatic/0.39 hisat2/2.2.0

echo 'Aligning selected ChipSeq data from PRJNA549284 to TAIR10 reference genome'

ref=/blue/jkim6/share/braskey/data/TAIR10/
index=/blue/jkim6/share/braskey/data/TAIR10/HISAT2-index/
reads=/blue/jkim6/share/braskey/data/PRJNA549284/
aln=/blue/jkim6/share/braskey/data/PRJNA549284/HISAT2/

# Copy reference fasta into index folder, and generate index
#cp ${ref}TAIR10.fa ${index}
#hisat2-build ${index}TAIR10.fa ${index}TAIR10

for id in SRR9312923 SRR9312924 SRR10811458 SRR9312925 SRR9312926 SRR10811459 SRR9312927 SRR9312928 SRR9312922 SRR10811456 SRR10811457
do
  if [ "${id}" == "SRR9312924" ] || [ "${id}" == "SRR10811456" ]
  then
    minLen=80
  else
    minLen=36
  fi

  # Apply trimmomatic to clip adapter sequences and filter low quality reads
  trimmomatic SE -threads 1 \
    ${reads}${id}.fastq ${reads}${id}_trimmed.fastq \
    ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${minLen}

  # Apply HISAT2 to align reads to TAIR10 reference genome
  hisat2 -U ${reads}${id}_trimmed.fastq \
    -x ${index}TAIR10 \
    -S ${aln}${id}.sam
done