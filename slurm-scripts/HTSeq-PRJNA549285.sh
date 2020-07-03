#!/bin/bash
#SBATCH --job-name=HTSeq-PRJNA549285        # Job name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --account=jkim6                     # Group providing CPU and memory resources
#SBATCH --qos=jkim6                         # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                          # Number of CPU cores to use
#SBATCH --mem=2gb                           # Job memory request
#SBATCH --time=48:00:00                     # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=HTSeq-PRJNA549285_%j.log   # Standard output and error log

pwd; hostname; date

module load samtools/1.10 python/3.8 htseq/0.11.2

echo 'Counting mapped reads from RNAseq data in PRJNA549285 aligned to TAIR10 reference genome'

aln=/ufrc/jkim6/share/braskey/data/PRJNA549285/TopHat/
gff_file=/ufrc/jkim6/share/braskey/data/TAIR10/TAIR10.gff
counts=/ufrc/jkim6/share/braskey/data/PRJNA549285/TopHat/expression-counts/

# Alignment files need to be converted from .bam to .sam format first
# Then use htseq-count to measure gene expression
for id in SRR9313209 SRR9313210 SRR9313211 SRR9313212 SRR9313213 SRR9313214 SRR9313223 SRR9313224 SRR9313225 SRR9313226 SRR9313227 SRR9313228
do
  #samtools view -h ${aln}${id}/accepted_hits.bam > ${aln}alignments/${id}.sam
  python -m HTSeq.scripts.count -m intersection-strict -s reverse --idattr=Parent \
    ${aln}alignments/${id}.sam ${gff_file} > ${counts}${id}_counts.txt
done