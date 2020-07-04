#!/bin/bash
#SBATCH --job-name=featureCounts-PRJNA549285        # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                 # Where to send mail	
#SBATCH --account=jkim6                             # Group providing CPU and memory resources
#SBATCH --qos=jkim6                                 # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                                  # Number of CPU cores to use
#SBATCH --mem=2gb                                   # Job memory request
#SBATCH --time=48:00:00                             # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=featureCounts-PRJNA549285_%j.log   # Standard output and error log

pwd; hostname; date

module load subread/2.0.0 samtools/1.10

echo 'Counting mapped reads from RNAseq data in PRJNA549285 aligned to TAIR10 reference genome'
TopHat=/ufrc/jkim6/share/braskey/data/PRJNA549285/TopHat/
aln=/ufrc/jkim6/share/braskey/data/PRJNA549285/TopHat/alignments/
gff_file=/ufrc/jkim6/share/braskey/data/TAIR10/TAIR10_modified.gff
counts=/ufrc/jkim6/share/braskey/data/PRJNA549285/TopHat/featureCounts/

# Use featureCounts to measure gene expression
# SRR9313209 SRR9313210 SRR9313211 SRR9313212 SRR9313213 SRR9313214 SRR9313223 SRR9313224 SRR9313225 SRR9313226 SRR9313227 SRR9313228
for id in SRR9313213
do
  samtools view -h ${TopHat}${id}/accepted_hits.bam > ${aln}${id}.sam
  featureCounts -T 1 -t exon -g gene_id -s 2 \
    -a ${gff_file} \
    -o ${counts}${id}_counts.txt \
    ${aln}${id}.sam
done

# Extract gene expression counts from output files
#mkdir ${counts}just-counts
# SRR9313209 SRR9313210 SRR9313211 SRR9313212 SRR9313213 SRR9313214 SRR9313223 SRR9313224 SRR9313225 SRR9313226 SRR9313227 SRR9313228
for id in SRR9313213
do
  cut -f 1,7 ${counts}${id}_counts.txt > ${counts}just-counts/${id}_counts.txt
  tail -n +2 ${counts}just-counts/${id}_counts.txt > ${counts}just-counts/${id}_counts.txt.tmp && mv ${counts}just-counts/${id}_counts.txt.tmp ${counts}just-counts/${id}_counts.txt
done

# Extract gene lengths - needed to calculate RPKM, FPKM, and TPM
cut -f 1,6 ${counts}SRR9313209_counts.txt > ${counts}gene-lengths.txt
tail -n +2 ${counts}gene-lengths.txt > ${counts}gene-lengths.txt.tmp && mv ${counts}gene-lengths.txt.tmp ${counts}gene-lengths.txt