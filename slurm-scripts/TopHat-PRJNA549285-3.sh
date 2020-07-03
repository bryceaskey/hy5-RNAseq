#!/bin/bash
#SBATCH --job-name=TopHat-PRJNA549285-3       # Job name
#SBATCH --mail-type=END,FAIL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu           # Where to send mail	
#SBATCH --account=jkim6                       # Group providing CPU and memory resources
#SBATCH --qos=jkim6                           # QOS to run job on (investment or burst)
#SBATCH --ntasks=1                            # Number of CPU cores to use
#SBATCH --mem=2gb                             # Job memory request
#SBATCH --time=48:00:00                       # Time limit hrs:min:sec (max is 744:00:00)
#SBATCH --output=TopHat-PRJNA549285-3_%j.log  # Standard output and error log

pwd; hostname; date

module load trimmomatic/0.39 bowtie2/2.3.5 tophat/2.1.2

echo 'Trimming, filtering, and aligning selected RNAseq data in PRJNA549285 to TAIR10 reference genome'

ref=/ufrc/jkim6/share/braskey/data/TAIR10/
index=/ufrc/jkim6/share/braskey/data/TAIR10/bowtie2-index/
reads=/ufrc/jkim6/share/braskey/data/PRJNA549285/
aln=/ufrc/jkim6/share/braskey/data/PRJNA549285/TopHat/

# Fix chromosome naming in TAIR10.fa to match TAIR10.gff
#sed -i 's|1 CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02|Chr1|g' ${ref}TAIR10.fa
#sed -i 's|2 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02|Chr2|g' ${ref}TAIR10.fa
#sed -i 's|3 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02|Chr3|g' ${ref}TAIR10.fa
#sed -i 's|4 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02|Chr4|g' ${ref}TAIR10.fa
#sed -i 's|5 CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02|Chr5|g' ${ref}TAIR10.fa
#sed -i 's|mitochondria CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2005-06-03|ChrM|g' ${ref}TAIR10.fa
#sed -i 's|chloroplast CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2005-06-03|ChrC|g' ${ref}TAIR10.fa

# Copy reference and annotation files into index folder, and generate index
#cp ${ref}TAIR10.fa ${index}
#cp ${ref}TAIR10.gff ${index}
#bowtie2-build ${index}TAIR10.fa ${index}TAIR10

for id in SRR9313210
do
  # Apply trimmomatic to remove adapter sequences and filter low quality reads
  #trimmomatic SE -threads 1 \
  #  ${reads}${id}.fastq ${reads}${id}_trimmed.fastq \
  #  ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  # Apply TopHat to align trimmed and filtered reads to TAIR10 reference genome
  mkdir ${aln}${id}
  tophat -o ${aln}${id} \
    -G ${index}TAIR10.gff \
    --library-type fr-firststrand \
    --min-intron-length 40 --max-intron-length 2000 \
    ${index}TAIR10 \
    ${reads}${id}_trimmed.fastq
done