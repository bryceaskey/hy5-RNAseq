# Script to modify .gff annotation file so that each exon has a gene_id attribute that matches 
# that of the parent gene.
# If not installed, gffutils must first be installed with command: py -3 -m pip install gffutils
import gffutils
import os.path

# Set paths to input gff file, output filtered gff file that will be created, and 
# gffutils database (if it exists).
gff_file = "D:/bioinformatics/TAIR10.gff"
modified_gff_file = "D:/bioinformatics/TAIR10_filtered.gff"
database = "D:/bioinformatics/TAIR10.db"

# If gffutils database exists, use it. Otherwise make a new database from the input gff file.
if os.path.exists(database):
    db = gffutils.FeatureDB(database)
else:
    db = gffutils.create_db(gff_file, dbfn=database, force=True, keep_order=True)

# Loop through each feature in the database object, and write to a new gff file only those features
# which are genes. Any mRNA/exon/other objects should not be included.
with open(modified_gff_file, "w") as fout:
    for d in db.directives:
        fout.write('##{0}\n'.format(d))
    for feature in db.all_features():
        if feature.featuretype == "gene":
            fout.write(str(feature) + "\n")