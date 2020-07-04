# Script to modify .gff annotation file so that each exon has a gene_id attribute that matches 
# that of the parent gene.
# If not installed, gffutils must first be installed with command: py -3 -m pip install gffutils
import gffutils
import os.path

# Set paths to input gff file, output modified gff file that will be created, and 
# gffutils database (if it exists).
gff_file = "C:/Users/bca08_000/Documents/hy5-RNAseq/data/TAIR10.gff"
modified_gff_file = "C:/Users/bca08_000/Documents/hy5-RNAseq/data/TAIR10_modified.gff"
database = "C:/Users/bca08_000/Documents/hy5-RNAseq/data/TAIR10.db"

# If gffutils database exists, use it. Otherwise make a new database from the input gff file.
if path.exists(database):
    db = gffutils.FeatureDB(database)
else:
    db = gffutils.create_db(gff_file, dbfn=database, force=True, keep_order=True)

# Loop through each feature in the database object, and write to a new gff file an exact copy of
# any feature which is not an exon. For features which are exons, write a modified copy of the
# feature which has a new attribute, "gene_id", with value matching the ID of its parent gene.
with open(modified_gff_file, "w") as fout:
    for d in db.directives:
        fout.write('##{0}\n'.format(d))
    for feature in db.all_features():
        if feature.featuretype != "exon":
            fout.write(str(feature) + "\n")
            if feature.featuretype == "gene":
                for exon in db.children(feature.id, featuretype="exon", order_by="start"):
                    exon.attributes["gene_id"] = feature.id
                    fout.write(str(exon) + "\n")