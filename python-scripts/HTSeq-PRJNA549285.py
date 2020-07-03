# Script to modify .gff annotation file so that each exon has a gene_id attribute.
# If running local (i.e. not on hipergator), gffutils must first be installed with
# command: py -3 -m pip install gffutils

import gffutils

gff_file = "C:/Users/Bryce/Documents/hy5-RNAseq/data/TAIR10.gff"
database = "C:/Users/Bryce/Documents/hy5-RNAseq/data/TAIR10.db"

db = gffutils.create_db(gff_file, dbfn=database, force=True, keep_order=True)
#db = gffutils.FeatureDB(database)
exons = db.children('AT1G01010', featuretype='exon')

print(exons)