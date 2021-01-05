import argparse
from Bio import Seq
from BCBio import GFF
import pandas as pd

def add_locus_tag_dir(gff, layout, seq_prod):

    # Dictionary of direction of locus tags
    limit_info = dict(gff_type = seq_prod)
    dct = {}
    with open(gff) as gff_file:
        for rec in GFF.parse(gff_file, limit_info=limit_info):
            for f in rec.features:
                locus_tag = f.qualifiers["locus_tag"][0]
                strand = f.qualifiers["strand"][0]
                dct[tag.strip('\n')] = {locus_tag: strand}
    
    # Write each to layout file
    gff_name = gff.strip(".gff")
    with open(layout, 'r+') as out_file:
        out_file_reader = csv.reader(out_file, sep = '\t')
            for row in out_file_reader:
                if row[0] == "//NODECLASS":
                    if row[1] in dct:
                        dct.get(row[1])
                        print('//NODECLASS\t\"{row[0]}\"\t\"{row[1]}\"\t\"gff_name\"')     


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gff", type = str, help = "input GFF3 file")
    parser.add_arguemnt("-l", "--layout", type = str, help = ".layout file to be written to")
    parser.add_argument("-q", "--seq-prod", type = str, required = False, default = "CDS", help = "Seq product type - e.g. cds, rRNA, tRNA. Defaults to CDS")
    args = parser.parse_args()