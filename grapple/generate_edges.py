import argparse
import csv
import subprocess
from os import path

def generate_edges (input, output_dir, path, threads):
    
    if not path.exists(args.input):
        print("ERROR: Input file not found")

    # Set 
    out_file=os.path.splitext(args.input)[0]
    out_file_tsv=out_file + ".allelesAsGeneFamily.tsv"
    out_file_edges=out_file + ".allelesAsGeneFamily.edges"

    # 1. Generate alleles as gene family files
    with open(out_file_tsv, 'w') as out_file:
        writer = csv.writer(out_file, delimiter = '\t')
        with open(args.input) as in_file:
            reader = csv.reader(in_file, delimiter = '\t')
            header = next(reader)
            writer.writerow(header)
            if header != None:
                for row in reader:
                    new_row = [] 
                    new_row.append(row[1])
                    new_row.append(row[1]) #I promise this isn't a typo!
                    new_row.extend(row[2:])
                    writer.writerow(new_row) 

    # 2. Generate edge file using PIRATE adapter script
    script = "scripts\pangenome_graph.pl"
    output = os.path.join(output_dir, out_file_edges)
    pipe = subprocess.POpen(["perl",
                            os.path.join(args.path, script),
                            "--input", out_file_tsv,
                            "--output", output,
                            "--gffs", "./modified_gffs/",
                            "--no-clutser",
                            "--gfa1"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, required = True, help = "Input file - e.g. PIRATE.all_alleles.wp.90.tsv")
    parser.add_argument("-o", "--output_dir", type = str, required = False, default = "pangenome_maps", help = "Output directory name")
    parser.add_argument("-p", "--path", type = str, required = True, help = "Path to PIRATE folder")
    parser.add_argument("-n", "--threads", type = int, required = False, default = 2, help = "Number of threads to use")
    args = parser.parse_args()

    generate_edges(args.input, args.output, args.path, args.threads)