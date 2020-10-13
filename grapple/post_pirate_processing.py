import argparse
import csv
import os
import subprocess

def post_pir_proc():

    #1. recreate PIRATE.all_alleles.tsv
    script_1 = "scripts/link_clusters_runner.pl"
    if path.exists("/PIRATE_.all_alleles.tsv"):
        print("PIRATE.all_alleles.tsv already exists")
    else
        pipe_1 = subprocess.Popen(["perl",
                                os.path.join(args.pirate_path, script_1),
                                "-l", "./loci_list.tab",
                                "-t", args.thresholds,
                                "-o", "./",
                                "-c", "./co-ords",
                                "-parallel", args.threads,
                                "--all-alleles"
                                ])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--thresholds", type = str, required = True, help = "list of thresholds to run these scripts over. Must be subset of those run in PIRATE")
    parser.add_argument("-p", "--split_paralogs", action = store_false, default = "pangenome_maps", help = "Output directory name")
    parser.add_argument("-q", "--pirate_path", type = Path, required = True, help = "Path to PIRATE folder")
    parser.add_argument("-d", "--paralog_dir", type = Path, required = False, default = "with-paralogs/", help = "Set name of output directory")
    parser.add_argument("-n", "--threads", type = int, required = False, default = 2, help = "Number of threads to use")
    args = parser.parse_args()

    post_pir_proc(args.thresholds, args.split_paralogs, args.pirate_path, args.paralog_dir, args.threads)
