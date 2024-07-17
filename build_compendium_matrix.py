#!/usr/bin/env python3
import os
import argparse
from argparse import RawTextHelpFormatter
import json
import util
import expression

def main():
    # Set up arguments and config 
    description='''
Build an expression matrix out of rsem_genes.results files.

'''
    # Mandatory arguments
    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument("--name", dest="compendium_name", action='store', required=True,
        help="Base name of the compendium; to be used in the output filenames.")
    parser.add_argument("--input", dest="inputfile", action="store", required=True,
        help="Path to two-column TSV file listing sample ID and path to the rsem_genes.results file for that sample.")

    metric = parser.add_mutually_exclusive_group(required=True)
    metric.add_argument("--tpm", action='store_true',
        help="Produce a compendium which uses Hugo gene names and the log2(TPM+1) metric.") 
    metric.add_argument("--expected-count", action='store_true', dest="do_expected_count",
        help="Produce a compendium which uses Ensembl gene IDs and the Expected Count metric.") 

    # Optional arguments
    parser.add_argument("--output", dest="outputdir", action="store", default="output",
        help="Path to the output directory for compendium files. Will be created.")
    parser.add_argument("--mapfile", dest="ensembl_hugo_mapping_file", action="store",
        default="EnsGeneID_Hugo_Observed_Conversions.txt",
        help="Path to the file that maps between Ensemble Gene IDs and Hugo gene names.")

    arguments = vars(parser.parse_args())

    # Retrieve expression-only default config and update with user-supplied values
    config = util.default_config(arguments["compendium_name"])
    config.update(arguments)
    # If we're doing expected count, rename the output files to say ensembl_expCount instead of hugo_log2tpm
    if config["do_expected_count"]:
        expression_key="full_expression_tsv"
        # replace rightmost text only in case its in the basename
        config[expression_key] = "ensembl_expCount".join(config[expression_key].rsplit("hugo_log2tpm", 1))

    # Make outputdir or complain
    if os.path.isdir(config["outputdir"]):
        os.write(2, f"Error: Output directory '{config['outputdir']}' "
            "already exists; please delete or rename.\n".encode('utf-8'))
        exit(1)
    else:
        print("Creating new output dir '{}'".format(config["outputdir"]))
        os.mkdir(config["outputdir"])
    
    # Convert the sample paths file to json
    with open(config["inputfile"], "r") as f:
        samples_and_paths = map(lambda x: x.rstrip(), f.readlines())
        paths_real = {}
        for row in samples_and_paths:
            try:
                (sample, path) = row.split("\t")
            except ValueError as e:
                os.write(2, f"Error: Couldn't split '{row}' from input file '{config['inputfile']}'"
                "into sample ID and path.\n".encode('utf-8'))
                raise e
            paths_real[sample] = path

    if not paths_real:
        os.write(2, b"Error: Input file must list as least one sample.\n")
        exit(1)

    sample_paths_file = os.path.join(config["outputdir"], config["sample_pathsfile"])
    with open(sample_paths_file, "w") as f:
        json.dump(paths_real, f, indent=2)

    #for k,v in config.items():
    #    print( f"{k}:{v}")

    # Create the expression files
    expression.main(config)
    
if __name__ == "__main__":
    main()
