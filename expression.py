#!/usr/bin/env python3

# Expression.py
# Create expression files for a compendium
# Takes a long time to run
#
# Created from v4.3/create-expression-files.py
# Which was derived from compendium_creation/v4/expression/v4-TCGA_and_TARGET.ipynb

# Outputs:
# - cohort.tsv.gz                   all samples, all genes, 12 significant digits
#                                uniq hugo, log2 TPM + 1 (or, Ensembl and Expected Count)

import numpy as np
import pandas as pd
import time
import math
import json
import csv
import sys
import os

# Cheap version of timeit with a message
def twrite(text):
    sys.stdout.write(text)
    sys.stdout.flush()
    return time.time()
def tend(start):
    end = time.time()
    print( " Elapsed: {} s".format(round(end - start, 1)))

# Returns dictionary { sample_id : path}
def get_sample_paths(sp_file):
    start = twrite("Getting paths to samples...")
    with open(sp_file, "r") as f:
        sp = json.load(f)
    tend(start)
    return sp

# Get a series of the RSEM expression of a single sample, either TPM or expected_count
def get_single_expression(expression_file, column):
    expression_column = pd.read_csv(expression_file, sep="\t", index_col=0)[column]
    return expression_column

# Takes the dict of sample id -> path
# Creates a single DF of the expression
# column headers are unique sample IDs; rows are ensembl;
# values are un-modified from the rsem_genes_results.file, either TPM or expected_counts
def make_df_from_individual(samplepaths, metric):
    start = twrite("Getting {} expression from individual files...".format(metric))
    ensembl_df = pd.DataFrame()
    c = 0
    for sample_id, path in samplepaths.items():
        ensembl_df[sample_id] = get_single_expression(path, metric)
        c += 1
        sys.stdout.write(" {}".format(c))
        sys.stdout.flush()
        #if c > 10:
        #    break
    tend(start)
    return ensembl_df

# Takes ensembl df in TPM space
# returns hugo df with unique gene ids in TPM space
# Duplicate gene IDs have values summed rather than averaged.
def convert_to_hugo(ensembl_df, ensembl_hugo_mapping_file, ensembl_hugo_NA_key):
    start = twrite("Converting to Hugo...")
    ensembl_df_local = ensembl_df.copy() # Don't mutate original DF!
    with open(ensembl_hugo_mapping_file, mode='r') as infile:
        ensembl_to_hugo = dict((rows[1],rows[0]) for rows in csv.reader(infile, delimiter='\t'))
    ensembl_df_local.index = ensembl_df_local.index.map(lambda x: ensembl_to_hugo[x])
    # Drop the NA key. if it wasn't present, then hugo_df is same as this modified ensembl df.
    try:
        hugo_df = ensembl_df_local.drop(ensembl_hugo_NA_key)
    except KeyError:
        print("Didn't find any {} in index to drop--please confirm this is expected.".format(ensembl_hugo_NA_key))
        hugo_df = ensembl_df_local

    tend(start)
    start = twrite("Uniquifying gene names...")
    uniq_hugo_df = hugo_df.groupby(hugo_df.index).sum()
    tend(start)
    return uniq_hugo_df

def log2_normalize(hugo_df):
    start = twrite("Log2 +1 normalizing...")
    norm_df = np.log2(hugo_df + 1)
    tend(start)
    return norm_df

######## Main ########
def main(config):

    workdir=config["outputdir"]
    sample_paths_file= os.path.join(workdir, config["sample_pathsfile"])
    ensembl_hugo_mapping_file= config["ensembl_hugo_mapping_file"]
    ensembl_hugo_NA_key= config["ensembl_hugo_NA_key"]

    # Outputs
    full_tsv=os.path.join(workdir, config["full_expression_tsv"])

    print( "Running..." )
    # Get sample lists and figure out which TH and TCGA samples we want

    sample_paths = get_sample_paths(sample_paths_file)
    individual_samples = { sid : sample_paths[sid]  for sid in sample_paths.keys() if sample_paths[sid] }

    print( "Will load {} individual samples.".format(
        len(individual_samples)))

    # Next, assemble the df of samples retrieved from individual rsem_genes.results files
    # unique sample IDs, ensembl gene IDs, & either plain TPM or plain expected_counts
    if config["do_expected_count"]:
        # expected_counts, ensembl
        individual_ensembl_plain_df = make_df_from_individual(individual_samples, metric="expected_count")
    elif config["tpm"]:
        # TPM, ensembl
        individual_ensembl_plain_df = make_df_from_individual(individual_samples, metric="TPM")
    else:
        print("Unknown metric! Please choose Expected Count or TPM.")
        exit()

    # If we have any individual samples, process them  & merge with bulk
    if(individual_ensembl_plain_df.empty):
        print("No individual sample expression files were requested.")
        exit()
    else:
        if(config["do_expected_count"]):
            # Expected counts stays in plain ensembl, so just merge
            combined_df = individual_ensembl_plain_df
        elif config["tpm"]:
            # Convert to hugo while in TPM space.
            individual_hugo_tpm_df = convert_to_hugo(individual_ensembl_plain_df, ensembl_hugo_mapping_file, ensembl_hugo_NA_key)

            # Put into log2(TPM+1) space
            individual_hugo_log2_tpm_df = log2_normalize(individual_hugo_tpm_df)

            combined_df = individual_hugo_log2_tpm_df

    # So now we have the combined_df which is either plain expected counts ensembl or log2(tpm + 1) hugo
    # Write to TSV in 12 digits
    start = twrite("Writing to tsv.gz...")
    combined_df.to_csv(full_tsv,
                       sep="\t",
                       compression='gzip',
                       header=True,
                       index_label="Gene",
                       float_format="%.12g")
    tend(start)

    print("Done.")

if __name__ == "__main__":
    print("Please import this script and pass it your config dict:\n  import expression as exp\n  exp.main(config)")
