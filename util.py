#!/usr/bin/env python3
import os
import datetime

# Helper functions

# default_config -- returns the base json configuration
# Takes the compendium name to incorporate into output files
def default_config(name):
    defaults= {
        "ensembl_hugo_NA_key" : "NA",
        "sample_pathsfile" : "samples_vs_rgr_files_{}_{}.json".format(name, today()),
        "full_expression_tsv" : "{}_hugo_log2tpm_{}.tsv.gz".format(name, today()),
    }    
    return defaults

# Todays date in the format YYYY-MM-DD
def today():
    return str(datetime.date.today())
