#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Build ORFeus model
# Author: Mary Richardson
# Date: 2023.03.03
# -----------------------------------------------------------------------------

import argparse
import logging
import os
import numpy as np

from import_data import read_data_file
from hmm_run import run_hmm
from output_text import write_results_fasta, write_results_table
from output_logo import logo


def main():

# ------------------------------------------------------------------------------
# Accept command line arguments
# ------------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description='Run ORFeus model')

    logo()

    # Input data
    parser.add_argument('data_file', type=str,
                        help=('processed data file path '
                            '(generated by orfeus_build_model.py)'))
    parser.add_argument('parameters_file_h1', type=str,
                        help=('model parameters file path '
                            '(generated by orfeus_build_model.py)'))
    parser.add_argument('parameters_file_h0', type=str,
                        help=('null model parameters file path '
                            '(generated by orfeus_build_model.py)'))

    # ORFeus options
    parser.add_argument('-o', '--outdir', type=str, default=os.getcwd(),
                        help=('output directory for all intermediate and '
                        'results files'))
    parser.add_argument('-c', '--coverage', type=float, default=0.1,
                        help=('mean riboseq coverage threshold'))
    parser.add_argument('--threads', type=int, default=8,
                        help=(('number of processes to run in parallel when '
                            'running the model on all transcripts')))

    args = parser.parse_args()

# ------------------------------------------------------------------------------
# Prepare output directory
# ------------------------------------------------------------------------------

    # Create output directory
    if args.outdir == os.getcwd():
        outdir = os.path.join(args.outdir, r'orfeus')
        os.makedirs(outdir, exist_ok=True)
    else:
        if not os.path.isdir(args.outdir):
            print('Output directory does not exist. Please provide a valid path.')
        outdir = args.outdir

    # Set up log file
    #log_file = os.path.join(outdir, 'orfeus_build.log')
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        filemode = 'w',
                        level = logging.DEBUG)

# ------------------------------------------------------------------------------
# Import data and parameters
# ------------------------------------------------------------------------------

    logging.info('Importing the processed data...')

    # Import the data (specify dtypes where necessary)
    data_df = read_data_file(args.data_file)

    logging.info('Finished importing the data.')


    logging.info('Importing the parameters...')

    parameters_h1 = np.load(args.parameters_file_h1, allow_pickle=True)
    parameters_h0 = np.load(args.parameters_file_h0, allow_pickle=True)

    logging.info('Finished importing the parameters.')

# ------------------------------------------------------------------------------
# Run model
# ------------------------------------------------------------------------------

    logging.info('Running the model...')

    # Calculate the predicted viterbi path for each sequence
    results = run_hmm(data_df,
                      args.coverage,
                      parameters_h1,
                      parameters_h0,
                      args.threads)

    # Write results to files
    write_results_fasta(results,
                        parameters_h1,
                        os.path.join(outdir, 'predicted_mRNA.fn'),
                        os.path.join(outdir, 'predicted_protein.fa'))

    write_results_table(results,
                        parameters_h1,
                        os.path.join(outdir, 'results.txt'))

    logging.info('Finished running the model.')


if __name__== "__main__":
    main()
