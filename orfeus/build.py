#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Build ORFeus model
# Author: Mary Richardson
# Date: 2023.03.03
# -----------------------------------------------------------------------------

import itertools as it
import argparse
import logging
import os
import numpy as np

from import_data import process_data, read_data_file
from parameters_set import set_parameters
from coverage_threshold import coverage_threshold_simulation
from output_logo import logo

NUCLEOTIDES = ['A','C','G','T']
CODONS = [''.join(i) for i in it.product(NUCLEOTIDES, repeat = 3)]


def main():

# ------------------------------------------------------------------------------
# Accept command line arguments
# ------------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description='Build ORFeus model')

    logo()

    # Input data
    parser.add_argument('plus_file', type=str,
                        help='riboseq plus strand bg/wig file path')
    parser.add_argument('minus_file', type=str,
                        help='riboseq minus strand bg/wig file path')
    parser.add_argument('seqs_file', type=str,
                        help='genome or transcriptome fasta file path')
    parser.add_argument('annotations_file', type=str,
                        help='annotations gtf/gff file path')
    parser.add_argument('-s', '--skip', nargs='+', default=[],
                        help=('list known frameshifted genes that should be '
                        'filtered out before codon assignments'))

    # ORFeus options
    parser.add_argument('--data_file', type=str,
                        help='processed data file path')
    parser.add_argument('-o', '--outdir', type=str, default=os.getcwd(),
                        help=('output directory for all intermediate and '
                        'results files'))

    # Transition probabilities
    parser.add_argument('-a', '--alpha', type=float, default=1e-5,
                        help='probability of PRF')
    parser.add_argument('-b', '--beta', type=float, default=0.5,
                        help='probability of -1 PRF given PRF')
    parser.add_argument('-g', '--gamma', type=float, default=1e-4,
                        help='probability of stop codon readthrough')
    parser.add_argument('-d', '--delta', type=float, default=1e-3,
                        help='probability of upstream or downstream short ORFs')
    parser.add_argument('-z', '--zeta', type=float, default=1e-10,
                        help='probability of multiple non-overlapping ORFs')

    # Feature frequencies
    parser.add_argument('--f5', type=float, default=None,
                        help=('fraction of transcripts with 5\'UTRs '
                        '(default: calculated from annotations)'))
    parser.add_argument('--f3', type=float, default=None,
                        help=('fraction of transcripts with 3\'UTRs '
                        '(default: calculated from annotations)'))

    # Feature lengths
    parser.add_argument('--utr5', type=int, default=None,
                        help=('mean nucleotide length of 5\'UTRs '
                        '(default: calculated from annotations)'))
    parser.add_argument('--utr3', type=int, default=None,
                        help=('mean nucleotide length of 3\'UTRs '
                        '(default: calculated from annotations)'))
    parser.add_argument('--orf', type=int, default=None,
                        help=('mean nucleotide length of main ORFs '
                        '(default: calculated from annotations)'))
    parser.add_argument('--uorf', type=int, default=50,
                        help=('mean nucleotide length of upstream short ORFs '
                        '(suggested: 50)'))
    parser.add_argument('--dorf', type=int, default=50,
                        help=('mean nucleotide length of downstream short ORFs '
                        '(suggested: 50)'))

    # Feature values
    parser.add_argument('--start_codons', type=str, nargs='*',
                        default=[],
                        help=('list of valid start codons '
                        '(default: infer from the annotations)'))
    parser.add_argument('--sense_codons', type=str, nargs='*',
                        default=list(filter(lambda i: i not in ['TAA','TAG','TGA'], CODONS)),
                        help=('list of valid sense codons '
                        '(default: all codons except TAA, TAG, TGA)'))
    parser.add_argument('--stop_codons', type=str, nargs='*',
                        default=['TAA','TAG','TGA'],
                        help=('list of valid stop codons '
                        '(default: TAA, TAG, TGA)'))

    # Emission probabilities
    parser.add_argument('-r', '--bins', type=int, default=25,
                        help=('number of bins to group emission probabilities '
                        '(suggested: 25)'))
    parser.add_argument('-l', '--log', default=False, action='store_true',
                        help=('whether to bin the riboseq values in logspace '
                        '(suggested: false)'))
    parser.add_argument('-m', '--min', type=int, default=-1,
                        help=('minimum reads per transcript '
                        '(default: -1)'))
    parser.add_argument('-M', '--max', type=int, default=np.inf,
                        help=('maximum reads per transcript '
                        '(default: infinity)'))
    parser.add_argument('-p', '--pool', default=False, action='store_true',
                        help=('whether to pool the observed riboseq values '
                        'from all codons of a given state type '
                        '(suggested: false, unless you have a small '
                        'transcriptome with <=20 ORFs)'))
    parser.add_argument('-f', '--fit', default=False, action='store_true',
                        help=('whether to fit a log-normal distribution to the '
                        'observed riboseq values or use the raw frequencies '
                        '(suggested: false, unless you have a small '
                        'transcriptome with <=20 ORFs)'))

    # Simulation parameters
    parser.add_argument('-c', '--coverage', default=False, action='store_true',
                        help='whether to run the coverage threshold simulation')
    parser.add_argument('--iters', type=int, default=10,
                        help='number of iterations to simulate')
    parser.add_argument('--window', type=int, default=10,
                        help=('window (in nts) around true event to count as '
                              'correct altORF prediction'))
    parser.add_argument('--threads', type=int, default=16,
                        help=(('number of processes to run in parallel when '
                            'simulating transcripts')))

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

    # Set up logging
    #log_file = os.path.join(outdir, 'orfeus_build.log')
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        filemode = 'w',
                        level = logging.DEBUG)

# ------------------------------------------------------------------------------
# Process data
# ------------------------------------------------------------------------------

    if args.data_file:
        logging.info('Importing the processed data...')

        # Import the data (specify dtypes where necessary)
        data_df = read_data_file(args.data_file)

        logging.info('Finished importing the data.')

    else:
        logging.info('Processing the data...')

        data_file = os.path.join(outdir, 'data.txt.gz')

        # Check if file already exists
        if os.path.exists(data_file):
            logging.error(('Processed data file already exists. '
                        'If you\'d like to import the data from this file, '
                        'use the --data_file option. '
                        'Otherwise, provide a new output directory path.'))

        # If the file does not exist, create it
        else:

            # Process the data
            data_df = process_data(args.plus_file,
                                   args.minus_file,
                                   args.seqs_file,
                                   args.annotations_file,
                                   skip=args.skip)

            # Export to a tsv file
            data_df.to_csv(data_file, sep='\t', compression='gzip')

            logging.info('Finished processing the data.')

# ------------------------------------------------------------------------------
# Set parameters
# ------------------------------------------------------------------------------

    logging.info('Setting the model parameters...')

    # Set parameters for H1 model (altORFs and canonical ORFs)
    parameters_h1 = set_parameters(data_df,
                                   args.f5,
                                   args.f3,
                                   args.utr5,
                                   args.utr3,
                                   args.orf,
                                   args.uorf,
                                   args.dorf,
                                   args.start_codons,
                                   args.sense_codons,
                                   args.stop_codons,
                                   args.bins,
                                   args.log,
                                   args.min,
                                   args.max,
                                   args.pool,
                                   args.fit,
                    [args.alpha, args.beta, args.gamma, args.delta, args.zeta])

    logging.info('Setting the null model parameters...')

    # Set parameters for H0 model (canonical ORFs only)
    parameters_h0 = set_parameters(data_df,
                                   args.f5,
                                   args.f3,
                                   args.utr5,
                                   args.utr3,
                                   args.orf,
                                   args.uorf,
                                   args.dorf,
                                   args.start_codons,
                                   args.sense_codons,
                                   args.stop_codons,
                                   args.bins,
                                   args.log,
                                   args.min,
                                   args.max,
                                   args.pool,
                                   args.fit,
                                  [0,0,0,0,0])

    # Save the parameters
    np.save(os.path.join(outdir, 'parameters_h1'),
                parameters_h1)
    np.save(os.path.join(outdir, 'parameters_h0'),
                parameters_h0)

    logging.info('Finished setting parameters')

# ------------------------------------------------------------------------------
# Coverage threshold
# ------------------------------------------------------------------------------

    if args.coverage:

        logging.info('Simulating sequences...')

        coverage_outdir = os.path.join(outdir, 'coverage')
        os.makedirs(coverage_outdir, exist_ok=True)

        # Set coverage thresholds for testing
        coverages = np.arange(0.01, 1, 0.01)

        coverage_threshold_simulation(data_df,
                                      parameters_h1,
                                      coverages,
                                      args.iters,
                                      args.window,
                                      coverage_outdir,
                                      args.threads)

        logging.info('Finished simulating sequences')




if __name__== "__main__":
    main()
