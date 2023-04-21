#!/bin/env python

# -----------------------------------------------------------------------------
# HMM run model functions
# Author: Mary Richardson
# Date: 2021.10.25
# -----------------------------------------------------------------------------

import logging
import numpy as np
import multiprocessing as mp

logger = logging.getLogger(__name__)

from hmm_helper import convert_emissions
from hmm_algorithms import viterbi, score_path

NUCLEOTIDES = ['A','C','G','T']


def get_annotated_path(path, states):
    """
    Convert a path of states to a path of state indices
    Args:
        path: path of states
        states: list of state names
    Returns:
        path: path of state indices
    Raises:

    """
    # Get the annotated state path
    try:
        path = [states.index(s) for s in path]
    except:
        path = []

    return path


def get_transcripts(data_df,
                    states,
                    coverage=0):
    """
    Generates a list of transcripts with at least the specified mean coverage
    Args:
        data_df: pandas DataFrame with all transcript info
        states: list of state names
        coverage: (optional) minimum mean coverage threhold
    Returns:
        transcripts: dict of transcripts with at least the minimum coverage
    Raises:

    """
    # Get the transcript properties of interest
    df = data_df[['transcript_id',
                  'transcript_name',
                  'seq',
                  'log_reads',
                  'state',
                  'mean_reads']]

    # Filter down to transcripts above the coverage threshold
    df = df[df['mean_reads'] > coverage]

    # Group by transcript
    df['transcript_id'] = df['transcript_id'].apply(str)
    df['transcript_name'] = df['transcript_name'].apply(str)
    df['mean_reads'] = df['mean_reads'].apply(float)

    df = df.groupby(['transcript_id',
                     'transcript_name',
                     'mean_reads'],
                      as_index=False) \
           .agg({'seq': lambda x: list(x),
                 'log_reads': lambda x: list(x),
                 'state': lambda x: get_annotated_path(list(x), states)})

    # Save transcript
    transcripts = df.values.tolist()

    return transcripts


def hmm_single(x, r,
               path,
               parameters_h1,
               parameters_h0):
    """
    Run the mdoel on a single transcript
    Args:
        x: nucleotide sequence
        r: rho sequence
        path: annotated path
        parameters_h1: list of all parameter matrices for the current model
        parameters_h0: list of all parameter matrices for the null model
    Returns:
        result: list of results of running the model on this transcript
    Raises:

    """
    # Unpack the arguments
    states_h1, states_idx_h1, altorfs_idx_h1, \
    p_emit_nt_h1, p_emit_rho_h1, rho_ranges_h1, p_trans_h1, p_trans_rev_h1, \
    prev_states_options_h1, next_states_options_h1 = parameters_h1

    states_h0, states_idx_h0, altorfs_idx_h0, \
    p_emit_nt_h0, p_emit_rho_h0, rho_ranges_h0, p_trans_h0, p_trans_rev_h0, \
    prev_states_options_h0, next_states_options_h0 = parameters_h0

    log_emit_nt_h1 = np.log(p_emit_nt_h1)
    log_emit_rho_h1 = np.log(p_emit_rho_h1)
    log_emit_nt_h0 = np.log(p_emit_nt_h0)
    log_emit_rho_h0 = np.log(p_emit_rho_h0)
    log_trans_h1 = np.log(p_trans_h1)
    log_trans_h0 = np.log(p_trans_h0)

    # Generate the one-hot encoding and reformatted emission probabilities
    # for this sequence
    log_emit_h1 = convert_emissions(x,
                                    r,
                                    log_emit_nt_h1,
                                    log_emit_rho_h1,
                                    NUCLEOTIDES,
                                    rho_ranges_h1)
    log_emit_h0 = convert_emissions(x,
                                    r,
                                    log_emit_nt_h0,
                                    log_emit_rho_h0,
                                    NUCLEOTIDES,
                                    rho_ranges_h0)

    # Predict the viterbi state path
    path_h1, P_h1 = viterbi(log_emit_h1,
                            log_trans_h1,
                            prev_states_options_h1,
                            next_states_options_h1)
    path_h0, P_h0 = viterbi(log_emit_h0,
                            log_trans_h0,
                            prev_states_options_h0,
                            next_states_options_h0)

    # Calculate the probability of the annotated (h2), predicted (h1), and
    # null (h0) paths
    path_score_h2 = score_path(path, log_emit_h1, log_trans_h1)
    path_score_h1 = score_path(path_h1, log_emit_h1, log_trans_h1)
    path_score_h0 = score_path(path_h0, log_emit_h1, log_trans_h1)

    # Calculate the log-odds score
    log_odds_score = path_score_h1 - path_score_h0

    # Make a list of all values to return for this transcript
    return path, path_h1, path_h0, \
           path_score_h2, path_score_h1, path_score_h0, log_odds_score, \
           log_emit_h1, log_trans_h1


def hmm_parallel(args):
    """
    Run the mdoel on a single transcript
    Args: (all packaged into the list args)
        transcripts: dict of transcripts
        parameters_h1: list of all parameter matrices for the current model
        parameters_h0: list of all parameter matrices for the null model
    Returns:
        result: list of results of running the model on this transcript
    Raises:

    """
    # Unpack the arguments
    transcript, parameters_h1, parameters_h0 = args

    transcript_id, transcript_name, coverage, x, r, path = transcript

    path_h2, path_h1, path_h0, \
    path_score_h2, path_score_h1, path_score_h0, log_odds_score, \
    log_emit_h1, log_trans_h1 = hmm_single(x, r,
                                           path,
                                           parameters_h1,
                                           parameters_h0)

    # Make a list of all values to return for this transcript
    result = [transcript_id, transcript_name, coverage, x, r, \
              path_h2, path_h1, path_h0, \
              path_score_h2, path_score_h1, path_score_h0, \
              log_odds_score]

    logging.info('%s (%s)' % (transcript_id,transcript_name))

    return result


def run_hmm(data_df,
            coverage,
            parameters_h1,
            parameters_h0,
            threads):
    """
    Run the model on multiple transcripts in parallel
    Args:
        data_df: processed data
        coverage: mean reads per transcript coverage threshold
        parameters_h1: list of all parameter matrices for the current model
        parameters_h0: list of all parameter matrices for the null model
    Returns:
        results: array containing all results of running the model
    Raises:

    """
    # Unpack the parameters
    states, states_idx, altorfs_idx, \
    p_emit_nt, p_emit_rho, rho_ranges, p_trans, p_trans_rev, \
    prev_states_options, next_states_options = parameters_h1

    # Get transcripts
    transcripts = get_transcripts(data_df,
                                  states,
                                  coverage)

    with mp.Pool(threads) as pool:

        # Number of sequences
        C = len(transcripts)

        # For each sequence, calculate the predicted viterbi path
        args = [(transcripts[c], parameters_h1, parameters_h0) \
                 for c in range(C)]
        results_paths = pool.map(hmm_parallel, args)

    pool.join()

    # Get results as an array of prediction dictionaries
    results = np.array(results_paths, dtype=object)

    return results
