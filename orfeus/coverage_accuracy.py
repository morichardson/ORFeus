#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Accuracy simulation
# Author: Mary Richardson
# Date: 2023.03.20
# -----------------------------------------------------------------------------

import logging
import numpy as np
import multiprocessing as mp

logger = logging.getLogger(__name__)

from hmm_simulation import simulation


def accuracy_calculation(true_path,
                         pred_path):
    """
    Compares the true and predicted paths to determine the accuracy of the
        prediction
    Args:
        true_path: simulated path
        pred_path: predicted path
    Returns:
        accuracy: mean fraction of states that are correctly assigned in the
            predicted path
    Raises:

    """
    T=0
    # Check how many positions are correctly predicted
    for i in range(len(true_path)):
        if true_path[i] == pred_path[i]:
            T += 1
    accuracy = T/len(true_path)

    return accuracy


def accuracy(p_emit_nt,
             p_emit_rho,
             p_trans,
             p_trans_rev,
             log_emit_nt,
             log_emit_rho,
             log_trans,
             prev_states_options,
             next_states_options,
             nts,
             rho_ranges,
             altorfs_idx,
             coverage,
             event):
    """
    Simulates a path and corresponding sequence of nts and rho values,
    then runs Viterbi on the sequence to predict the path
    Args:
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        p_trans_rev: the reverse transition probability matrix
        log_emit_nt: the nucleotide log emission probability matrix
        log_emit_rho: the riboseq log emission probability matrix
        log_trans: the log transition probability matrix
        prev_state_options: list of indices of all possible
            (nonzero probability) previous states
        next_state_options: list of indices of all possible
            (nonzero probability) next states
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        coverage: riboseq coverage (mean reads/nt)
    Returns:
        accuracy: mean fraction of states that are correctly assigned in the
            predicted path
    Raises:

    """
    # Simulate a sequence
    true_path, pred_path = simulation(p_emit_nt,
                                      p_emit_rho,
                                      p_trans,
                                      p_trans_rev,
                                      log_emit_nt,
                                      log_emit_rho,
                                      log_trans,
                                      prev_states_options,
                                      next_states_options,
                                      nts,
                                      rho_ranges,
                                      altorfs_idx,
                                      coverage,
                                      event)

    # Calculate the accuracy of the prediction
    accuracy = accuracy_calculation(true_path, pred_path)

    return accuracy


def accuracy_simulation(parameters,
                        nts,
                        coverages,
                        iters,
                        event,
                        threads,
                        file):
    """
    Simulates sequences and compares the true path to the predicted path
        in parallel to estimate specificity of event detection
    Args:
        parameters: all model transition and emission parameters
        nts: list of valid nucleotides
        coverages: mean riboseq coverages to simulate when simulating sequences
        iters: number of sequences to simulate at each coverage value
        threads: number of processes to run in parallel (when simulating
            sequences)
        file: output file for results
    Returns:

    Raises:

    """
    states, states_idx, altorfs_idx, \
    p_emit_nt, p_emit_rho, rho_ranges, p_trans, p_trans_rev, \
    prev_states_options, next_states_options = parameters

    log_emit_nt = np.log(p_emit_nt)
    log_emit_rho = np.log(p_emit_rho)
    log_trans = np.log(p_trans)

    try:
        with open(file,'w') as f:
            f.write('Simulating any ORF or altORF \n')
    except IOError:
        logging.error('Could not open coverage ouput file %s' % file)

    result = {}

    for coverage in coverages:
        # Simulate sequences in parallel
        with mp.Pool(processes=threads) as pool:
            args = [(p_emit_nt,
                     p_emit_rho,
                     p_trans,
                     p_trans_rev,
                     log_emit_nt,
                     log_emit_rho,
                     log_trans,
                     prev_states_options,
                     next_states_options,
                     nts,
                     rho_ranges,
                     altorfs_idx,
                     coverage,
                     event)]
            results = pool.starmap(accuracy, args * iters)
        pool.join()
        res = np.array(results)

        acc = np.mean(res)
        result[coverage] = acc

        try:
            with open(file,'a') as f:
                f.write('%f coverage accuracy: %f \n' % (coverage, acc))
        except IOError:
            logging.error('Could not open coverage ouput file %s' % file)

    return result
