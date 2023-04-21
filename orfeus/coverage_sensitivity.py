#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Sensitivity simulation
# Author: Mary Richardson
# Date: 2023.03.20
# -----------------------------------------------------------------------------

import logging
import numpy as np
import multiprocessing as mp

logger = logging.getLogger(__name__)

from hmm_simulation import simulation


def find_event(path,
               event_idx,
               i,
               window):
    """
    Determine whether the specified event is present at i +/- window in the
    path
    Args:
        path: state path
        event_idx: state indices of event
        i: position in path
        window: number of nucleotides +/- i to look
    Returns:
        True or False: whether the event is present
    Raises:

    """
    L = len(path)

    # Window to look for event
    window_min = max(0, i-window)
    window_max = min(L, i+window+1)

    if any(j in path[window_min:window_max] for j in event_idx):
       return True
    else: return False


def find_orf(true_path,
             pred_path,
             idx,
             start_idx,
             stop_idx,
             window):
    """
    Determine whether the specified event is present at i +/- window in the
    path
    Args:
        true_path: simulated path
        pred_path: predicted path
        idx: state indices of ORF codons
        start_idx: state indices of ORF start codons
        stop_idx: state indices of ORF stop codons
        window: number of nucleotides +/- i to look
    Returns:
        found: whether the ORF is present within the searched window
    Raises:

    """
    # Find the true site(s) of the event
    i = list(np.where(np.in1d(true_path, idx))[0])

    # Check the predicted path for the event
    found_start = find_event(pred_path, start_idx, min(i), window)
    found_stop = find_event(pred_path, stop_idx, max(i), window)

    found = found_start and found_stop

    return found


def find_recoding(true_path,
                  pred_path,
                  idx,
                  window):
    """
    Determine whether the specified event is present at i +/- window in the
    path
    Args:
        true_path: simulated path
        pred_path: predicted path
        idx: state indices of the recoding event
        window: number of nucleotides +/- i to look
    Returns:
        found: whether the recoding event is present within the searched window
    Raises:

    """
    # Find the true site(s) of the event
    i = list(np.where(np.in1d(true_path, idx))[0])

    # Check the predicted path for the event
    found = find_event(pred_path, idx, max(i), window)

    return found


def sensitivity_search(true_path,
                       pred_path,
                       altorfs_idx,
                       window,
                       event):
    """
    Compares the true and predicted paths to determine whether the simulated
    altORF event is correctly inferred within +/- the specified window of
    nucleotides
    Args:
        true_path: simulated path
        pred_path: predicted path
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        window: number of nucleotides +/- true site to look
        event: event type to look for
            ORF = canonical ORF
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
    Returns:
        True or False: whether the event is correctly inferred within +/- the
            specified window in the predicted path
    Raises:

    """
    if event=='uORF':
        idx = altorfs_idx['uORF']
        start_idx = altorfs_idx['start_uORF']
        stop_idx = altorfs_idx['stop_uORF']
        return find_orf(true_path, pred_path, idx, start_idx, stop_idx, window)

    elif event=='dORF':
        idx = altorfs_idx['dORF']
        start_idx = altorfs_idx['start_dORF']
        stop_idx = altorfs_idx['stop_dORF']
        return find_orf(true_path, pred_path, idx, start_idx, stop_idx, window)

    elif event=='uORFdORF':
        idx = altorfs_idx['uORF'] + altorfs_idx['dORF']
        start_idx = altorfs_idx['start_uORF'] + altorfs_idx['start_dORF']
        stop_idx = altorfs_idx['stop_uORF'] + altorfs_idx['stop_dORF']
        return find_orf(true_path, pred_path, idx, start_idx, stop_idx, window)

    elif event=='pPRF':
        idx = altorfs_idx['PRF']
        return find_recoding(true_path, pred_path, idx, window)

    elif event=='mPRF':
        idx = altorfs_idx['PRF']
        return find_recoding(true_path, pred_path, idx, window)

    elif event=='SCR':
        idx = altorfs_idx['SCR']
        return find_recoding(true_path, pred_path, idx, window)

    elif event=='ORF':
        idx = altorfs_idx['start'] + altorfs_idx['sense'] + altorfs_idx['stop']
        start_idx = altorfs_idx['start']
        stop_idx = altorfs_idx['stop']
        return find_orf(true_path, pred_path, idx, start_idx, stop_idx, window)

    else:
        logger.error(('Invalid event type specified for sensitivity calculation. '
                      'Choose from: uORF, dORF, uORFdORF, pPRF, mPRF, SCR, ORF'))


def sensitivity(p_emit_nt,
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
                window,
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
        window: number of nucleotides +/- true site to look
        coverage: riboseq coverage (mean reads/nt)
        event: event type to look for
            ORF = canonical ORF
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
    Returns:
        found: 1 if the true event was found in the predicted sequence, 0 if not
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

    # Check whether the event was correctly predicted
    found = sensitivity_search(true_path,
                               pred_path,
                               altorfs_idx,
                               window,
                               event)

    return int(found)


def sensitivity_simulation(parameters,
                           nts,
                           coverages,
                           iters,
                           window,
                           event,
                           threads,
                           file):
    """
    Simulates sequences and compares the true path to the predicted path
        in parallel to estimate sensitivity of event detection
    Args:
        parameters: all model transition and emission parameters
        nts: list of valid nucleotides
        coverages: mean riboseq coverages to simulate when simulating sequences
        iters: number of sequences to simulate at each coverage value
        window: number of nucleotides +/- true site to look
        event: event type to look for
            ORF = canonical ORF
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
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
            f.write('Simulating ' + event + '\n')
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
                     window,
                     coverage,
                     event)]
            results = pool.starmap(sensitivity, args * iters)
        pool.join()
        res = np.array(results)

        sens = sum(res)/len(res) # TP/(TP+FN)
        result[coverage] = sens

        try:
            with open(file,'a') as f:
                f.write('%f coverage sensitivity: %f \n' % (coverage, sens))
        except IOError:
            logging.error('Could not open coverage ouput file %s' % file)

    return result
