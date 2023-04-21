#!/bin/env python

# -----------------------------------------------------------------------------
# HMM simulation functions
# Author: Mary Richardson
# Date: 2020.01.08
# -----------------------------------------------------------------------------

import logging
import numpy as np

logger = logging.getLogger(__name__)

from hmm_helper import convert_emissions
from hmm_algorithms import viterbi


def sample_transitions(p_trans):
    """
    Sample a path from the model, using the model transitions
    Args:
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
    Returns:
        path: simulated path
    Raises:

    """
    M = len(p_trans) # Number of states

    # Choose a state path using transition probabilities
    path = []

    # Start at the begin state and correct for the offset (1-indexed)
    next_state = np.random.choice(M, p=p_trans[0]) - 1

    # Continue adding next states until you reach the end state
    while next_state != -1:
        curr_state = next_state
        path.append(curr_state)
        next_state = np.random.choice(M, p=p_trans[curr_state+1]) - 1

    return path


def sample_emissions(path,
                 p_emit_nt,
                 p_emit_rho,
                 nts,
                 rho_ranges):
    """
    Sample a sequence for the path from the model, using the model emissions
    Args:
        path: simulated path
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
    Returns:
        x: simulated nt sequence
        r: simulated rho sequence
    Raises:

    """
    # Set possible emission values (zero if zero bin, mean of bin otherwise)
    rhos = np.append([0], np.mean([rho_ranges[1:-1],
                                   rho_ranges[2:]],
                                   axis=0))

    # Choose the nucleotide and riboseq observation for each state using
    # emission probabilities
    x = []
    r = []
    for s in path:

        # Assign the nt value
        x.append(np.random.choice(nts, p=p_emit_nt[s]))

        # Assign the rho value
        r.append(np.random.choice(rhos, p=p_emit_rho[s]))

    x = np.array(x)
    r = np.array(r)
    L = len(r)

    # Normalize rho by mean rho over the simulated sequence
    if np.mean(r) == 0: r = np.zeros(L)
    else: r = r/np.mean(r)

    return x, r


def simulate_seq_any(p_emit_nt,
                     p_emit_rho,
                     p_trans,
                     nts,
                     rho_ranges):
    """
    Sample a sequence from the model (can be canonical or noncanonical)
    Args:
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
    Returns:
        x: simulated nt sequence
        r: simulated rho sequence
        path: simulated path
    Raises:

    """
    # Sample a path using the transitions
    path = sample_transitions(p_trans)

    # Sample a nt and rho sequence using the emissions
    x, r = sample_emissions(path,
                            p_emit_nt,
                            p_emit_rho,
                            nts,
                            rho_ranges)

    return x, r, path


def simulate_seq_canonical(p_emit_nt,
                           p_emit_rho,
                           p_trans,
                           nts,
                           rho_ranges,
                           altorfs_idx):
    """
    Sample a canonical sequence from the model
    Args:
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
    Returns:
        x: simulated nt sequence
        r: simulated rho sequence
        path: simulated path
    Raises:

    """
    altorfs = altorfs_idx['uORF'] + \
              altorfs_idx['dORF'] + \
              altorfs_idx['PRF'] + \
              altorfs_idx['SCR']

    # If the simulated sequence contains non-canonical states,
    # regenerate until it doesn't
    noncanonical = True
    while noncanonical:

        # Sample a path using the transitions
        path = sample_transitions(p_trans)

        # Check if non-canonical
        noncanonical = any([path[i] in altorfs for i in range(len(path))])

    # Sample a nt and rho sequence using the emissions
    x, r = sample_emissions(path,
                            p_emit_nt,
                            p_emit_rho,
                            nts,
                            rho_ranges)

    return x, r, path


def check_unique(path,
                 altorfs_idx,
                 event):
    """
    Determine whether there is more than one simulated altORF in a path
    Args:
        path: simulated path
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        event: altORF event type to simulate
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
    Returns:
        unique: boolean indicating whether there is only one altORF
    Raises:

    """
    # Check that the altORF event is only one of its kind in the path
    unique = True
    if event=='uORF':
        if len([i for i in path if i in altorfs_idx['start_uORF']]) > 3:
            unique = False
    elif event=='dORF':
        if len([i for i in path if i in altorfs_idx['start_dORF']]) > 3:
            unique = False
    elif event=='uORFdORF':
        if len([i for i in path if i in (altorfs_idx['start_uORF'] + \
                                         altorfs_idx['start_dORF'])]) > 3:
            unique = False
    elif event=='pPRF':
        if len([i for i in path if i==altorfs_idx['X1']]) > 1:
            unique = False
    elif event=='mPRF':
        if len([i for i in path if i==altorfs_idx['X2']]) > 1:
            unique = False
    elif event=='SCR':
        if len([i for i in path if i in altorfs_idx['SCR']]) > 3:
            unique = False

    return unique


def simulate_seq_noncanonical(p_emit_nt,
                              p_emit_rho,
                              p_trans,
                              p_trans_rev,
                              nts,
                              rho_ranges,
                              altorfs_idx,
                              event):
    """
    Sample a noncanonical sequence from the model
    Args:
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        p_trans_rev: the reverse transition probability matrix
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        event: altORF event type to simulate
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
    Returns:
        x: simulated nt sequence
        r: simulated rho sequence
        path: simulated path
    Raises:

    """
    M = len(p_trans) # Number of states

    # Choose a state path using transition probabilities
    # Start at the noncanonical state
    if event=='uORF':
        idx = np.random.choice(altorfs_idx['start_uORF'])
    elif event=='dORF':
        idx = np.random.choice(altorfs_idx['start_dORF'])
    elif event=='uORFdORF':
        idx = np.random.choice(altorfs_idx['start_uORF'] + \
                               altorfs_idx['start_dORF'])
    elif event=='pPRF':
        idx = altorfs_idx['X1']
    elif event=='mPRF':
        idx = altorfs_idx['X2']
    elif event=='SCR':
        idx = np.random.choice(altorfs_idx['SCR'])
    else:
        logger.error(('Invalid event type specified for simulating sequence. '
                      'Choose from: uORF, dORF, uORFdORF, pPRF, mPRF, SCR'))
    path = [idx]

    # Add the next state and check
    next_state = np.random.choice(M, p=p_trans[idx+1]) - 1
    # For +PRF make sure the next state is not a -PRF
    if event=='pPRF':
        while next_state == altorfs_idx['X2']:
            next_state = np.random.choice(M, p=p_trans[idx+1]) - 1

    # Continue adding next states until you reach the end state
    while next_state != -1:
        curr_state = next_state
        path.append(curr_state)
        next_state = np.random.choice(M, p=p_trans[curr_state+1]) - 1

    # Continue adding previous states until you reach the begin state
    prev_state = np.random.choice(M, p=p_trans_rev[idx+1]) - 1
    while prev_state != -1:
        curr_state = prev_state
        path.insert(0, prev_state)
        prev_state = np.random.choice(M, p=p_trans_rev[curr_state+1]) - 1

    # Check that the altORF is only one of its kind in the simulated sequence
    # If not, resimulate until it is
    unique = check_unique(path, altorfs_idx, event)
    if not unique:
        return simulate_seq_noncanonical(p_emit_nt,
                                         p_emit_rho,
                                         p_trans,
                                         p_trans_rev,
                                         nts,
                                         rho_ranges,
                                         altorfs_idx,
                                         event)

    # Sample a nt and rho sequence using the emissions
    x, r = sample_emissions(path,
                            p_emit_nt,
                            p_emit_rho,
                            nts,
                            rho_ranges)

    return x, r, path


def scale_coverage(r, coverage):
    """
    Scale riboseq read coverage to reach target coverage (reads/nt)
    Args:
        r: simulated rho sequence
        coverage: mean transcript coverage to simulate (reads/nt)
    Returns:
        r: updated simulated rho sequence
    Raises:

    """
    # Scale the coverage as necessary
    if coverage > 0:

        # Convert riboseq values to a probability distribution
        r_prob = r/sum(r)

        # Calculate how many reads to assign to the transcript to generate
        # the target coverage
        reads = int(coverage * len(r))

        # Assign reads to the simulated sequence based on the probability
        # distribution and the target coverage
        L = len(r)
        r = np.zeros(L)
        for _ in range(reads):
            i = np.random.choice(L, p=r_prob)
            r[i] += 1

        # Renormalize the reads
        if np.mean(r) == 0: r = np.zeros(L)
        else: r = r/np.mean(r)

    else:
        logger.warning(('Cannot scale coverage for sequence. '
                        'Target coverage is less than or equal to zero.'))

    return r


def simulate_seq(p_emit_nt,
                 p_emit_rho,
                 p_trans,
                 p_trans_rev,
                 nts,
                 rho_ranges,
                 altorfs_idx,
                 event,
                 coverage=0):
    """
    Sample a sequence of the given event type from the model
    Args:
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        p_trans_rev: the reverse transition probability matrix
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        event: event type to simulate
            any = any ORF(s), canonical or noncanonical
            ORF = canonical ORF
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
        coverage: (optional) mean transcript coverage to simulate (reads/nt),
            default is not to scale the rho values at all
    Returns:
        x: simulated nt sequence
        r: simulated rho sequence
        path: simulated path
    Raises:

    """
    # Simulate a sequence
    if event=='any':
        x, r, path = simulate_seq_any(p_emit_nt,
                                      p_emit_rho,
                                      p_trans,
                                      nts,
                                      rho_ranges)
    elif event=='ORF':
        x, r, path = simulate_seq_canonical(p_emit_nt,
                                            p_emit_rho,
                                            p_trans,
                                            nts,
                                            rho_ranges,
                                            altorfs_idx)

    elif event in ['uORF','dORF','uORFdORF','pPRF','mPRF','SCR']:
        x, r, path = simulate_seq_noncanonical(p_emit_nt,
                                               p_emit_rho,
                                               p_trans,
                                               p_trans_rev,
                                               nts,
                                               rho_ranges,
                                               altorfs_idx,
                                               event)

    else:
        logger.error(('Invalid event type specified for simulating sequence. '
                      'Choose from: uORF, dORF, pPRF, mPRF, SCR, any, ORF'))

    # If there are no reads, run again recursively
    if (sum(r)==0):
        logger.warning('No reads simulated. Simulating new sequence...')
        return simulate_seq(p_emit_nt,
                            p_emit_rho,
                            p_trans,
                            p_trans_rev,
                            nts,
                            rho_ranges,
                            altorfs_idx,
                            event,
                            coverage)

    # Scale coverage to reach target coverage
    r = scale_coverage(r, coverage)

    return x, r, path


def simulation(p_emit_nt,
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
               event='any'):
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
        event: (optional) event type to look for
            ORF = canonical ORF
            uORF = uORF only
            dORF = dORF only
            uORFdORF = uORF or dORF
            pPRF = +1 frameshift
            mPRF = -1 frameshift
            SCR = stop codon readthrough
            default is any event
    Returns:
        true_path: simulated path
        pred_path: predicted path
    Raises:

    """
    # Reset the seed for each run (to avoid multiprocessing purgatory where
    # it reuses the exact same seed for multiple runs...)
    np.random.seed()

    # Simulate sequence
    x, r, true_path = simulate_seq(p_emit_nt,
                                   p_emit_rho,
                                   p_trans,
                                   p_trans_rev,
                                   nts,
                                   rho_ranges,
                                   altorfs_idx,
                                   event,
                                   coverage)

    # Generate the emission probabilities for this sequence
    log_emit = convert_emissions(x,
                                 r,
                                 log_emit_nt,
                                 log_emit_rho,
                                 nts,
                                 rho_ranges)

    # Predict the viterbi state path
    pred_path, P = viterbi(log_emit,
                           log_trans,
                           prev_states_options,
                           next_states_options)

    return true_path, pred_path
