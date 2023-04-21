#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Set model states
# Author: Mary Richardson
# Date: 2020.07.07, last update 2023.03.02
# -----------------------------------------------------------------------------

import numpy as np


def set_states(start_codon_freqs,
               sense_codon_freqs,
               stop_codon_freqs,
               parameters):
    """
    Generate a list of all states for this model
    Args:
        start_codon_freqs: dict mapping each start codon to its frequency
        stop_codons_freqs: dict mapping each stop codon to its frequency
        sense_codons_freqs: dict mapping each sense codon to its frequency
        parameters: list of parameter values for
            alpha, beta, gamma, delta, zeta
    Returns:
        states: list of state names
        states_idx: dictionary mapping state type to
            state indices (1-indexed)
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
    Raises:
        ValueError: Invalid parameter values. Make sure all parameter values
        are greater than or equal to zero.
    """
    # Unpack parameters
    [alpha, beta, gamma, delta, zeta] = parameters

    # List states and their corresponding indices in the
    # transition and emission matrices
    states = []
    states_idx = {}

    # Create dictionary to store all indices for each type of state
    altorfs_idx = {'5UTR':[],
                   '3UTR':[],
                   'ORF':[],
                   'start':[],
                   'sense':[],
                   'stop':[],
                   'X1':[],
                   'X2':[],
                   'PRF':[],
                   'SCR':[],
                   'uORF':[],
                   'dORF':[]}

    # Begin counting indices
    idx = 0

    # Add UTR states
    for state in ['5UTR','3UTR']:
        states += [state]
        states_idx[state] = idx
        altorfs_idx[state] = idx
        idx += 1

    # Add start codon states
    start_states = [(codon + pos) for codon in sorted(start_codon_freqs) \
                                  for pos in ['1','2','3']]
    states += [(state + '#') for state in start_states]
    n = len(start_states)
    states_idx['start'] = list(range(idx, idx+n))
    altorfs_idx['start'] = list(range(idx, idx+n))
    altorfs_idx['ORF'] = list(range(idx, idx+n))
    idx += n

    # Add sense codon states
    sense_states = [(codon + pos) for codon in sorted(sense_codon_freqs) \
                                  for pos in ['1','2','3']]
    states += [(state) for state in sense_states]
    n = len(sense_states)
    states_idx['sense'] = list(range(idx, idx+n))
    altorfs_idx['sense'] = list(range(idx, idx+n))
    altorfs_idx['ORF'] += list(range(idx, idx+n))
    idx += n

    # Add stop codon states
    stop_states = [(codon + pos) for codon in sorted(stop_codon_freqs) \
                                 for pos in ['1','2','3']]
    states += [(state + '*') for state in stop_states]
    n = len(stop_states)
    states_idx['stop'] = list(range(idx, idx+n))
    altorfs_idx['stop'] = list(range(idx, idx+n))
    altorfs_idx['ORF'] += list(range(idx, idx+n))
    idx += n

    # Add PRF states
    if alpha > 0:

        # Add X1 state
        states += ['X1']
        states_idx['X1'] = idx
        altorfs_idx['X1'] = idx
        altorfs_idx['PRF'] = [idx]
        idx += 1

    if alpha > 0 and beta > 0:

        # Add X2 state
        states += ['X2']
        states_idx['X2'] = idx
        altorfs_idx['X2'] = idx
        altorfs_idx['PRF'] += [idx]
        idx += 1

    # Add SCR states
    if gamma > 0:

        # Add stop codon readthrough states
        stop_states = [(codon + pos) for codon in sorted(stop_codon_freqs) \
                                     for pos in ['1','2','3']]
        states += [(state + '**') for state in stop_states]
        n = len(stop_states)
        states_idx['stop_SCR'] = list(range(idx, idx+n))
        altorfs_idx['SCR'] = list(range(idx, idx+n))
        idx += n

    # Add uORF and dORF states
    if delta > 0:

        # Add uORF start codon states
        start_states = [(codon + pos) for codon in sorted(start_codon_freqs) \
                                      for pos in ['1','2','3']]
        states += [(state + '#u') for state in start_states]
        n = len(start_states)
        states_idx['start_uORF'] = list(range(idx, idx+n))
        altorfs_idx['start_uORF'] = list(range(idx, idx+n))
        altorfs_idx['uORF'] = list(range(idx, idx+n))
        idx += n

        # Add uORF sense codon states
        sense_states = [(codon + pos) for codon in sorted(sense_codon_freqs) \
                                      for pos in ['1','2','3']]
        states += [(state + 'u') for state in sense_states]
        n = len(sense_states)
        states_idx['sense_uORF'] = list(range(idx, idx+n))
        altorfs_idx['uORF'] += list(range(idx, idx+n))
        idx += n

        # Add uORF stop codon states
        stop_states = [(codon + pos) for codon in sorted(stop_codon_freqs) \
                                     for pos in ['1','2','3']]
        states += [(state + '*u') for state in stop_states]
        n = len(stop_states)
        states_idx['stop_uORF'] = list(range(idx, idx+n))
        altorfs_idx['stop_uORF'] = list(range(idx, idx+n))
        altorfs_idx['uORF'] += list(range(idx, idx+n))
        idx += n

        # Add dORF start codon states
        start_states = [(codon + pos) for codon in sorted(start_codon_freqs) \
                                      for pos in ['1','2','3']]
        states += [(state + '#d') for state in start_states]
        n = len(start_states)
        states_idx['start_dORF'] = list(range(idx, idx+n))
        altorfs_idx['start_dORF'] = list(range(idx, idx+n))
        altorfs_idx['dORF'] = list(range(idx, idx+n))
        idx += n

        # Add dORF sense codon states
        sense_states = [(codon + pos) for codon in sorted(sense_codon_freqs) \
                                      for pos in ['1','2','3']]
        states += [(state + 'd') for state in sense_states]
        n = len(sense_states)
        states_idx['sense_dORF'] = list(range(idx, idx+n))
        altorfs_idx['dORF'] += list(range(idx, idx+n))
        idx += n

        # Add dORF stop codon states
        stop_states = [(codon + pos) for codon in sorted(stop_codon_freqs) \
                                     for pos in ['1','2','3']]
        states += [(state + '*d') for state in stop_states]
        n = len(stop_states)
        states_idx['stop_dORF'] = list(range(idx, idx+n))
        altorfs_idx['stop_dORF'] = list(range(idx, idx+n))
        altorfs_idx['dORF'] += list(range(idx, idx+n))
        idx += n

    return states, states_idx, altorfs_idx
