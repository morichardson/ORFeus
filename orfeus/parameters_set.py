#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Set initial transition and emission probs
# Author: Mary Richardson
# Date: 2020.07.07, last update 2023.03.02
# -----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import itertools as it

from parameters_states import *
from parameters_transitions import *
from parameters_emissions import *


NUCLEOTIDES = ['A','C','G','T']
CODONS = [''.join(i) for i in it.product(NUCLEOTIDES, repeat = 3)]


def set_parameters(data_df,
                   f_5UTR,
                   f_3UTR,
                   len_5UTR,
                   len_3UTR,
                   len_orf,
                   len_uorf,
                   len_dorf,
                   start_codons,
                   sense_codons,
                   stop_codons,
                   num_bins,
                   logspace,
                   min_coverage,
                   max_coverage,
                   pool,
                   fit,
                   parameters):
    """
    Set the model transition and emission probabilities
    Args:
        data_df: pandas DataFrame with all transcript info
        f_5UTR: float indicating the fraction of ORFs having 5'UTRs
        f_3UTR: float indicating the fraction of ORFs having 3'UTRs
        len_5UTR: float indicating the mean length of 5'UTRs
        len_3UTR: float indicating the mean length of 3'UTRs
        len_orf: float indicating the mean length of ORFs
        len_uorf: float indicating the mean length of uORFs
        len_dorf: float indicating the mean length of dORFs
        start_codons: list of valid start codons
        sense_codons: list of valid sense codons
        stop_codons: list of valid stop codons
        num_bins: int number of bins to group the riboseq data into
        logspace: boolean indicating whether to bin in logspace
        min_coverage: int minimum mean reads per transcript
        max_coverage: int maximum mean reads per transcript
        pool: boolean indicating whether to pool all codons of a given type
            (start, stop, sense) when calculating riboseq emissions
        fit: boolean indicating whether to fit a log-normal distribution to
            the observed non-zero values when calculating riboseq emissions
        parameters: list of parameter values for
            alpha, beta, gamma, delta, zeta
    Returns:
        states: list of state names
        states_idx: dictionary mapping state type to
            state indices (1-indexed)
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        p_emit_nt: the nucleotide emission probability matrix
        p_emit_rho: the riboseq emission probability matrix
        p_trans: the transition probability matrix
        p_trans_rev: the reverse transition probability matrix
        prev_state_options: list of indices of all possible
            (nonzero probability) previous states
        next_state_options: list of indices of all possible
            (nonzero probability) next states
    """
    # Only train on annotated regions of transcripts
    data_df.drop(data_df[~data_df['feature'] \
           .isin(['5UTR', 'ORF', '3UTR'])].index, inplace=True)

    # Set constants
    f_5UTR, f_3UTR, len_5UTR, len_3UTR, len_orf, len_uorf, len_dorf, \
    start_codon_freqs, sense_codon_freqs, stop_codon_freqs = \
        set_constants(data_df,
                    CODONS,
                    f_5UTR,
                    f_3UTR,
                    len_5UTR,
                    len_3UTR,
                    len_orf,
                    len_uorf,
                    len_dorf,
                    start_codons,
                    sense_codons,
                    stop_codons)

    # Set model states
    states, states_idx, altorfs_idx = \
        set_states(start_codon_freqs,
                    sense_codon_freqs,
                    stop_codon_freqs,
                    parameters)

    # Set emission probabilities
    p_emit_nt = \
        set_nt_emission_probs(states,
                    states_idx,
                    NUCLEOTIDES)

    p_emit_rho, rho_ranges = \
        set_rho_emission_probs(data_df,
                    states,
                    states_idx,
                    num_bins,
                    logspace,
                    min_coverage,
                    max_coverage,
                    pool,
                    fit)

    # Set transition probabilities
    p_trans = \
        set_transition_probs(states,
                    states_idx,
                    f_5UTR,
                    f_3UTR,
                    len_3UTR,
                    len_5UTR,
                    len_orf,
                    len_uorf,
                    len_dorf,
                    start_codon_freqs,
                    sense_codon_freqs,
                    stop_codon_freqs,
                    parameters)

    p_trans_rev = rev_transition_probs(p_trans)

    prev_states_options = get_prev_states(p_trans)
    next_states_options = get_next_states(p_trans)

    all_parameters = np.array([states, states_idx, altorfs_idx, \
           p_emit_nt, p_emit_rho, rho_ranges, p_trans, p_trans_rev, \
           prev_states_options, next_states_options], dtype=object)

    return all_parameters
