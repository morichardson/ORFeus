#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Set emission probabilities
# Author: Mary Richardson
# Date: 2020.07.07, last update 2023.03.02
# -----------------------------------------------------------------------------

import logging
import numpy as np
import pandas as pd
import scipy.stats as ss

logger = logging.getLogger(__name__)

from parameters_helper import get_codon, get_nt


def set_nt_emission_probs(states,
                          states_idx,
                          nts):
    """
    Calculate the probability of observing each nucleotide in each state
    Args:
        states: list of all states in the model
        states_idx: dictionary mapping each general state type to its
            state index (1-indexed)
        nts: list of nucleotides
    Returns:
        p_emit_nt: the nucleotide emission probability matrix,
            i.e. if the state is ATG, the prob of A at the first nucleotide is 1
            and the prob of all other nucleotides is 0 (note that this array
            has one row per nucleotide, in the same order as the nts array)
    Raises:
        ValueError: Invalid nt emission probability matrix
    """
    # Prepare the nucleotide emission probabilities matrix
    M = len(states)
    N = len(nts)
    p_emit_nt = np.zeros((M,N))

    # Get a list of noncodon states (can be any nt)
    noncodon_states = [states_idx['5UTR'], states_idx['3UTR']]
    if 'X1' in states:
        noncodon_states.append(states_idx['X1'])
    if 'X2' in states:
        noncodon_states.append(states_idx['X2'])

    # Iterate through row for each possible observed nt
    for i in range(N):

        # For each possible state
        for s in range(M):

            if s in noncodon_states:
                # Any nt could be observed, so set value to 1/4
                p_emit_nt[s, i] = 1/4

            else:
                # Only the nt at this position in this codon could be observed
                codon = get_codon(s, states) # Sequence of this codon
                nt_pos = get_nt(s, states)   # Subcodon nucleotide

                # If this nt is in this state, change the value to 1
                if nts[i] == codon[nt_pos]:
                    p_emit_nt[s, i] = 1

    # Check that emission probabilities for each state sum to 1
    if not np.allclose(np.sum(p_emit_nt, axis=1), np.ones(len(p_emit_nt))):
        logging.error(('Invalid nucleotide emission probability matrix. '
                    'Nucleotide emissions for each state do not sum to 1.'))

    return p_emit_nt


def rho_bins(df, num_bins, logspace=False):
    """
    Fit a distribution to the observed riboseq values for each state
    Args:
        df: pandas DataFrame with all riboseq info
        num_bins: int number of bins to group the riboseq data into
        logspace: boolean indicating whether to bin in logspace
    Returns:
        rho_ranges: array of riboseq emission bin ranges

    """
    # List of ranges for discretizing riboseq observations
    if logspace:
        min_reads = np.nanmin(df['log_reads'])
        max_reads = np.nanmax(df['log_reads'])
    else:
        min_reads = np.nanmin(df['norm_reads'])
        max_reads = np.nanmax(df['norm_reads'])

    # Create R-1 discrete bins for non-zero values
    # (evenly spaced in log or linear space, depending on limits)
    rho_ranges = np.linspace(min_reads, max_reads, num=num_bins)

    # Make sure bin boundaries are specified in linear space
    if logspace:
        rho_ranges = np.exp(rho_ranges)

    # Add a bin for zeros (dummy lower boundary)
    rho_ranges[0] = 0
    rho_ranges = np.insert(rho_ranges, 0, -.0001) # Add a bin for zeros (dummy lower boundary)

    return rho_ranges


def fit_rho_emissions(df,
                      states,
                      states_idx,
                      rho_ranges,
                      pool):
    """
    Fit a distribution to the frequency of observed riboseq values for
        each state
    Args:
        df: pandas DataFrame with all riboseq info
        states: list of state names
        states_idx: dictionary mapping state type to
            state indices (1-indexed)
        rho_ranges: array of riboseq emission bin ranges
        pool: boolean indicating whether to pool all codons of a given type
            (start, stop, sense) when calculating riboseq emissions
    Returns:
        p_emit_rho: the riboseq emission probability matrix
    Raises:

    """
    # Isolate state labels
    df['nt_pos'].fillna(-1, inplace=True)
    df_states = df[['state','codon_type','nt_pos']].drop_duplicates()

    # Update SENSE states (aggregate of all SENSE codon states if pool)
    if pool:
        df_grouped = df.groupby(['codon_type','nt_pos']) \
                       .agg({'log_reads': lambda x: x.dropna().tolist(),
                             'norm_reads': lambda x: x.dropna().tolist()})
        df_grouped = df_states.join(df_grouped,
                                    on=['codon_type','nt_pos'],
                                    how='left')

    # Update START and STOP states
    # (aggregate all START or STOP codon states)
    else:
        df_grouped1 = df.groupby(['state']) \
                        .agg({'log_reads': lambda x: x.dropna().tolist(),
                              'norm_reads': lambda x: x.dropna().tolist()})
        df_grouped2 = df.groupby(['codon_type','nt_pos']) \
                        .agg({'log_reads': lambda x: x.dropna().tolist(),
                              'norm_reads': lambda x: x.dropna().tolist()})
        df_grouped1 = df_states.join(df_grouped1,
                                     on='state',
                                     how='left')
        df_grouped2 = df_states.join(df_grouped2,
                                     on=['codon_type','nt_pos'],
                                     how='left')

        df_grouped1 = df_grouped1[~df_grouped1['codon_type'] \
                                 .isin(['START','STOP'])]
        df_grouped2 = df_grouped2[df_grouped2['codon_type'] \
                                 .isin(['START','STOP'])]
        df_grouped = pd.concat([df_grouped1, df_grouped2])

    # Fit a log normal distribution to the nonzero reads
    df_grouped['p'] = df_grouped.apply(lambda row: \
                                 ss.norm.fit(row['log_reads']),
                                 axis=1)

    # Calculate the the fraction of zero reads
    df_grouped['q'] = df_grouped.apply(lambda row: \
                                 float((len(row['norm_reads']) -
                                 np.count_nonzero(row['norm_reads'])) /
                                 len(row['norm_reads'])),
                                 axis=1)

    # Make sure all states are in order
    df_grouped.set_index('state', inplace=True)
    df_grouped = df_grouped.reindex(states)

    # Keep only these parameters for each state
    p = df_grouped['p'].tolist()
    q = df_grouped['q'].tolist()

    # Set representative emission values (mean each bin greater than zero)
    emit = np.mean([rho_ranges[:-1], rho_ranges[1:]], axis=0)

    # Discretize the emissions
    M = len(states)
    R = len(rho_ranges)-1
    p_emit_rho = np.zeros([M,R])

    # Skip altORF states for now
    for i in (states_idx['5UTR'] + states_idx['ORF'] + states_idx['3UTR']):

        # Add the zero emissions
        if np.any(np.isnan(q[i])):
            p_emit_rho[i,0] = None
        else:
            p_emit_rho[i,0] = q[i]

        # Add the nonzero emissions
        for j in range(len(emit)-1):
            if np.any(np.isnan(p[i])):
                p_emit_rho[i,j+1] = None
            else:
                p_emit_rho[i,j+1] = ss.norm(*p[i]).pdf(emit[j]) * (1-q[i])

    return p_emit_rho


def nofit_rho_emissions(df,
                        states,
                        states_idx,
                        rho_ranges,
                        pool):
    """
    Calculate the frequency observed riboseq values for each state
    Args:
        df: pandas DataFrame with all riboseq info
        states: list of state names
        states_idx: dictionary mapping state type to
            state indices (1-indexed)
        rho_ranges: array of riboseq emission bin ranges
        pool: boolean indicating whether to pool all codons of a given type
            (start, stop, sense) when calculating riboseq emissions
    Returns:
        p_emit_rho: the riboseq emission probability matrix
    Raises:

    """# Group the log riboseq reads by range
    # (include the upper limit of the range)
    df['rho_ranges'] = pd.cut(df['norm_reads'],
                                    rho_ranges,
                                    right=True)

    # Get an index of all states and log read ranges
    df_states = list(set(df['state']) & set(states))
    idx = pd.MultiIndex \
            .from_product([df_states, df['rho_ranges'].cat.categories],
                          names=['state', 'rho_ranges'])

    # For each state, count the reads in each range
    df = df.groupby(['state', 'rho_ranges']) \
           .size()

    # Reset the row names to make sure ALL log reads ranges are included
    # (fill in zeros for any missing ranges)
    df = df.reindex(index=idx, fill_value=0) \
           .reset_index()
    df.columns = ['state', 'rho_ranges', 'counts'] # Rename columns

    # Pivot to a MxR matrix of counts for each state
    df = df.pivot(index='state',
                  columns='rho_ranges',
                  values='counts')

    # Reorder the states
    df = df.reindex(states)

    # Emission prob is the frequency of counts in each state
    p_emit_rho = df.to_numpy()

    # Update SENSE states (mean of all SENSE codon states if pool)
    if pool:
        for nt in [0,1,2]:
            mean_sense = \
                np.nanmean(p_emit_rho[nt+states_idx['sense'][0]:states_idx['sense'][-1]+1:3],
                        axis=0)
            p_emit_rho[nt+states_idx['sense'][0]:states_idx['sense'][-1]+1:3] \
                = mean_sense

    return p_emit_rho


def set_rho_emission_probs(data_df,
                           states,
                           states_idx,
                           num_bins,
                           logspace,
                           min_coverage,
                           max_coverage,
                           pool,
                           fit):
    """
    Calculate the probability of observing observing each range of riboseq
        values in each state
    Args:
        data_df: pandas DataFrame with all transcript info
        states: list of all states in the model
        states_idx: dictionary mapping each general state type to its
            state index (1-indexed)
        num_bins: int number of bins to group the riboseq data into
        logspace: (optional) boolean indicating whether to bin in logspace
        min_coverage: (optional) int minimum mean reads per transcript
        max_coverage: (optional) int maximum mean reads per transcript
        pool: (optional) boolean indicating whether to pool all codons of a given type
            (start, stop, sense) when calculating riboseq emissions
        fit: (optional) boolean indicating whether to fit a log-normal distribution to
            the observed non-zero values when calculating riboseq emissions
    Returns:
        p_emit_rho: the riboseq emission probability matrix
        rho_ranges: array of riboseq emission bin ranges
    Raises:
        ValueError: Invalid riboseq emission probability matrix
    """
    # Filter down to states that are in the list of states
    df = data_df[data_df['state'].isin(states)]

    # Filter down to transcripts above min coverage and below max coverage
    df = df[(df['mean_reads']>=min_coverage) & (df['mean_reads']<=max_coverage)]

    # Generate riboseq emission bin ranges
    rho_ranges = rho_bins(df, num_bins, logspace)

    # Calculate the emission probabilities from the observed riboseq data
    # for all canonical states
    if fit:
        p_emit_rho = fit_rho_emissions(df,
                                       states,
                                       states_idx,
                                       rho_ranges,
                                       pool)
    else:
        p_emit_rho = nofit_rho_emissions(df,
                                       states,
                                       states_idx,
                                       rho_ranges,
                                       pool)

    # Update START states (mean of all START codon states)
    for nt in [0,1,2]:
        mean_start = np.nanmean(p_emit_rho[nt+states_idx['start'][0]:states_idx['start'][-1]+1:3], \
                              axis=0)
        p_emit_rho[nt+states_idx['start'][0]:states_idx['start'][-1]+1:3] \
                = mean_start
        if 'start_uORF' in states_idx:
            p_emit_rho[nt+states_idx['start_uORF'][0]:states_idx['start_uORF'][-1]+1:3] \
                = mean_start
        if 'start_dORF' in states_idx:
            p_emit_rho[nt+states_idx['start_dORF'][0]:states_idx['start_dORF'][-1]+1:3] \
                = mean_start

    # Update STOP states (mean of all STOP codon states)
    for nt in [0,1,2]:
        mean_stop = np.nanmean(p_emit_rho[nt+states_idx['stop'][0]:states_idx['stop'][-1]+1:3], \
                             axis=0)
        p_emit_rho[nt+states_idx['stop'][0]:states_idx['stop'][-1]+1:3] \
                = mean_stop
        if 'stop_uORF' in states_idx:
            p_emit_rho[nt+states_idx['stop_uORF'][0]:states_idx['stop_uORF'][-1]+1:3] \
                = mean_stop
        if 'stop_dORF' in states_idx:
            p_emit_rho[nt+states_idx['stop_dORF'][0]:states_idx['stop_dORF'][-1]+1:3] \
                = mean_stop
        if 'stop_SCR' in states_idx:
            p_emit_rho[nt+states_idx['stop_SCR'][0]:states_idx['stop_SCR'][-1]+1:3] \
                = mean_stop

    # Add uORF and dORF states
    if 'sense_uORF' in states_idx:
        p_emit_rho[states_idx['sense_uORF'][0]:states_idx['sense_uORF'][-1]+1] = \
            p_emit_rho[states_idx['sense'][0]:states_idx['sense'][-1]+1]
    if 'sense_dORF' in states_idx:
        p_emit_rho[states_idx['sense_dORF'][0]:states_idx['sense_dORF'][-1]+1] = \
        p_emit_rho[states_idx['sense'][0]:states_idx['sense'][-1]+1]

    # Add PRF states
    mean_orf = np.nanmean(p_emit_rho[states_idx['sense'][0]:states_idx['sense'][-1]+1],
                          axis=0)
    if 'X1' in states:
        p_emit_rho[states_idx['X1']] = mean_orf # X1
        if 'X2' in states:
            p_emit_rho[states_idx['X2']] = mean_orf # X2

    # Replace missing states
    for row in np.where(np.isnan(p_emit_rho).any(axis=1))[0]:
        logging.warning('Emissions missing for ' + states[row])
        p_emit_rho[row] = mean_orf

    # Add a pseudocount to all counts
    p_emit_rho += 1

    # Normalize to make sure the probs for each state sum to 1
    p_emit_rho = p_emit_rho / p_emit_rho.sum(axis=1, keepdims=True)

    # Check that emission probabilities for each state sum to 1
    if not np.allclose(np.sum(p_emit_rho, axis=1), np.ones(len(p_emit_rho))):
        logging.error(('Invalid riboseq emission probability matrix. '
                        'Riboseq emissions for each state do not sum to 1.'))

    return p_emit_rho, rho_ranges
