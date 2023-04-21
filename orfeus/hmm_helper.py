#!/usr/bin/env python

# -----------------------------------------------------------------------------
# HMM helper functions
# Author: Mary Richardson
# Date: 2020.01.08
# -----------------------------------------------------------------------------

import logging
import numpy as np
from scipy.special import logsumexp

logger = logging.getLogger(__name__)

# Ignore divide by zero error for log
np.seterr(divide='ignore')


def logdotexp(A, B):
    """
    Code inspiration: https://stackoverflow.com/questions/23630277/numerically-stable-way-to-multiply-log-probability-matrices-in-numpy
    Given matrices A (one-hot encoded) and B (log space), calculate the
    dot product computed in logspace using the logsumexp function such that
    the log values of A are multiplied by B and then summed (in linear space)
    Args:
        A: numpy array (one-hot encoded)
        B: numpy array (log space)
    Returns:
        C: numpy array resulting from the dot product of A and B

    """
    max_B = np.max(A)
    C = np.log(np.dot(A, np.exp(B - max_B)))
    C += max_B

    return C


def one_hot_nt(x, nts):
    """
    Generate the one-hot encoding for this nucleotide sequence
    Args:
        x: the nucleotide sequence
        nts: list of valid nucleotides
    Returns:
        one_hot_xt: the one-hot-encoded matrix for this nucleotide sequence
    Raises:
        ValueError: Invalid nucleotide sequence
    """
    # Length of sequence
    L = len(x)

    # Number of nucleotides
    N = len(nts)

    # Convert the sequence to numbers (A=0, C=1, G=2, T=3)
    one_hot_x = np.zeros((L,N))
    for i,nt in enumerate(x):
        if nt in nts:
            one_hot_x[i, nts.index(nt)] = 1
        elif nt=='N': # Any nucleotide
            one_hot_x[i, :] = np.ones(N)
        elif nt=='R': # Purine
            one_hot_x[i, nts.index('A')] = 1
            one_hot_x[i, nts.index('G')] = 1
        elif nt=='Y': # Pyrimidine
            one_hot_x[i, nts.index('C')] = 1
            one_hot_x[i, nts.index('T')] = 1
        else:
            logging.error(('Invalid nucleotide sequence. '
                           'Unexpected character in sequence: %s' % i))

    return one_hot_x


def one_hot_riboseq(r, rho_ranges):
    """
    Generate the one-hot encoding for this rho sequence
    Args:
        r: the rho sequence
        rho_ranges: array of riboseq emission bin ranges
    Returns:
        one_hot_rho: the one-hot-encoded matrix for this rho sequence
    Raises:
        ValueError: Invalid ribosome profiling values
    """
    # Length of sequence
    L = len(r)

    # Number of discretized riboseq ranges
    R = len(rho_ranges)-1

    try:
        # Convert the riboseq sequence to discretized ranges
        # Include the upper limit of the range and shift to be zero-indexed
        y = np.digitize(r, rho_ranges, right=True) - 1

        # For anything beyond the max observed value,
        # assign it to the uppermost bin
        y[y >= R] = R-1

        # Create a one-hot encoding for the sequence
        one_hot_r = np.eye(R)[np.array(y)]

    except:
        logging.error('Invalid ribosome profiling values.')

    return one_hot_r


def convert_emissions(x, r, log_emit_nt, log_emit_rho, nts, rho_ranges):
    """
    Generate the log emission probabilities at each position of a sequence
    Note: Returns the emission probabilities at each sequence position as
    an L x M matrix. This allows matrix addition down the line, since the
    forward and backward matrices are L x M
    Args:
        x: the nucleotide sequence
        r: the rho sequence
        log_emit_nt: log nucleotide emission probabilities
        log_emit_rho: log riboseq emission probabilities
        nts: list of valid nucleotides
        rho_ranges: array of riboseq emission bin ranges
    Returns:
        one_hot_x: the one-hot-encoded matrix for this nucleotide sequence
        one_hot_r: the one-hot-encoded matrix for this rho sequence
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
    """
    # Calculate the nucleotide emissions for this sequence
    one_hot_x = one_hot_nt(x, nts)
    log_emit_x = logdotexp(one_hot_x, log_emit_nt.T)

    # Create a one-hot encoding for this riboseq sequence
    one_hot_r = one_hot_riboseq(r, rho_ranges)
    log_emit_r = logdotexp(one_hot_r, log_emit_rho.T)

    # Calculate the total emissions for this sequence
    log_emit = log_emit_x + log_emit_r

    return log_emit
