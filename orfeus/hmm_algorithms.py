#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Hidden Markov Model Algorithms
# Author: Mary Richardson
# Date: 2020.01.08
# -----------------------------------------------------------------------------

import logging
import numpy as np
from scipy.special import logsumexp

logger = logging.getLogger(__name__)

# Ignore divide by zero error for log
np.seterr(divide='ignore')


def viterbi(log_emit,
            log_trans,
            prev_states_options,
            next_states_options):
    """
    Determine the best path for a sequence using the Viterbi algorithm
    A note about dimensions:
        transition matrix includes begin state (1-indexed)
        next and prev state options include begin state (1-indexed)
        emission matrix excludes begin state (0-indexed)
        viterbi matrix excludes begin state (0-indexed)
    Args:
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
        log_trans: log transition probability matrix
        prev_state_options: list of indices of all possible
            (nonzero probability) previous states
        next_state_options: list of indices of all possible
            (nonzero probability) next states
    Returns:
        path: Viterbi path
        P: total probability of Viterbi path
    Raises:

    """
    L = len(log_emit)                   # Length of observed sequence
    M = len(log_trans) - 1              # Number of states
    V = np.log(np.zeros([L,M]))         # Viterbi score matrix
    ptr = np.zeros([L,M])               # Viterbi traceback matrix
    path = np.zeros(L, dtype=np.int32)  # Most probable viterbi path

    # Initialization

    # Only compute transitions to nonzero probability next states
    first_states = next_states_options[0]

    # All sequences start in the begin state
    for k in first_states:
        V[0,k-1] = log_emit[0,k-1] + log_trans[0,k]

    # Recursion

    # For each position in the observed sequence 1..L
    for i in range(1, L):

        # Only compute transitions from nonzero probability current states
        curr_states = np.argwhere(np.isfinite(V[i-1,:])).flatten() + 1
        next_states = list(set([k for j in curr_states \
                                  for k in next_states_options[j] if k != 0]))

        # For each possible next state
        for k in next_states:
            Z = [V[i-1,j-1] + log_trans[j,k] for j in curr_states]

            # Update the viterbi score matrix
            V[i,k-1] = log_emit[i,k-1] + np.max(Z)

            # Update the viterbi traceback matrix
            ptr[i,k-1] = curr_states[np.argmax(Z)] - 1

    # Termination

    # Only compute transitions from nonzero probability prev states
    last_states = prev_states_options[0]

    # All sequences end in the end state
    Z = [V[L-1,j-1] + log_trans[j,0] for j in last_states]

    # Total probability of sequence
    P = np.max(Z)

    # Traceback

    # End in the end state
    path[L-1] = last_states[np.argmax(Z)] - 1

    # For each position in the observed sequence L-1..0
    for i in reversed(range(L-1)):

        # Traceback to the most probable previous state
        path[i] = ptr[i+1, path[i+1]]

    return path, P


def forward(log_emit,
            log_trans,
            prev_states_options,
            next_states_options):
    """
    Calculate the probability of a sequence using the Forward algorithm
    Args:
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
        log_trans: log transition probability matrix
        prev_state_options: list of indices of all possible
            (nonzero probability) previous states
        next_state_options: list of indices of all possible
            (nonzero probability) next states
    Returns:
        F: Forward probabiliy matrix
        P: total probability
    Raises:

    """
    L = len(log_emit)            # Length of observed sequence
    M = len(log_trans) - 1       # Number of states
    F = np.log(np.zeros([L,M]))  # Forward score matrix

    # Initialization

    # Only compute transitions to nonzero probability next states
    first_states = next_states_options[0]

    # All sequences start in the begin state
    for k in first_states:
        F[0,k-1] = log_emit[0,k-1] + log_trans[0,k]

    # Recursion

    # For each position in the observed sequence 1..L
    for i in range(1, L):

        # Only compute transitions from nonzero probability current states
        curr_states = np.argwhere(np.isfinite(F[i-1,:])).flatten() + 1
        next_states = list(set([k for j in curr_states \
                                  for k in next_states_options[j] if k != 0]))

        # For each possible next state
        for k in next_states:
            Z = [F[i-1,j-1] + log_trans[j,k] for j in curr_states]

            # Update the forward score matrix
            F[i,k-1] = log_emit[i,k-1] + logsumexp(Z)

    # Termination

    # Only compute transitions from nonzero probability prev states
    last_states = prev_states_options[0]

    # All sequences end in the end state
    Z = [F[L-1,j-1] + log_trans[j,0] for j in last_states]

    # Total probability of sequence
    P = logsumexp(Z)

    return F, P


def backward(log_emit,
             log_trans,
             prev_states_options,
             next_states_options):
    """
    Calculate the probability of a sequence using the Backward algorithm
    Args:
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
        log_trans: log transition probability matrix
        prev_state_options: list of indices of all possible
            (nonzero probability) previous states
        next_state_options: list of indices of all possible
            (nonzero probability) next states
    Returns:
        B: Backward probabiliy matrix
        P: total probability
    Raises:

    """
    L = len(log_emit)            # Length of observed sequence
    M = len(log_trans) - 1       # Number of states
    B = np.log(np.zeros([L,M]))  # Backward score matrix

    # Initialization

    # Only compute transitions from nonzero probability prev states
    last_states = prev_states_options[0]

    # Only compute transitions from nonzero probability prev states
    for k in last_states:

        # All sequences end in the end state
        B[L-1,k-1] = log_trans[k,0]

    # Recursion

    # For each position in the observed sequence L-1 to 1
    for i in reversed(range(L-1)):

        # Only compute transitions from nonzero probability current states
        curr_states = np.argwhere(np.isfinite(B[i+1,:])).flatten() + 1
        prev_states = list(set([j for k in curr_states \
                                  for j in prev_states_options[k] if j != 0]))

        # For each possible previous state
        for j in prev_states:
            Z = [log_emit[i+1,k-1] + B[i+1,k-1] + log_trans[j,k] \
                 for k in curr_states]

            # Update the backward score matrix
            B[i,j-1] =  logsumexp(Z)

    # Termination

    # Only compute transitions from nonzero probability prev states
    first_states = next_states_options[0]

    # All sequences start in the begin state
    Z = [log_emit[0,j-1] + log_trans[0,j] + B[0,j-1] for j in first_states]

    # Total probability of sequence
    P = logsumexp(Z)

    return B, P


def posterior_decoding(log_emit,
                       log_trans,
                       prev_states_options,
                       next_states_options):
    """
    Determine the best path for a sequence using the Posterior Decoding
    algorithm
    Args:
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
        log_trans: log transition probability matrix
        prev_state_options: list of indices of all possible
            (nonzero probability) previous states
        next_state_options: list of indices of all possible
            (nonzero probability) next states
    Returns:
        path: most probable Posterior Decoding path
        P: total probability of chosen Posterior Decoding path
    Raises:

    """
    # Compute the forward and backward scores for the sequence
    # (total probability is the same when computed with forward or backward)
    F, P_F = forward(log_emit,
                     log_trans,
                     prev_states_options,
                     next_states_options)
    B, P_B = backward(log_emit,
                     log_trans,
                     prev_states_options,
                     next_states_options)

    # Calculate the posterior decoding result at each position
    log_P = (F + B) - P_F
    P = np.exp(log_P)

    # Calculate the path using posterior decoding
    path = np.argmax(P, axis=1)

    return path, P


def score_path(path,
               log_emit,
               log_trans):
    """
    Calculates the score of the given state path, using the given model
    parameters
    Args:
        path: state path
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
        log_trans: log transition probability matrix
    Returns:
        P: total probability of the path
    Raises:

    """
    # Length of observed sequence
    L = len(path)

    if L==0:
        return np.nan

    # Add the probability of the initial transition
    # from the begin state (in log space)
    P = log_trans[0, path[0]+1]

    # For each position in the sequence 0..L-1
    # add the probability of the transition to this residue and
    # the emission for this residue (in log space)
    for i in range(0, L-1):
        P += log_emit[i, path[i]] + log_trans[path[i]+1, path[i+1]+1]

    # Add the probability of the final emission and
    # the transition to the end state (in log space)
    P += log_emit[L-1, path[L-1]] + log_trans[path[L-1]+1, 0]

    return P
