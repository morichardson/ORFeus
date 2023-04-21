#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Set parameters helper functions
# Author: Mary Richardson
# Date: 2020.07.07, last update 2023.03.02
# -----------------------------------------------------------------------------



def get_codon(i, states):
    """
    Get the codon sequence for this state
    Args:
        i: index of state
        states: list of all states in the model
    Returns:
        codon: codon sequence
    """
    state = states[i] # Label for this state (e.g. 'ATG1#')
    codon = state[:3] # Sequence of this codon
    return codon


def get_nt(i, states):
    """
    Get the subcodon nucleotide position for this state
    Args:
        i: index of state
        states: list of all states in the model
    Returns:
        nt: nucleotide position in codon (1, 2, or 3)
    """
    state = states[i] # Label for this state (e.g. 'ATG1#')
    nt = int(state[3])  # Subcodon nucleotide of this codon
    nt = nt-1
    return nt


def normalize(d):
    """
    Normalize the values of a dictionary
    Args:
        d: dictionary with numerical values
    Returns:
        normalized_d: dictionary with values normalized so they sum to 1
    """
    factor = 1.0/sum(d.values())
    normalized_d = {k: v*factor for k, v in d.items()}
    return normalized_d
