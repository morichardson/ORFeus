#!/usr/bin/env python

# -----------------------------------------------------------------------------
# HMM output functions
# Author: Mary Richardson
# Date: 2021.10.25
# -----------------------------------------------------------------------------

import logging
import numpy as np
import itertools as it
import multiprocessing as mp

logger = logging.getLogger(__name__)


AMINO_ACIDS = {'TTT':'F','TTC':'F',
               'TTA':'L','TTG':'L',
               'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
               'TAT':'Y','TAC':'Y',
               'TAA':'*','TAG':'*','TGA':'*',
               'TGT':'C','TGC':'C',
               'TGG':'W',
               'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
               'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
               'CAT':'H','CAC':'H',
               'CAA':'Q','CAG':'Q',
               'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
               'ATT':'I','ATC':'I','ATA':'I',
               'ATG':'M',
               'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
               'AAT':'N','AAC':'N',
               'AAA':'K','AAG':'K',
               'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
               'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
               'GAT':'D','GAC':'D',
               'GAA':'E','GAG':'E',
               'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}


def get_ranges(values):
    """
    Code source: https://stackoverflow.com/questions/4628333/converting-a-list-of-integers-into-range-in-python
    Yields ranges of values present in a list
    (e.g. [2,3,4,8,9] would give ranges 2-4 and 8-9)
    Args:
        values: list of integer values
    Returns:

    Raises:

    """
    for a, b in it.groupby(enumerate(values), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def get_events(path, altorfs_idx):
    """
    Finds all altORF events and their indices in the given path
    Args:
        path: the state path
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
    Returns:
        altorfs: list of all ORF and altORF events and corresponding indices
            in the path, in string format for printing
    Raises:

    """
    altorfs = []

    if len(path)==0:
        return altorfs

    for event in ['PRF','SCR','uORF','dORF','ORF']:
        if event in altorfs_idx:

            event_idx = set(altorfs_idx[event])
            path_idx = [i for i,s in enumerate(path) if s in event_idx]

            ranges = list(get_ranges(path_idx))
            for range in ranges:

                # Print type of PRF and interval of positions
                if event=='PRF':
                    if range[1]-range[0]==0: # +1 PRF
                        altorfs.append('+1 %s (%i-%i)' % (event, range[0]+1, range[1]+1))
                    elif range[1]-range[0]==1: # -1 PRF
                        altorfs.append('-1 %s (%i-%i)' % (event, range[0]+1, range[1]+1))

                # Print event and interval of positions
                else:
                    altorfs.append('%s (%i-%i)' % (event, range[0]+1, range[1]+1))

    return altorfs


def readable_sequence(x,
                      path,
                      states_idx):
    """
    Generates a human-readable representation of the predicted state path for
    the given nucleotide sequence:
        [] = bounds of main ORFs, uORFs, dORFs
         | = SCR site
        () = PRF skipped nucleotide (+1 PRF) or repeated nucleotides (-1 PRF)
         . = UTR nucleotide
    Args:
        x: the nucleotide sequence
        path: the state path
        states_idx: dictionary mapping each general state type to its
            state index (1-indexed)
    Returns:
        readable_seq: the state path in a readable format
    Raises:

    """
    readable_seq = ''

    for i in range(len(x)):

        # UTR
        if path[i]==states_idx['5UTR'] or path[i]==states_idx['3UTR']:
            readable_seq += '.'

        # uORF start codon (1st nt)
        elif path[i] in range(states_idx['start_uORF'][0],
                              states_idx['start_uORF'][-1]+1,
                              3):
            readable_seq += '[' + x[i]

        # dORF start codon (1st nt)
        elif path[i] in range(states_idx['start_dORF'][0],
                              states_idx['start_dORF'][-1]+1,
                              3):
            readable_seq += '[' + x[i]

        # ORF start codon (1st nt)
        elif path[i] in range(states_idx['start'][0],
                              states_idx['start'][-1]+1,
                              3):
            readable_seq += '[' + x[i]

        # uORF stop codon (3rd nt)
        elif path[i] in range(states_idx['stop_uORF'][0]+2,
                              states_idx['stop_uORF'][-1]+1,
                              3):
            readable_seq += x[i] + ']'

        # dORF stop codon (3rd nt)
        elif path[i] in range(states_idx['stop_dORF'][0]+2,
                              states_idx['stop_dORF'][-1]+1,
                              3):
            readable_seq += x[i] + ']'

        # ORF stop codon (3rd nt)
        elif path[i] in range(states_idx['stop'][0]+2,
                              states_idx['stop'][-1]+1,
                              3):
            readable_seq += x[i] + ']'

        # SCR
        elif (path[i] in range(states_idx['stop_SCR'][0]+2,
                               states_idx['stop_SCR'][-1]+1,
                               3)):
            readable_seq += x[i] + '|'

        # PRF
        elif path[i]==states_idx['X1'] or path[i]==states_idx['X2']:
            readable_seq += '(' + x[i] + ')'

        # Sense codon
        else:
            readable_seq += x[i]

    return readable_seq


def get_orfs(readable_seq):
    """
    Extracts the nucleotide sequence for each ORF from a prediction result
    Args:
        readable_seq: the state path in a readable format
    Returns:
        seqs: nucleotide sequences for each ORF
    Raises:

    """
    # Get the sequence and info for each ORF
    seqs = []
    for x in readable_seq:
        if x=='.':
            continue
        elif x=='[':
            seq = ''
        elif x==']':
            seqs.append(seq)
        elif x=='(':
            continue
        elif x==')': # Skip the frameshift nucleotide
            seq = seq[:-1] # (really it could be repeated for a -1 or -2 PRF)
        elif x=='|': # Skip the readthrough stop
            seq = seq[:-3] # (really it could be replaced by a near-cognate codon)
        else:
            seq += x

    return seqs


def translate(x, amino_acids):
    """
    Translate a nucleotide sequence to protein
    Args:
        x: nucleotide sequence
        amino_acids: dict mapping each codon to an amino acid (codon table)
    Returns:
        y: protein sequence
    Raises:

    """
    # Translate the nt sequence using the codon table
    L = len(x)
    y = [amino_acids[''.join(x[i:i+3])] for i in range(0,L-2,3)]

    # If there is no stop codon return an empty sequence
    if '*' not in y:
        logging.warning('Couldn\'t translate sequence (missing stop codon): %s' % x)
        return ''

    # If there is a stop codon, truncate the protein sequence
    else:
        # Find the first stop codon
        idx_stop = y.index('*')

        # Terminate the protein sequence right before the first stop codon
        return ''.join(y[:idx_stop])


def trace_dump(x,
               r,
               path,
               states,
               log_emit,
               log_trans):
    """
    Prints the path score at each position of the sequence
    Args:
        x: the nucleotide sequence
        r: the rho sequence
        path: the state path
        states: list of all states in the model
        log_emit: the total emission probabilities at each position for
            this sequence (combined nucleotide and rho emissions)
        log_trans: log transition probability matrix
    Returns:

    Raises:

    """
    # Length of observed sequence
    L = len(path)

    if L==0:
        return np.nan

    # Print the header
    print('{:5s} {:5s} {:15s} {:5s} {:10s} {:15s} {:15s} {:15s}'\
          .format('pos',
                  'nt',
                  'rho',
                  'state_idx',
                  'state',
                  'log_transition',
                  'log_emission',
                  'log_posterior'))

    # Add the probability of the initial transition
    # from the begin state (in log space)
    P = log_trans[0, path[0]+1]

    # For each position in the sequence 0..L-1
    # add the probability of the transition to this residue and
    # the emission for this residue (in log space)
    for i in range(0, L-1):
        P += log_emit[i, path[i]] + log_trans[path[i]+1, path[i+1]+1]
        print('{:5.0f} {:5s} {:15.10f} {:5.0f} {:10s} {:15.10f} {:15.10f} {:15.10f}'\
              .format(i,
                      x[i],
                      r[i],
                      path[i],
                      states[path[i]],
                      log_trans[path[i]+1, path[i+1]+1],
                      log_emit[i, path[i]],
                      P))

    # Add the probability of the final emission and
    # the transition to the end state (in log space)
    P += log_emit[L-1, path[L-1]] + log_trans[path[L-1]+1, 0]
    print('{:5.0f} {:5s} {:15.10f} {:5.0f} {:10s} {:15.10f} {:15.10f} {:15.10f}'\
          .format(L-1,
                  x[L-1],
                  r[L-1],
                  path[L-1],
                  states[path[L-1]],
                  log_trans[path[i]+1,
                  path[i+1]+1],
                  log_emit[i, path[i]],
                  P))


def write_results_fasta(results,
                        parameters_h1,
                        fn_file,
                        fa_file,
                        chunk=60):
    """
    Writes a fasta file of predicted nucleotide and protein sequences
    Args:
        results: array containing all results of running the model
        parameters_h1: list of all parameter matrices for the current model
        file: path to output fasta file
        fn_file: path to output nucleotide sequence fasta file
        fa_file: path to output protein sequence fasta file
        chunk: max line length for sequence entries, add line breaks to
            maintain this max length
    Returns:

    Raises:
        IOError: an error occurred accessing the file
    """
    # Unpack the parameters
    states, states_idx, altorfs_idx, \
    p_emit_nt, p_emit_rho, rho_ranges, p_trans, p_trans_rev, \
    prev_states_options, next_states_options = parameters_h1

    # Open the fasta files for writing
    try:
        fn = open(fn_file,'w')
    except IOError:
        logging.error('Could not open fasta file %s for writing' % file)

    try:
        fa = open(fa_file,'w')
    except IOError:
        logging.error('Could not open fasta file %s for writing' % file)

    else:

        with fn, fa:

            for result in results:

                # Unpack the transcript attributes
                transcript_id, transcript_name, coverage, x, r, \
                          path_h2, path_h1, path_h0, \
                          path_score_h2, path_score_h1, path_score_h0, \
                          log_odds_score = result

                # Write the nucleotide sequence for this transcript
                id = '%s (%s)' % (transcript_id, transcript_name)
                seq = readable_sequence(x, path_h1, states_idx)

                # Write the seq id
                fn.write('>%s orfeus prediction\n' % id)

                # Write the sequence with a fixed max line length
                i = 0
                while i < len(seq)//chunk:
                    fn.write('%s\n' % seq[i*chunk:(i+1)*chunk])
                    i += 1
                fn.write('%s\n' % seq[i*chunk:])

                # Write the protein sequence(s) for this transcript
                xs = get_orfs(seq)
                k = 0
                for x in xs:
                    y = translate(x, AMINO_ACIDS)

                    # Write the seq id
                    fa.write('>%s orfeus prediction %i\n' % (id,k))
                    k += 1

                    # Write the sequence with a fixed max line length
                    i = 0
                    while i < len(y)//chunk:
                        fa.write('%s\n' % y[i*chunk:(i+1)*chunk])
                        i += 1
                    fa.write('%s\n' % y[i*chunk:])


def write_results_table(results,
                        parameters_h1,
                        file):
    """
    Writes a text file of predicted ORF and altORF events
    Args:
        results: array containing all results of running the model
        parameters_h1: list of all parameter matrices for the current model
        file: path to output fasta file
    Returns:

    Raises:
        IOError: an error occurred accessing the file
    """
    # Unpack the parameters
    states, states_idx, altorfs_idx, \
    p_emit_nt, p_emit_rho, rho_ranges, p_trans, p_trans_rev, \
    prev_states_options, next_states_options = parameters_h1

    # Open the fasta file for writing
    try:
        f = open(file,'w')

    except IOError:
        logging.error('Could not open output file %s for writing' % file)

    else:

        # Write the header
        f.write('{:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s}\n' \
                 .format('transcript_id',
                         'transcript_name',
                         'mean_coverage',
                         'predicted=annotated',
                         'annotated_events',
                         'predicted_events',
                         'log_odds_score',
                         'annotated_path_score',
                         'predicted_path_score',
                         'null_path_score'))

        # Sort by log-odds score
        results = sorted(results, key=lambda x: x[-1], reverse=True)

        with f:

            for result in results:

                # Unpack the transcript attributes
                transcript_id, transcript_name, coverage, x, r, \
                          path_h2, path_h1, path_h0, \
                          path_score_h2, path_score_h1, path_score_h0, \
                          log_odds_score = result

                # Get the events
                events_h2 = ', '.join(get_events(path_h2, altorfs_idx))
                events_h1 = ', '.join(get_events(path_h1, altorfs_idx))

                # Check whether annotated ORF = predicted ORF
                if path_score_h2==path_score_h1: canonical = 'true'
                else: canonical = 'false'

                # Write the results
                f.write('{:20s} {:20s} {:20.4f} {:20s} {:20s} {:20s} {:20.3f} {:20.3f} {:20.3f} {:20.3f}\n' \
                        .format(transcript_id,
                                transcript_name,
                                coverage,
                                canonical,
                                events_h2,
                                events_h1,
                                log_odds_score,
                                path_score_h2,
                                path_score_h1,
                                path_score_h0))
