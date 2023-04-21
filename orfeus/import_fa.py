#!/usr/bin/env python

# -----------------------------------------------------------------------------
# FASTA READER/WRITER
# Author: Mary Richardson
# Date: 2022.05.25
# -----------------------------------------------------------------------------

import gzip
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def reverse_complement(seq):
    """
    Generates the reverse complement of a nucleotide sequence
    Args:
        seq: string representing a nucleotide sequence
    Returns:
        rev_comp: string representing the reverse complement nucleotide
            sequence
    Raises:
        KeyError: an error occurred finding the reverse complement
    """
    # Define the complement of each nucleotide
    complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'R':'Y', 'Y':'R'}

    # Complement the nucleotide sequence
    try:
        comp = [complement[nt] for nt in seq]
    except KeyError:
        logging.error(('Invalid nucleotide in transcript sequence, '
            'cannot find reverse complement for sequence %s') % seq)

    # Reverse the nucleotide sequence
    rev_comp = comp[::-1]

    return ''.join(rev_comp)


def read_fasta(file, rev_comp=False):
    """
    Fetches sequence(s) from a fasta file
    Args:
        file: path to a fasta file
        rev_comp: (optional) whether to calculate the reverse complement
    Returns:
        seqs: dictionary mapping each sequence ID to its sequence
    Raises:
        IOError: an error occurred accessing the file
        IndexError, TypeError, KeyError: an error occurred reading the file
    """
    # Open the fasta file for reading
    try:
        if file.endswith('.gz'): f = gzip.open(file,'r')
        else: f = open(file,'r')

    except IOError:
        logging.error('Could not open fasta file %s for reading' % file)

    else:

        # Prepare a dictionary to store the seqs
        seqs = {}

        with f:
            for line in f:

                try:
                    # If this is a seq id line, start a new sequence
                    if line[0]=='>':
                        id = line.strip()[1:].split()[0]
                        seqs[id] = ''

                    # Otherwise, append to the current sequence
                    else:

                        # For negative strand, sequence is reverse complement
                        if rev_comp:
                            seqs[id] = reverse_complement(line.strip()) + seqs[id]

                        # For positive strand, sequence is unchanged
                        else:
                            seqs[id] = seqs[id] + line.strip()

                except:
                    logging.error('File %s is not in expected fasta format' \
                          % file)

    return seqs


def write_fasta(seqs: dict, file: str, chunk=60):
    """
    Writes a fasta file of sequences
    Args:
        seqs: dictionary mapping sequence IDs to sequences
        file: path to output fasta file
        chunk: max line length for sequence entries, add line breaks to
            maintain this max length
    Returns:

    Raises:
        IOError: an error occurred accessing the file
    """
    # Open the fasta file for writing
    try:
        f = open(file,'w')

    except IOError:
        logging.error('Could not open fasta file %s for writing' % file)

    else:

        with f:
            # Write all sequences
            for id in seqs:

                # Write the seq id
                f.write('>%s\n' % id)

                # Write the sequence with a fixed max line length
                i = 0
                while i < len(seqs[id])//chunk:
                    f.write('%s\n' % seqs[id][i*chunk:(i+1)*chunk])
                    i += 1
                f.write('%s\n' % seqs[id][i*chunk:])


def fasta_to_df(file, rev_comp=False):
    """
    Generate a DataFrame containing the nucleotide sequence at each
    position of each transcript
    Args:
        file: path to a fasta file containing a nucleotide sequence
            for each transcript
        rev_comp: (optional) whether to calculate the reverse complement
    Returns:
        pandas DataFrame with the transcript nucleotide sequences
    Raises:
        IOError: an error occurred accessing the transcriptome file
        IndexError, TypeError, KeyError: an error occurred reading the
            transcriptome file
    """
    # Create a dict of transcript sequences
    seqs = read_fasta(file, rev_comp)

    # Convert to a dataframe
    df = pd.DataFrame(seqs.items(), columns=['seq_id','seq'])
    df = df.assign(seq = df['seq'].str.split('')).explode('seq')
    df = df[df['seq'].astype(bool)]
    df['seq_pos'] = df.groupby('seq_id').cumcount()+1

    return df
