#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Data processing
# Author: Mary Richardson
# Date: 2020.01.08, last update 2023.03.01
# -----------------------------------------------------------------------------

import logging
import numpy as np
import pandas as pd

from import_wig import *
from import_bg import *
from import_fa import *
from import_gtf import *

logger = logging.getLogger(__name__)

# Ignore divide by zero error for log
np.seterr(divide='ignore')



def read_counts_data(file, strand):
    """
    Generate a DataFrame containing the ribo-seq read counts at each
    position of the genome
    Args:
        file: path to a read counts Wiggle or BedGraph file
        strand: strand of reads specified in the input file ('+' or '-)'
    Returns:
        df: pandas DataFrame with the read counts and strand
    Raises:
        IOError: an error occurred accessing the read counts file
        IndexError, TypeError, KeyError: an error occurred reading the
            read counts file
    """
    # Parse Wiggle or BedGraph file
    try: # Try to parse as bedgraph
        df = bg_to_df(file)
    except:
        try: # Try to parse as wiggle
            df = wig_to_df(file)
        except:
            logging.error(('Could not parse read counts file %s. '
                    'Make sure file is a Wiggle (.wig) or BedGraph (.bg) file.')
                    % file)
            raise

    df['strand'] = strand
    return df


def seq_data(file, strand):
    """
    Generate a DataFrame containing the nucleotide sequence at each
    position of the genome or transcriptome
    Args:
        file: path to a fasta file
        strand: strand of seqs specified in the input file ('+' or '-)'
    Returns:
        df: pandas DataFrame with the seq and strand
    Raises:
        IOError: an error occurred accessing the seq file
        IndexError, TypeError, KeyError: an error occurred reading the
            seq counts file
    """
    # Parse fasta file
    try:
        if strand=='+':
            df = fasta_to_df(file, rev_comp=False)
        elif strand=='-':
            df = fasta_to_df(file, rev_comp=True)

    except:
        logging.error(('Could not parse read counts file %s. '
                'Make sure file is a Fasta file.')
                % file)
        raise

    df['strand'] = strand
    return df


def feature_label(feature):
    """
    Updates feature labels to shortened feature names
    Args:
        feature: current feature name
    Returns:
        feature: updated feature name
    """
    if feature=='five_prime_utr':
            return '5UTR'
    elif feature=='exon':
            return 'ORF'
    elif feature=='three_prime_utr':
            return '3UTR'
    else: return feature


def annotations_data(file):
    """
    Parse a GTF/GFF file to extract all annotated transcripts
    Args:
        file: path to an annotations gff/gtf file
    Returns:
        df: pandas DataFrame of all annotated transcripts
    Raises:
        IOError: an error occurred accessing the file
    """
    # Parse input info
    try:
        df = gtf_to_df(file)
    except:
        logging.error(('Could not parse annotations file %s. '
                'Make sure file is a GTF or GFF file.')
                % file)
        raise

    # Filter down to columns of interest
    df = df[['chrom',
             'transcript_id',
             'transcript_name',
             'transcript_biotype',
             'feature',
             'strand',
             'start',
             'end']]

    # Filter down to features of interest (protein-coding UTRs and ORFs)
    df = df[df['transcript_biotype'].isin(['protein_coding'])]
    df = df[df['feature'].isin(['five_prime_utr',
                                'exon',
                                'three_prime_utr'])]
    df['feature'] = pd.Categorical(df['feature'].map(feature_label),
                                   categories=['5UTR', 'ORF', '3UTR'],
                                   ordered=True)

    # Sort by chrom then transcript then feature
    df.sort_values(by=['chrom',
                       'transcript_id',
                       'feature'],
                   inplace=True)

    # Add chromosome position
    df['chrom_pos'] = df.apply(lambda x: list(range(int(x['start']),
                                                   int(x['end'])+1)),
                                                   1)
    df = df.explode('chrom_pos')
    df = df.drop(columns=['start', 'end'])

    # Add transcript position
    df['transcript_pos'] = df.groupby('transcript_id') \
                             .cumcount()+1

    # Add feature position
    df['feature_pos'] = df.groupby(['transcript_id','feature']) \
                          .cumcount()+1

    # Add transcript length
    df['transcript_len'] = df.groupby(['transcript_id']) \
                                      ['transcript_pos'].transform(len)

    # Add feature length
    df['feature_len'] = df.groupby(['transcript_id','feature']) \
                                   ['feature_pos'].transform(len)

    if df.empty:
        logging.error(('Processed data is empty. Check formatting of GTF/GFF'
                'file %s') % file)
        raise

    return df


def process_data(plus_file,
                 minus_file,
                 seqs_file,
                 annotations_file,
                 skip=[]):
    """
    Merge riboseq data and annotations data into one DataFrame
    Args:
        plus_file: plus strand read counts bg/wig file
        minus_file: minus strand read counts bg/wig file
        seqs_file: genome or transcriptome nucleotide sequence fasta file
        annotations_file: annotations gff/gtf file
        skip: (optional) names of frameshifted genes to exclude
            when adding codon/state labels (recommended to prevent out-of-frame
            and partial codons from being added)
    Returns:
        data_df: pandas DataFrame of all annotated transcripts
    Raises:
        IOError: an error occurred accessing the file
        ValueError: could not merge sequence and annotation data
    """
    # Annotations
    annot_df = annotations_data(annotations_file)

    # Read counts
    reads_plus_df = read_counts_data(plus_file, '+')
    reads_minus_df = read_counts_data(minus_file, '-')
    reads_df = pd.concat([reads_plus_df, reads_minus_df], axis=0)

    # Sequences
    seqs_plus_df = seq_data(seqs_file, '+')
    seqs_minus_df = seq_data(seqs_file, '-')
    seqs_df = pd.concat([seqs_plus_df, seqs_minus_df], axis=0)

    logging.info('annotations: \n{}'.format(annot_df.head()))
    logging.info('reads: \n{}'.format(reads_df.head()))
    logging.info('genome sequence: \n{}'.format(seqs_df.head()))

    # Merge read counts and annotations
    df = pd.merge(annot_df, reads_df,
                  on=['chrom','chrom_pos','strand'],
                  how='left')

    # Fill any missing read read_counts positions with zero
    df['reads'].fillna(0, inplace=True)

    # Get the IDs for merging
    seq_ids = set(seqs_df['seq_id'].unique())
    annot_chroms = set(annot_df['chrom'].unique())
    annot_txs = set(annot_df['transcript_id'].unique())

    # Delete old dataframes
    del annot_df
    del reads_df

    # Merge sequence by chromosomes
    if seq_ids.issubset(annot_chroms) or annot_chroms.issubset(seq_ids):

        # Add chromosome length
        chrom_len = seqs_df['seq_id'].value_counts().to_dict()
        df['chrom_len'] = df['chrom'].map(chrom_len)//2

        # Correct the genome position to account for strand
        df['chrom_pos'] = np.where(df['strand']=='-', \
                          (df['chrom_len']-df['chrom_pos']+1), \
                          df['chrom_pos'])

        # Correct the transcript position to account for strand
        df['feature_pos'] = np.where(df['strand']=='-', \
                            (df['feature_len']-df['feature_pos']+1), \
                            df['feature_pos'])

        # Re-sort by corrected feature position
        df.sort_values(by=['chrom',
                           'transcript_id',
                           'feature',
                           'feature_pos'],
                       inplace=True)

        # Correct transcript position
        df['transcript_pos'] = df.groupby('transcript_id') \
                                 .cumcount()+1

        df = pd.merge(df, seqs_df,
                      left_on=['chrom','chrom_pos','strand'],
                      right_on=['seq_id','seq_pos','strand'],
                      how='left') \
               .drop(columns=['seq_id','seq_pos'])

        # Delete old dataframe
        del seqs_df

    # Or merge sequence by transcripts
    elif seq_ids.issubset(annot_txs) or annot_txs.issubset(seq_ids):
        df = pd.merge(df, seqs_df,
                      left_on=['transcript_id','transcript_pos','strand'],
                      right_on=['seq_id','seq_pos','strand'],
                      how='left') \
               .drop(columns=['seq_id','seq_pos'])

        # Delete old dataframe
        del seqs_df

    else:
        logging.error(('Could not merge sequence and annotation data. '
                'Make sure sequence IDs (chromosomes or transcripts) are the '
                'same in the sequence and annotation files.'))

    # Add nt position (for ORF only)
    df['nt_pos'] = np.where(df['feature'] == 'ORF',
                           ((df['feature_pos']+2)%3+1),
                           '')

    # Add codon position (for ORF only)
    df['codon_pos'] = np.where(df['feature'] == 'ORF',
                              (df['feature_pos']-1)//3,
                              -1)

    # Add codon type (for ORF only)
    conditions = [(df['feature'] == 'ORF') & \
                  (df['codon_pos'] == 0),
                  (df['feature'] == 'ORF') & \
                  (df['codon_pos'] == df['feature_len']//3-1)]
    choices = ['START', 'STOP']
    df['codon_type'] = np.select(conditions,
                                 choices,
                                 default=df['feature'])

    # Normalize riboseq reads by the mean riboseq reads per transcript
    df['mean_reads'] = df.groupby('transcript_id') \
                                   ['reads'].transform('mean')
    #df['mean_orf_reads'] = df.groupby(['gene_id','feature']) \
    #                                    ['reads'].transform('mean')
    #df.loc[data_df['feature']!='ORF', 'mean_orf_reads'] = np.nan
    #df['mean_orf_reads'] = df.groupby('gene_id') \
    #                                    ['mean_orf_reads'].transform('first')

    # Normalize by mean reads over transcript
    df['norm_reads'] = (df['reads']/df['mean_reads']) \
                                .replace(np.nan, 0)
    df['log_reads'] = np.log(df['norm_reads'] + 1)

    # Find the codon sequence at each genome position that contains a ORF
    codons_df = df[df['feature']=='ORF'].copy()

    # Ignore known frameshifted genes,
    # since the codon sequence will be incorrectly inferred
    codons_df = codons_df[~codons_df['transcript_name'].isin(skip)]

    # Codons are just every three nucleotides
    codons_df['codon_seq'] = codons_df.groupby(['transcript_id','codon_pos']) \
                                      ['seq'].transform(lambda x: ''.join(x))

    # Add the codon sequence back to the final dataframe
    df = df.merge(codons_df[['transcript_id',
                             'transcript_pos',
                             'codon_seq']],
                  on=['transcript_id','transcript_pos'],
                  how='left')

    # Delete old dataframe
    del codons_df

    # Add the state to the dataframe
    conditions = [(df['codon_type'] == 'START'),
                  (df['codon_type'] == 'STOP'),
                  (df['codon_type'] == 'ORF'),
                  (df['feature'] == '5UTR') | (df['feature'] == '3UTR')]
    choices = [df['codon_seq'] + df['nt_pos'].astype(str) + '#',
               df['codon_seq'] + df['nt_pos'].astype(str) + '*',
               df['codon_seq'] + df['nt_pos'].astype(str),
               df['feature']]
    df['state'] = np.select(conditions, choices, default='')

    # Clean up the columns
    df.drop(columns=['chrom_pos', 'codon_pos', 'chrom_len'],
                     inplace=True)
    dtype0 = {'transcript_id': 'str',
             'transcript_name': 'str',
             'transcript_pos': 'int64',
             'transcript_len': 'int64',
             'feature': 'str',
             'feature_pos': 'int64',
             'feature_len': 'int64',
             'seq': 'str',
             'codon_seq': 'str',
             'codon_type': 'str',
             'nt_pos': 'str',
             'state': 'str',
             'reads': 'float64',
             'mean_reads': 'float64',
             'norm_reads': 'float64',
             'log_reads': 'float64'}
    df = df.astype(dtype0, copy=False)

    return df


def read_data_file(file):
    """
    Read processed data from a file into a DataFrame
    Args:
        file: path of csv file containing exported processed data
    Returns:
        data_df: pandas DataFrame of all annotated transcripts
    Raises:

    """
    # Import the data (specify dtypes where necessary)
    dtype0 = {'transcript_id': 'str',
             'transcript_name': 'str',
             'transcript_pos': 'int64',
             'transcript_len': 'int64',
             'feature': 'str',
             'feature_pos': 'int64',
             'feature_len': 'int64',
             'seq': 'str',
             'codon_seq': 'str',
             'codon_type': 'str',
             'nt_pos': 'str',
             'state': 'str',
             'reads': 'float64',
             'mean_reads': 'float64',
             'norm_reads': 'float64',
             'log_reads': 'float64'}

    data_df = pd.read_csv(file,
                          sep='\t',
                          index_col=0,
                          compression='gzip',
                          dtype=dtype0)

    return data_df
