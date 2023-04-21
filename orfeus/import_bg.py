#!/usr/bin/env python

# -----------------------------------------------------------------------------
# BEDGRAPH READER
# Author: Mary Richardson
# Date: 2022.11.02
# -----------------------------------------------------------------------------

import re
import gzip
import logging
import pandas as pd
from collections import defaultdict

logger = logging.getLogger(__name__)


READS_HEADER  = ['chrom', 'chrom_pos', 'reads']
PARAMETERS = {'step_type':'', 'chrom':'', 'start':0, 'step':0, 'span':1}

R_CHROM  = re.compile(r'(\s*chr\s*)')


def bg_to_df(file):
    """
    Open a BedGraph file and return a pandas DataFrame
    BedGraph file format: ensembl.org/info/website/upload/bed.html
    Args:
        file: path to a read counts wig file
    Returns:
        pandas DataFrame with the read counts
    Raises:
        IOError: an error occurred accessing the read counts file
        IndexError, TypeError, KeyError: an error occurred reading the
            read counts file
    """
    # Each column is a list stored as a value in this dict
    reads = defaultdict(list)

    for i, line in enumerate(lines(file)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines
            if key not in reads:
                reads[key] = [None] * i

        # Ensure this row has some value for each column
        for key in reads.keys():
            reads[key].append(line.get(key, None))

    return pd.DataFrame(reads)



def lines(file):
    # Open the read counts wig file for reading
    try:
        if file.endswith('.gz'): f = gzip.open(file,'r')
        else: f = open(file,'r')

    except IOError:
        logging.error('Could not open read counts file %s' % file)

    else:
        with f:
            try:
                for line in f:

                    # Skip header lines
                    if line.startswith('#'):
                        continue
                    elif line.startswith('track'):
                        continue

                    # Parse read count data lines
                    else:
                        fields = line.rstrip().split()

                        chrom = _get_value(fields[0])
                        start = int(fields[1])
                        end = int(fields[2])
                        reads = float(fields[3])

                        for pos in range(start,end):
                            yield set_read_counts(chrom, pos, reads)

            except (IndexError, TypeError, KeyError):
                logging.warning('Could not parse %s as BedGraph file.' % file)
                raise SyntaxError


def set_read_counts(chrom, pos, reads):
    read_counts = {}

    read_counts['chrom'] = chrom
    read_counts['chrom_pos'] = pos
    read_counts['reads'] = reads

    return read_counts



def _get_value(value):
    if not value:
        return None

    # Return chromosome number
    if 'chr' in value:
        value = re.split(R_CHROM, value)[1]

    # These values are equivalent to None
    elif value in ['', '.', 'NA']:
        return None

    return value
