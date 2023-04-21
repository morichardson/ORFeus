#!/usr/bin/env python

# -----------------------------------------------------------------------------
# WIG READER
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

R_KEYVALUE = re.compile(r'(\s*=\s*)')
R_CHROM  = re.compile(r'(\s*chr\s*)')


def wig_to_df(file):
    """
    Open a Wiggle file and return a pandas DataFrame
    Wiggle file format: ensembl.org/info/website/upload/wig.html
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

                    # Set variables from info lines
                    elif line.startswith('fixedStep'):
                        set_variables(line)
                    elif line.startswith('variableStep'):
                        set_variables(line)

                    # Parse read count data lines
                    else:
                        fields = line.rstrip().split()

                        step_type = PARAMETERS['step_type']
                        chrom = PARAMETERS['chrom']
                        start = int(PARAMETERS['start'])
                        step = int(PARAMETERS['step'])
                        span = int(PARAMETERS['span'])

                        if step_type=='fixedStep':
                            for pos in range(start,start+span):
                                reads = float(fields[0])
                                yield set_read_counts(chrom, pos, reads)
                            PARAMETERS['start'] = start + step

                        elif PARAMETERS['step_type']=='variableStep':
                            start = int(fields[0])
                            reads = float(fields[1])
                            for pos in range(start,start+span):
                                yield set_read_counts(chrom, pos, reads)

            except (IndexError, TypeError, KeyError):
                logging.warning('Could not parse %s as Wiggle file.' % file)
                raise SyntaxError


def set_variables(line):
    fields = line.rstrip().split()

    PARAMETERS['step_type'] = fields[0]

    # Field consists of "key1=value key2=value ..."
    for i, info in enumerate(fields[1:]):
        # It should be key=value
        key, _, value = re.split(R_KEYVALUE, info, 1)
        # Ignore the field if there is no value
        if value:
            PARAMETERS[key] = _get_value(value)


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
