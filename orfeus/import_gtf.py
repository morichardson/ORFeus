#!/usr/bin/env python

# -----------------------------------------------------------------------------
# GFF/GTF READER
# Author: Adapted from code by Kamil Slowikowski
# available from https://gist.github.com/slowkow/8101481
# Date: 2022.11.02
# -----------------------------------------------------------------------------

import re
import gzip
import logging
import pandas as pd
from collections import defaultdict

logger = logging.getLogger(__name__)


GTF_HEADER  = ['chrom', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def gtf_to_df(file):
    """
    Open an optionally gzipped GTF file and return a pandas DataFrame
    GTF file format: ensembl.org/info/website/upload/gff.html
    Args:
        file: path to an annotations gff/gtf file
    Returns:
        pandas DataFrame with the annotations
    Raises:
        IOError: an error occurred accessing the annotations file
        IndexError, TypeError, KeyError: an error occurred reading the
            annotations file
    """
    # Each column is a list stored as a value in this dict
    annotations = defaultdict(list)

    for i, line in enumerate(lines(file)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines
            if key not in annotations:
                annotations[key] = [None] * i

        # Ensure this row has some value for each column
        for key in annotations.keys():
            annotations[key].append(line.get(key, None))

    return pd.DataFrame(annotations)


def lines(file):
    # Open the annotations gtf file for reading
    try:
        if file.endswith('.gz'): f = gzip.open(file,'r')
        else: f = open(file,'r')

    except IOError:
        logging.error('Could not open annotations file %s' % file)

    else:
        with f:
            try:
                for line in f:

                    # Skip header lines
                    if line.startswith('#'):
                        continue

                    # Parse annotation data lines
                    else:
                        yield parse(line)

            except (IndexError, TypeError, KeyError):
                logging.error(('Could not parse annotations file %s.') % file)
                raise SyntaxError


def parse(line):
    annotations = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        annotations[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;..."
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value"
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value"
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value
        if value:
            annotations[key] = _get_value(value)

    return annotations


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes
    value = value.strip('"\'')

    # Return a list if the value has a comma
    if ',' in value:
        value = re.split(R_COMMA, value)

    # These values are equivalent to None
    elif value in ['', '.', 'NA']:
        return None

    return value
