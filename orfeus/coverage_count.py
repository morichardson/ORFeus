#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Gene coverage
# Author: Mary Richardson
# Date: 2023.03.20
# -----------------------------------------------------------------------------

import logging
import numpy as np

logger = logging.getLogger(__name__)



def coverage_count(data_df,
                   coverages,
                   file):
    """
    Counts the number of transcripts in the dataset with at least the each
        specified coverage level
    Args:
        data_df: pandas DataFrame of all annotated transcripts
        coverages: mean coverages per transcript (reads/nt)
        file: output file for results
    Returns:

    Raises:

    """
    try:
        with open(file,'w') as f:
            f.write('Counting transcripts with at least each coverage level \n')
    except IOError:
        logging.error('Could not open coverage ouput file %s' % file)

    # Calculate the transcripts for each coverage
    coverages = np.insert(coverages,0,0)

    counts = {}
    for coverage in coverages:

        # Count the number of genes with at least the specified coverage
        count = len(data_df[data_df['mean_reads']>=coverage] \
                           ['transcript_id'].unique())

        counts[coverage] = count
        try:
            with open(file,'a') as f:
                f.write('%f coverage transcript count: %f \n' % (coverage, count))
        except IOError:
            logging.error('Could not open coverage ouput file %s' % file)

    return counts
