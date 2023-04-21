#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Coverage threshold
# Author: Mary Richardson
# Date: 2020.06.04
# -----------------------------------------------------------------------------

import os
from coverage_sensitivity import sensitivity_simulation
from coverage_specificity import specificity_simulation
from coverage_accuracy import accuracy_simulation
from coverage_count import coverage_count

NUCLEOTIDES = ['A','C','G','T']


def coverage_threshold_simulation(data_df,
                                  parameters,
                                  coverages,
                                  iters,
                                  window,
                                  outdir,
                                  threads):

    # For each event type
    for event in ['uORFdORF', 'SCR', 'pPRF', 'mPRF']:

        # Sensitivity simulation
        sens_file = os.path.join(outdir, 'sensitivity_' + event + '.txt')
        sensitivity = sensitivity_simulation(parameters,
                                             NUCLEOTIDES,
                                             coverages,
                                             iters,
                                             window,
                                             event,
                                             threads,
                                             sens_file)

        # Specificity simulation
        spec_file = os.path.join(outdir, 'specificity_' + event + '.txt')
        specificity = specificity_simulation(parameters,
                                             NUCLEOTIDES,
                                             coverages,
                                             iters,
                                             event,
                                             threads,
                                             spec_file)

    # Canonical ORF sensitivity simulation
    sens_file = os.path.join(outdir, 'sensitivity_ORF.txt')
    sensitivity = sensitivity_simulation(parameters,
                                         NUCLEOTIDES,
                                         coverages,
                                         iters,
                                         window,
                                         'ORF',
                                         threads,
                                         sens_file)

#    # Accuracy simulation
#    for event in ['ORF', 'any']:
#        acc_file = os.path.join(outdir, 'accuracy_' + event + '.txt')
#        accuracy = accuracy_simulation(parameters,
#                                       NUCLEOTIDES,
#                                       coverages,
#                                       iters,
#                                       event,
#                                       threads,
#                                       acc_file)

    # Coverage count
    cov_file = os.path.join(outdir, 'coverage_counts.txt')
    counts = coverage_count(data_df,
                            coverages,
                            cov_file)
