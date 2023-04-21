#!/usr/bin/env python

# -----------------------------------------------------------------------------
# HMM plotting functions
# Author: Mary Richardson
# Date: 2021.01.05
# -----------------------------------------------------------------------------

import numpy as np
import scipy.stats as ss
import itertools as it
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter as fmt

from hmm_run import hmm_single, get_annotated_path
from output_text import readable_sequence, trace_dump

# Customize plotting style
colors = ['#6F6F6F',
          '#DDDDDD',
          '#E0BA4D',
          '#9ABC61',
          '#E3995A',
          '#5385AA',
          '#C55A56',
          '#A472A5']
basics = ['#5385AA',
          '#C55A56',
          '#6F6F6F']


def calculate_metagene_signal(data_df, w):
    """
    Calculate the mean signal upstream/downstream of the start and stop codon
    Args:
        data_df: processed data
        w: window of nucleotides upstream/downstream of the start or stop codon
            to return
    Returns:
        signal_start: mean signal upstream/downstream of the first ORF nt
        signal_stop: mean signal upstream/downstream of the last ORF nt
    Raises:

    """
    # Default zeros
    signal_5UTR = np.zeros(w)
    signal_3UTR = np.zeros(w)

    # Group the data by transcript then feature
    df = data_df.groupby(['transcript_id', 'feature']) \
                .agg(list).reset_index()

    # Get the reads by feature
    data_ORF = df[df['feature']=='ORF']['norm_reads'].tolist()
    data_5UTR = df[df['feature']=='5UTR']['norm_reads'].tolist()
    data_3UTR = df[df['feature']=='3UTR']['norm_reads'].tolist()

    # Calculate the total number of transcripts
    tx_count = len(data_ORF)

    # START codon signal
    if len(data_5UTR) > 0:
        signal_5UTR = np.sum(np.array([([0]*w + r)[-w:] for r in data_5UTR]), \
                             axis=0)
    signal_ORF = np.sum(np.array([(r + [0]*w)[:w] for r in data_ORF]), \
                        axis=0)

    # Aggregate START signal
    signal_start = np.concatenate((signal_5UTR, signal_ORF))
    signal_start /= tx_count

    # STOP codon
    signal_ORF = np.sum(np.array([([0]*w + r)[-w:] for r in data_ORF]), \
                       axis=0)
    if len(data_3UTR) > 0:
        signal_3UTR = np.sum(np.array([(r + [0]*w)[:w] for r in data_3UTR]), \
                            axis=0)

    # Aggregate STOP signal
    signal_stop = np.concatenate((signal_ORF, signal_3UTR))
    signal_stop /= tx_count

    return signal_start, signal_stop


def plot_metagene_signal(signal_start,
                         signal_stop,
                         w,
                         file):
    """
    Plot the mean signal upstream/downstream of the start and stop codon
    Args:
        signal_start: mean signal upstream/downstream of the first ORF nt
        signal_stop: mean signal upstream/downstream of the last ORF nt
        w: window of nucleotides upstream/downstream of the start or stop codon
            to return
        file: output plot file path
    Returns:

    Raises:

    """
    # Set the style and axes
    fig, (ax1,ax2) = plt.subplots(1, 2,
                                  sharex=False,
                                  sharey=True,
                                  figsize=(4.5,2.75))
    ax1.yaxis.set_major_formatter(fmt('%.2f'))

    # Plot the start signal
    signal_start = signal_start[-(w+15):]
    ax1.bar(np.linspace(-14, w, w+15),
            signal_start,
            color=basics,
            width=0.9)
    ax1.set_ylabel('Average rho', fontsize=12)
    ax1.set_xlabel('Position relative \nto start (nt)', fontsize=12)

    # Plot the stop signal
    signal_stop = signal_stop[:(w+15)]
    ax2.bar(np.linspace(-w, 14, w+15),
            signal_stop,
            color=basics,
            width = 0.9)
    ax2.set_xlabel('Position relative \nto stop (nt)', fontsize=12)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig(file+'.pdf', format='pdf')
    plt.savefig(file+'.eps', format='eps')
    plt.show()


def plot_emissions(data_df,
                   states,
                   file):
    """
    Plots a histogram of rho values for each type of canonical state
        (5'UTR, 3'UTR,
        start 1st nt, start 2nd nt, start 3rd nt,
        sense 1st nt, sense 2nd nt, sense 3rd nt,
        stop 1st nt, stop 2nd nt, stop 3rd nt)
    Args:
        data_df: processed data
        states: list of state names
        file: output plot file path
    Returns:

    Raises:

    """
    # Regroup the data by subcodon nucleotide position
    data_df['nt_pos'] = data_df['nt_pos'].astype(str)
    df = data_df.groupby(['codon_type', 'nt_pos']).agg(list).reset_index()
    df = df.sort_values('codon_type', ascending=False)

    # Set up the figure
    fig = plt.figure(figsize=(9,10))

    # Iterate through all rows of the emission distribution dataframe
    for i, (name, row) in enumerate(df.iterrows()):
        ax = plt.subplot(4, 3, i+1)

        # Get the reads data for this row (excluding values that are -inf)
        reads = np.array(row['norm_reads'])

        # Plot the reads data
        if i<9: color = basics[i%3]
        else: color = 'grey'
        ax.hist(reads, bins=np.linspace(xlim[0], xlim[1], 30),
                alpha=0.8, density=True, color=color)

        # Label the state
        ax.text(3.5, ylim[1]-ylim[1]*.2, row['codon_type'],
                fontsize=12, weight='bold', color=color)

        # Label the axes
        ax.set_xlabel('Rho', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig(file+'.pdf', format='pdf')
    plt.savefig(file+'.eps', format='eps')
    plt.show()


def get_prediction(path_h1,
                   path_h2,
                   altorfs_idx):
    """
    Get the positions of each altORF event from the viterbi path
    Args:
        path_h1: predicted state path
        path_h2: annotated state path
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
    Returns:
        annot_ORF: boolean array representing annotated ORF positions
        UTR5: boolean array representing 5'UTR positions
        UTR3: boolean array representing 3'UTR positions
        ORF: boolean array representing ORF positions
        uORF: boolean array representing uORF positions
        dORF: boolean array representing dORF positions
        pPRF: boolean array representing +1 PRF positions
        mPRF: boolean array representing -1 PRF positions
        SCR: boolean array representing SCR positions
    Raises:

    """
    if len(path_h2)==0:
        annot_ORF = np.zeros(len(path_h1))
    else:
        annot_ORF = np.array([s in altorfs_idx['ORF'] for s in path_h2])
    UTR5 = np.array([s==altorfs_idx['5UTR'] for s in path_h1])
    UTR3 = np.array([s==altorfs_idx['3UTR'] for s in path_h1])
    uORF = np.array([(s in altorfs_idx['uORF']) for s in path_h1])
    dORF = np.array([(s in altorfs_idx['dORF']) for s in path_h1])
    pPRF = np.array([(path_h1[i] in altorfs_idx['PRF']) and \
                     (path_h1[i+1] not in altorfs_idx['PRF']) \
            for i,s in enumerate(path_h1[:-1])] + [0])
    mPRF = np.array([(path_h1[i] in altorfs_idx['PRF']) and \
                     (path_h1[i+1] in altorfs_idx['PRF']) \
            for i,s in enumerate(path_h1[:-1])] + [0])
    SCR = np.array([(path_h1[i] in altorfs_idx['SCR']) and \
                    (path_h1[i+1] not in altorfs_idx['SCR']) \
            for i,s in enumerate(path_h1[:-1])] + [0])
    ORF = np.invert(np.logical_or(UTR5, UTR3))

    return annot_ORF, UTR5, UTR3, ORF, uORF, dORF, pPRF, mPRF, SCR


def plot_prediction(r,
                    path_h1,
                    path_h2,
                    altorfs_idx,
                    transcript_name,
                    file,
                    idx=0,
                    subcodon=True):
    """
    Plots the model prediction as shaded regions over the observed rho
        values for the given transcript
    Args:
        r: rho values
        path_h1: predicted state path
        path_h2: annotated state path
        altorfs_idx: dictionary mapping more general state type to
            state indices (0-indexed)
        transcript_name: name of the transcript
        file: output plot file path
        idx: true altORF position to plot for comparison, if there is one
        subcodon: whether to plot the signal broken out into subcodons 1, 2, 3
    Returns:

    Raises:

    """
    # Get the predicted positions
    annot_ORF, UTR5, UTR3, ORF, uORF, dORF, pPRF, mPRF, SCR = \
        get_prediction(path_h1,
                       path_h2,
                       altorfs_idx)
    start = sum(UTR5)-1

    if subcodon:

        # Split into subcodon signals
        nt1 = r.copy()
        nt1[start+1::3] = [0] * len(nt1[start+1::3])
        nt1[start+1::-3] = [0] * len(nt1[start+1::-3])
        nt1[start+2::3] = [0] * len(nt1[start+2::3])
        nt1[start+2::-3] = [0] * len(nt1[start+2::-3])

        nt2 = r.copy()
        nt2[start+0::3] = [0] * len(nt2[start+0::3])
        nt2[start+0::-3] = [0] * len(nt2[start+0::-3])
        nt2[start+2::3] = [0] *len(nt2[start+2::3])
        nt2[start+2::-3] = [0] * len(nt2[start+2::-3])

        nt3 = r.copy()
        nt3[start+0::3] = [0] * len(nt3[start+0::3])
        nt3[start+0::-3] = [0] * len(nt3[start+0::-3])
        nt3[start+1::3] = [0] * len(nt3[start+1::3])
        nt3[start+1::-3] = [0] * len(nt3[start+1::-3])

        # Create the plot
        fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,
                                    gridspec_kw={'height_ratios': [1.5,1,1,1]},
                                    sharex=True,
                                    sharey=True,
                                    figsize=(5, 4))
        for ax in [ax2,ax3,ax4]:
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            # Format the axes
            ax.set_xlim([0, len(r)])
            ax.set_ylim([0, max(r)+5])
            ax.yaxis.set_major_formatter(fmt('%.1f'))

        # Plot the signal
        ax1.plot(r, color='k', alpha=0.5)
        ax2.plot(nt1, color=colors[0], alpha=0.8)
        ax3.plot(nt2, color=colors[0], alpha=0.8)
        ax4.plot(nt3, color=colors[0], alpha=0.8)

        # Label the axes
        fig.add_subplot(111, frame_on=False)
        plt.tick_params(labelcolor='none', bottom=False, left=False)
        plt.xlabel('Position (nt)', fontsize=12)
        plt.ylabel('Rho', fontsize=12)

    else:

        # Create the plot
        fig,(ax1) = plt.subplots(1,1, figsize=(4, 2))

        # Plot the signal
        ax1.plot(r, color='k', alpha=0.5)

        # Format the axes
        ax1.set_xlim([0, len(r)])
        ax1.set_ylim([0, max(r)+5])
        ax1.yaxis.set_major_formatter(fmt('%.1f'))

        # Label the axes
        ax1.set_xlabel('Position (nt)', fontsize=12)
        ax1.set_ylabel('Rho', fontsize=12, rotation='vertical')

    # Title the figure
    plt.title(transcript_name, fontsize=14, weight='bold')

    # Save the unshaded plot
    plt.savefig(file+'.eps', format='eps')

    # Shade the annotated and predicted path on the first set of axes
    y = max(r)+5
    x = np.linspace(0, len(r), len(r), endpoint=False)
    ax1.fill_between(x, ORF*y,
                     facecolor=colors[1],
                     edgecolor=colors[1],
                     alpha=0.5,
                     linewidth=1)
    ax1.fill_between(x, uORF*y,
                     facecolor=colors[2],
                     edgecolor=colors[2],
                     alpha=0.5,
                     linewidth=1)
    ax1.fill_between(x, dORF*y,
                     facecolor=colors[2],
                     edgecolor=colors[2],
                     alpha=0.5,
                     linewidth=1)
    ax1.fill_between(x, pPRF*y,
                     facecolor = 'none',
                     edgecolor=colors[5],
                     linewidth=1)
    ax1.fill_between(x, mPRF*y,
                     facecolor = 'none',
                     edgecolor=colors[6],
                     linewidth=1)
    ax1.fill_between(x, SCR*y,
                     facecolor = 'none',
                     edgecolor=colors[7],
                     linewidth=1)
    ax1.fill_between(x, annot_ORF*y,
                     facecolor = 'none',
                     edgecolor=colors[1],
                     linewidth=1,
                     linestyle=':')

    # Plot the true location, if there is one
    if idx!=0:
        idx = start + idx
        ax1.axvline(idx, 0, y, color='k', linestyle=':')

    # Save and show the plot
    plt.savefig(file+'.pdf', format='pdf')
    plt.show()


def orfeus(data_df,
           transcript_id,
           transcript_name,
           parameters_h0,
           parameters_h1,
           idx=0,
           suppress=True,
           debug=False,
           subcodon=True):
    """
    Displays the predicted path and plots
    Args:
        data_df: transcript data
        transcript_id: ID of the transcript
        transcript_name: name of the transcript
        parameters_h1: list of all parameter matrices for the current model
        parameters_h0: list of all parameter matrices for the null model
        idx: true altORF position to plot for comparison, if there is one
        suppress: whether to suppress the path output
        debug: whether to print the full traceback
        subcodon: whether to plot the signal broken out into subcodons 1, 2, 3
    Returns:

    Raises:

    """
    # Unpack the parameters
    states, states_idx, altorfs_idx, \
    p_emit_nt_h1, p_emit_rho_h1, rho_ranges, p_trans_h1, p_trans_rev_h1, \
    prev_states_options_h1, next_states_options_h1 = parameters_h1

    # Extract just the transcript of interest
    df = data_df[data_df['transcript_id']==transcript_id]

    # Get the transcript coverage, sequence, and path
    c = df['mean_reads'].tolist()[0]
    x = df['seq'].tolist()
    r = df['log_reads'].tolist()
    path = get_annotated_path(df['state'].tolist(), states)

    print('Mean reads per nucleotide: %f' % c)

    path_h2, path_h1, path_h0, \
    path_score_h2, path_score_h1, path_score_h0, log_odds_score, \
    log_emit_h1, log_trans_h1 = hmm_single(x, r,
                                           path,
                                           parameters_h1,
                                           parameters_h0)

    # Print the annotated path, predicted path, and null path
    if not suppress:

        #seq_h2 = readable_sequence(x, path_h2, states_idx)
        #print('Annotated:               ', ''.join(seq_h2))

        #seq_h0 = readable_sequence(x, path_h0, states_idx)
        #print('Null:                    ', ''.join(seq_h0))

        seq_h1 = readable_sequence(x, path_h1, states_idx)
        print('Predicted:               ', ''.join(seq_h1))

        #print('Annotated path score:    %f' % path_score_h2)
        #print('Null path score:         %f' % path_score_h0)
        #print('Predicted path score:    %f' % path_score_h1)

    # Print the full trace for debugging
    if debug:
        trace_dump(x, r, path_h1, states, log_emit_h1, log_trans_h1)

    # Plot the prediction
    r = df['norm_reads'].tolist()
    file = transcript_name
    plot_prediction(r,
                    path_h1,
                    path_h2,
                    altorfs_idx,
                    transcript_name,
                    file,
                    idx,
                    subcodon)


def plot_results(results, label, file):
    """
    Plots mean riboseq coverage per transcript vs mean sensitivity or specificity
    Args:
        results: dict mapping mean riboseq coverage to mean sensitivity or specificity
        label: title of plot (sensitivity or specificity)
        file: output file for plot
    Returns:

    Raises:

    """
    fig, ax = plt.subplots(figsize=(4,3))

    for event in results:

        if event=='uORFdORF': c = colors[2]
        elif event=='pPRF': c = colors[5]
        elif event=='mPRF': c = colors[6]
        elif event=='SCR': c = colors[7]
        elif event=='ORF': c = colors[1]
        elif event=='any': c = 'k'

        # Plot the results
        x = list(results[event].keys())
        y = list(results[event].values())
        ax.scatter(x, y, color=c, alpha=0.7, marker='o')

    # Label the axes
    plt.xlim([0,1])
    plt.ylim([-0.1,1.1])
    plt.xlabel('Mean ORF coverage (footprints/nt)', fontsize=12)
    plt.ylabel(label, fontsize=12, rotation='vertical')
    plt.title(label, fontsize=14, weight='bold')

    # Save the plot
    plt.tight_layout()
    plt.savefig(file+'.pdf', format='pdf')
    plt.savefig(file+'.eps', format='eps')
    plt.show()


def plot_transcript_coverage(results, file):
    """
    Plots number of transcripts with at least the given mean riboseq coverage
    Args:
        results: dict mapping mean riboseq coverage to number of transcripts
        label: title of plot (sensitivity or specificity)
        file: output file for plot
    Returns:

    Raises:

    """
    fig, ax = plt.subplots(figsize=(4.2,3))

    # Plot the results
    x = list(results.keys())
    y = list(results.values())
    ax.scatter(x, y, color=colors[0], alpha=0.8)

    # Label the axes
    plt.xlim([0,1])
    plt.ylim([0,7000])
    plt.xlabel('Mean ORF coverage (footprints/nt)', fontsize=12)
    plt.ylabel('Number of transcripts', fontsize=12, rotation='vertical')
    plt.title('Number of transcripts', fontsize=14, weight='bold')

    # Save the plot
    plt.tight_layout()
    plt.savefig(file+'.pdf', format='pdf')
    plt.savefig(file+'.eps', format='eps')
    plt.show()
