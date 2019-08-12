from classes.cell_trace_config import cell_trace_config_filter, get_event_start_rows

import logging
import matplotlib.pyplot as plt
from matplotlib import ticker

from functions import merge_dataframe_list


# Plotting - single sample for whole time
def layout_plot(plot, tick_spacing=100, fov=(0, 2400, 0, 1.2), legend=False):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


def aligned_layout_plot(
    plot, tick_spacing=0.5, fov=(-20, 50, -0.05, 1.9), legend=False
):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


# plots single event
def plot_single_event(filtered_data, sample_id, cell_trace_config, event_df):
    # TODO: Handle plotting multiple celltypes/filterpatterns with legend

    cell_subset_df = filtered_data.get(sample_id)
    # Single sample analysis
    # For single sample over the whole experimental time
    # Calculate min, max, avg, stddev, sem from cell_subset_df (defined earlier)
    cell_subset_df.index = cell_subset_df.time
    del cell_subset_df["time"]  # delete time_column

    if len(cell_subset_df.columns) == 0:
        logging.warning("Sample {} has no data to plot".format(sample_id))
        return None
    cell_subset_df.index.name = None  # delete index name
    event_avg_df = cell_subset_df.mean(axis=1)
    event_min_df = cell_subset_df.min(axis=1)
    event_max_df = cell_subset_df.max(axis=1)
    # Standard deviation (distribution)
    event_std_df = cell_subset_df.std(axis=1)
    # standard error of mean
    event_sem_df = cell_subset_df.sem(axis=1)
    # print(ctc.sample_id, cell_avg_df)

    fig = plt.figure()
    fig.suptitle(sample_id)
    fig.set_facecolor("white")

    # Plot all cells from cell_subset_df over entire time (specified in Cell_Trace_Config). <For a single sample>
    sub1 = fig.add_subplot(211)  # 211
    cell_subset_df.plot(ax=sub1)
    cell_trace_config.add_event_time_points_to_plot(event_df, sub1)
    layout_plot(sub1)

    # Avg, min, max, std-dev for multiple cells in single sample over whole time.
    # All cells are averaged. For cell and filterpattern see below.
    sub2 = fig.add_subplot(212)  # 212
    event_avg_df.plot(
        ax=sub2, color="g", label=cell_trace_config.cell_type, linewidth=1
    )
    cell_trace_config.add_event_time_points_to_plot(event_df, sub2)
    # cell_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
    # cell_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
    event_avg_df.plot.line(yerr=event_std_df, ax=sub2, color="r", alpha=0.1)
    # cell_avg_df.plot.line(yerr=cell_sem_df, ax=sub2, color = 'c', alpha = 0.1)
    layout_plot(sub2)
    plt.show()


def plot_all_events(sample_data, cell_trace_configs):
    cell_trace_config_id_map = {ctc.sample_id: ctc for ctc in cell_trace_configs}
    filtered_data = cell_trace_config_filter(sample_data, cell_trace_configs)
    event_start_dfs = get_event_start_rows(sample_data, cell_trace_configs)
    for sample_id in filtered_data:
        plot_single_event(
            filtered_data,
            sample_id,
            cell_trace_config_id_map[sample_id],
            event_start_dfs[sample_id],
        )


def plot_windowed_events(df):
    if len(df) == 0:
        logging.warning("No data to plot!")
        return
    # For single or multiple sample, aligned for certain event (no transitions taken into account).
    # Interpolated data used!
    # Averaged all events and all cells pro timepoint
    # Average for specific cell type filter-pattern (see below)

    all_cell_avg_df = df.mean(axis=1)
    all_cell_min_df = df.min(axis=1)
    all_cell_max_df = df.max(axis=1)
    all_cell_std_df = df.std(axis=1)
    all_cell_sem_df = df.sem(axis=1)
    logging.info("Ignores filterpattern in statistic calculations")

    fig = plt.figure()

    # Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
    sub1 = fig.add_subplot(211)
    df.plot(ax=sub1, marker="*", label="A00C")
    aligned_layout_plot(sub1)

    sub2 = fig.add_subplot(212)
    all_cell_avg_df.plot(
        ax=sub2, color="k", label="A00C"
    )  # use interpolated df to calculate average...
    # all_cell_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
    # all_cell_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
    all_cell_avg_df.plot.line(yerr=all_cell_std_df, ax=sub2, color="r", alpha=0.1)
    # all_cell_avg_df.plot.line(yerr=all_cell_sem_df, ax=sub2, color = 'c', alpha = 0.1)

    aligned_layout_plot(sub2)
    plt.show()


from logging_functions import log_overlapping_transition
import numpy as np


def plot_transition_gaps_hist(found_transitions):
    # get event diffs
    gap_Ptrans = []
    for sample in found_transitions:
        for found_transition in sample:
            gap_Ptrans.append(
                (found_transition["second_event_start"])
                - (found_transition["first_event_end"])
            )
            log_overlapping_transition(found_transition, seconds_overlap=0)

    # print(gap_Ptrans)

    avg_duration = np.mean(gap_Ptrans)
    max_duration = np.max(gap_Ptrans)
    min_duration = np.min(gap_Ptrans)

    logging.info("avg gap duration: {}".format(avg_duration))
    logging.info("max gap duration: {}".format(max_duration))
    logging.info("min gap duration: {}".format(min_duration))

    # Histogram
    fig = plt.figure()
    plt.hist(gap_Ptrans, bins=100, alpha=0.5)
    plt.show()


def plot_transitions(df, transition_types, use_sem=False):
    colors = ["cyan", "red", "violet", "red", "blue"]

    fig = plt.figure()
    sub2 = fig.add_subplot(111)
    for tt, color in zip(transition_types, colors):
        tt_df = df.filter(
            regex=tt.get_filter_regex(
                use_cell=True,
                use_filter_pattern=True,
                use_first_event=True,
                use_second_event=True,
            )
        )

        tt_avg_df = tt_df.mean(axis=1)
        tt_min_df = tt_df.min(axis=1)
        tt_max_df = tt_df.max(axis=1)
        tt_std_df = tt_df.std(axis=1)
        tt_sem_df = tt_df.sem(axis=1)

        tt_avg_df.plot(ax=sub2, label=tt.get_filter_regex(use_all=True), color=color)
        if use_sem:
            tt_avg_df.plot(yerr=tt_sem_df, ax=sub2, color="grey", alpha=0.5)
        else:
            tt_avg_df.plot(
                yerr=tt_std_df,
                ax=sub2,
                # label=tt.get_filter_regex(use_all=True),
                alpha=0.5,
                color="grey",
            )
    aligned_layout_plot(sub2, tick_spacing=5, fov=(-18.5, 42.4, 0.0, 1.0), legend=True)
    plt.show()


def plot_fold_change(df, transitions):

    # Plot fold change
    ax = df.plot.box()  # its a series (diff for data frame)
    ax.set_title("")
    ax.set_xlabel("Transitions")
    ax.set_ylabel("Fold change")
    ax.set_xticklabels(transitions)

    # ax = fold_change_df.plot.box() #single transition type
    plt.show()
