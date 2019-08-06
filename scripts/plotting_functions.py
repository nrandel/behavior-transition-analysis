from classes.cell_trace_config import cell_trace_config_filter, get_event_start_rows

import logging
import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd


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


def plot_windowed_events(all_events):
    all_df = all_events.pop(0)
    for right_df in all_events:
        print(right_df.columns)
        all_df = pd.merge_ordered(all_df, right_df, on="time", how="outer")

    # Resets the index as time and drops time column (sollte spaeter kommen)
    all_df.index = all_df["time"]
    del all_df["time"]
    # print(all_df)

    # Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
    int_all_df = all_df.interpolate(
        method="index", axis=0, limit=None, inplace=False, limit_direction="both"
    )
    print(int_all_df.columns)

    # For single or multiple sample, aligned for certain event (no transitions taken into account).
    # Interpolated data used!
    # Averaged all events and all cells pro timepoint
    # Average for specific cell type filter-pattern (see below)
    all_cell_avg_df = int_all_df.mean(axis=1)
    all_cell_min_df = int_all_df.min(axis=1)
    all_cell_max_df = int_all_df.max(axis=1)
    # Standard deviation (distribution)
    all_cell_std_df = int_all_df.std(axis=1)
    # standard error of mean
    all_cell_sem_df = int_all_df.sem(axis=1)
    # print(all_cell_avg_df)

    fig = plt.figure()

    # Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
    sub1 = fig.add_subplot(211)
    all_df.plot(ax=sub1, marker="*", label="A00C")
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
