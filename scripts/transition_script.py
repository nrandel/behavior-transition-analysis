# %%
# Imports and Load data

import logging
from functions import get_sample_data
from classes import CellTraceConfig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import tqdm

from logging_functions import count_cell_filterpatterns
from plotting_functions import plot_all_events
from classes.cell_trace_config import extract_windows, cell_trace_config_filter
from plotting_functions import plot_windowed_events

from classes import PostBehaviorTransition

from transition_functions import find_behavior_before

PLOT_EVERY_SAMPLE = True

# Logging info:
# Logging has 4 levels (DEBUG, INFO, WARNING, ERROR)
# appends to log file rather than overwritting it
logging.basicConfig(level=logging.INFO, filename="log.log")

# TODO: Cleanup comments in get_sample_data
sample_data = get_sample_data()

# %%
#########################################################################################################################
# transitions
#########################################################################################################################
# %%
behavior_transitions = [
    PostBehaviorTransition("17-08-26L2-cl", "turn", "fw", 4.9),
    PostBehaviorTransition("17-08-26L6-cl", "turn", "fw", 4.9),
]

# Open all samples single Transitions
# behavior_transitions = [
#    PostBehaviorTransition(name, 'stim', 'fw', 4.9) for name in sample_data]

"""
# Open all samples multiple Transitions
behavior_transitions = [
    PostBehaviorTransition(name, 'fw', 'stim', 4.9) for name in sample_data] + [
    PostBehaviorTransition(name, 'bw', 'stim', 4.9) for name in sample_data] + [
    PostBehaviorTransition(name, 'turn', 'stim', 4.9) for name in sample_data] + [
    PostBehaviorTransition(name, 'hunch', 'stim', 4.9) for name in sample_data] + [
    PostBehaviorTransition(name, 'other', 'stim', 4.9) for name in sample_data]
"""

from transition_functions import extract_transitions

found_transitions = extract_transitions(sample_data, behavior_transitions)

from logging_functions import log_num_transitions

log_num_transitions(found_transitions)

# print(found_transitions) # Transitions
# %%
# Duration between diff behavior from PostTransition (min, max, avg)
# Note: only works for one type of transitions
from plotting_functions import plot_transition_gaps_hist

plot_transition_gaps_hist(found_transitions)


raise Exception("STOP HERE!")

# %%
# Define celltype, filter-pattern for transitions. Optional do not run regex-filter. Aligning to first event start or
# second event start: depends on stimulus. Output: a) aligned raw values, b) aloigned interpolated values.
# tqdm = progress bar

from classes import CellTransConfig

cell_Ptrans_configs = []
all_Ptrans_events = []

# Turn list of list of found transition dictionaries, into list of CellTransConfigs
for sample in tqdm(found_transitions):
    for found_transition in sample:

        # For all behavior except stimulus as first event
        cell_Ptrans_configs.append(
            CellTransConfig(
                found_transition["sample_id"],
                "basin",
                found_transition["second_event_start"],
                filter_pattern="A1R",
                first_event=found_transition["first_event"],
                second_event=found_transition["second_event"],
            )
        )

        # For stimulus as first_event
        # cell_Ptrans_configs.append(CellTransConfig(found_transition["sample_id"], "basin",
        #                                           found_transition["first_event_start"], 'A6',
        #                                           first_event=found_transition["first_event"],
        #                                           second_event=found_transition["second_event"]))

# Extract for specific time window and align several events.
# Define timepoints pre and post an event (event_df).
# This works for single sample or multiple samples aligned
# Note: In cell_subset_df, time was set to index, because for the previous avg calculation
# Add index and time = column

# Set the window range left and right from the event (in seconds)
left_half_window_size = (
    18.5
)  # If negative it goes further to right half (Good for skipping stimulus)
right_half_window_size = 42.4

# trans_df defined in pargraph before
windows = []
n_behavior_per_sample = {}

# TODO: Split filter data and extract windows
for ctc in tqdm(cell_Ptrans_configs):
    sample_df = sample_data.get(ctc.sample_id)
    n_behavior = n_behavior_per_sample.get(ctc.sample_id, 1)
    n_behavior_per_sample.setdefault(ctc.sample_id, 1)
    if sample_df is None:
        raise ValueError("{}: could not find sample data".format(ctc.sample_id))
        continue

    # Extract columns matching our cell type and the optional filter pattern.
    # cell_subset_df = sample_df.filter(regex=ctc.get_filter_regex()) #Get subset of cells
    # cell_subset_df.set_index(sample_df.time, inplace=True) #Set time to index (essential for min/max...)
    # cell_subset_df.reset_index(inplace = True) # Add index and time = column

    # Don't apply filter regex, but take all cells from lm_data
    # Try and except for cases when time was added to lm_data before (by partially running the notebook)

    # USE LM DATA SINCE IT DOESNT HAVE BEHAVIOR COLUMNS (makes filtering easier)
    cell_subset_df = lm_data.get(ctc.sample_id)  # Get subset of cells
    try:
        cell_subset_df.set_index(
            sample_df.time, inplace=True
        )  # Set time to index (essential for min/max...)
        cell_subset_df.reset_index(
            inplace=True, drop=True
        )  # Add index and time = column
    except:
        pass

    n_behavior_per_sample[ctc.sample_id] += 1
    window_start = ctc.event_time - left_half_window_size
    window_end = ctc.event_time + right_half_window_size

    # Get subset of rows between window_start and window_end
    # Including event_start
    # trans = cell_subset_df[(cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)]
    # Excluding event start
    trans = cell_subset_df[
        (cell_subset_df.time > window_start) & (cell_subset_df.time < window_end)
    ]
    # Normalizing the data to align on beginning of selected
    # behavior (event_df = Zero) by substracting events in window
    # around start of event of interest from start of event interest.
    # Note: using ":" in event.loc[] will select "all rows" in our window.
    # trans.loc[:, 'time'] = trans['time'] - row['time']
    trans.loc[:, "time"] = trans["time"] - ctc.event_time

    # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
    trans.rename(
        lambda x: "{}_{}_{}_{}_{}".format(
            ctc.sample_id, x, n_behavior, ctc.first_event, ctc.second_event
        ),
        axis="columns",
        inplace=True,
    )

    # Rename time collum to time
    trans.rename(columns={trans.columns[0]: "time"}, inplace=True)
    all_Ptrans_events.append(trans)  # Append a list with all event
# print(all_Ptrans_events)

# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
all_Ptrans_df = all_Ptrans_events.pop(0)
for right_df in all_Ptrans_events:
    all_Ptrans_df = pd.merge_ordered(all_Ptrans_df, right_df, on="time", how="outer")
# print(all_Ptrans_df)

# Resets the index as time and drops time column
all_Ptrans_df.index = all_Ptrans_df["time"]
del all_Ptrans_df["time"]
# print(all_Ptrans_df)

# Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
int_all_Ptrans_df = all_Ptrans_df.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)
print(int_all_Ptrans_df.columns)

# %%
# Figure out how many cells of a determine type and filter pattern are involved in a specific transition
def extract_parts(column: str, indicies):
    parts = column.split("_")
    return tuple(parts[i] for i in indicies)


# indicies: 0: sample_id, 1: cell, 2: filter_pattern, 3: n_obs, 4: first_event, 5:second_event
parts = [extract_parts(column, [0, 1, 2]) for column in int_all_Ptrans_df.columns]
print(len(parts))
print(len(set(parts)))
# print(set(parts))
# %%
# For multiple transition events, group after transition (first, or second event) <most useful> with option to group
# after celltype, filterpattern, sample_id, observations, using class TransitionType

transition_types = [
    # TransitionType(first_event = 'fw'),
    # TransitionType(first_event = 'bw'),
    # TransitionType(first_event = 'turn'),
    # TransitionType(first_event = 'hunch'),
    # TransitionType(first_event = 'other'),
    TransitionType(first_event="stim")
]

print(transition_types)
# %%

# For multiple transition types!!
# If a dataframe with NANs is plotted (raw-data = non interpolated), use
# marker = '+', or 'o', since the line in the lineplot only connects
# consecutive data points
def aligned_layout_plot(
    plot, tick_spacing=5, fov=(-18.5, 42.4, 0.0, 1.0), legend=False
):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


colors = ["cyan", "red", "violet", "red", "blue"]
# colors = {tt.get_filter_regex(use_cell=True, use_filter_pattern=True, use_first_event=True, use_second_event=True) : color for tt, color in zip(transition_types, colors)}

print(colors)

fig = plt.figure()
sub2 = fig.add_subplot(111)
for tt, color in zip(transition_types, colors):
    int_some_Ptrans_df = int_all_Ptrans_df.filter(
        regex=tt.get_filter_regex(
            use_cell=True,
            use_filter_pattern=True,
            use_first_event=True,
            use_second_event=True,
        )
    )

    # Average and stddev, min, max, sem for post_behavior_transition events
    all_Ptrans_avg_df = int_some_Ptrans_df.mean(axis=1)  # Interpolated data used
    all_Ptrans_min_df = int_some_Ptrans_df.min(axis=1)
    all_Ptrans_max_df = int_some_Ptrans_df.max(axis=1)
    # Standard deviation (distribution)
    all_Ptrans_std_df = int_some_Ptrans_df.std(axis=1)
    # standard error of mean
    all_Ptrans_sem_df = int_some_Ptrans_df.sem(axis=1)

    all_Ptrans_avg_df.plot(
        ax=sub2, label=tt.get_filter_regex(use_all=True), color=color
    )  # use interpolated df to calculate average...
    # all_Ptrans_avg_df.plot(yerr=all_Ptrans_std_df, ax=sub2, label = tt.get_filter_regex(use_all=True), alpha = 0.005, color = color)
    all_Ptrans_avg_df.plot.line(
        yerr=all_Ptrans_sem_df, ax=sub2, color="grey", alpha=0.5
    )
aligned_layout_plot(sub2, legend=True)

# %%
# To make the plot in the notebook and not in an extra window

# Only for one type of transitions (I still need this!!)

# Average and stddev, min, max, sem for post_behavior_transition events
all_Ptrans_avg_df = int_all_Ptrans_df.mean(axis=1)  # Interpolated data used
all_Ptrans_min_df = int_all_Ptrans_df.min(axis=1)
all_Ptrans_max_df = int_all_Ptrans_df.max(axis=1)
# Standard deviation (distribution)
all_Ptrans_std_df = int_all_Ptrans_df.std(axis=1)
# standard error of mean
all_Ptrans_sem_df = int_all_Ptrans_df.sem(axis=1)
# wrong zur haelfte: Want to have avg per celltyp over time point,
# and not avg over all cells per timepoint (refer to Data_filter or Grouper)

# Plotting for multi-events (same_behavioral_transition)
# If a dataframe with NANs is plotted (raw-data = non interpolated), use
# marker = '+', or 'o', since the line in the lineplot only connects
# consecutive data points
def aligned_layout_plot(
    plot, tick_spacing=5, fov=(-18.5, 42.4, 0.0, 1.0), legend=False
):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
sub1 = fig.add_subplot(111)  # 211
all_Ptrans_df.plot(ax=sub1, marker="*", label=ctc.cell_type)
aligned_layout_plot(sub1)

# sub2 = fig.add_subplot(212) #212
# all_Ptrans_avg_df.plot(ax=sub2, color = 'c', label = ctc.cell_type) #use interpolated df to calculate average...
# all_Ptrans_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
# all_Ptrans_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
# all_Ptrans_avg_df.plot.line(yerr=all_Ptrans_std_df, ax=sub2, color = 'lightgrey', alpha = 0.1)
# all_Ptrans_avg_df.plot.line(yerr=all_Ptrans_sem_df, ax=sub2, color = 'grey', alpha = 0.1)
# aligned_layout_plot(sub2)
# %%
# all all_Ptrans_df with left and right window in same data frame and aligned to second_event_start
# (for stim first event on first_event_start_)
# Average over the frames within the same sample (eg average over 1 sec)
# extract all negative values for pre_event
# extract positiv values for post event, but in case of the stimulus, the first two seconds must not be included


pre_data = all_Ptrans_df[all_Ptrans_df.index < 0.0]
post_data = all_Ptrans_df[
    all_Ptrans_df.index > 2.0
]  # for stim-behavior >2, otherwise >=0

int_pre_data = pre_data.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)
int_post_data = post_data.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)

# Average over time for each cell type
pre_data_avg = pre_data.mean(axis=0)
post_data_avg = post_data.mean(axis=0)
# %%
# Fold change (bar plot) in cell activity post - pre/pre (merged orderd,No interpolation, NO avg)
# Fold change for all cells individually

# fold changeI: post/pre
# fold_change = (post_data_avg)/(pre_data_avg)

# fold changeII: post-pre/pre
fold_change = ((post_data_avg) - (pre_data_avg)) / (pre_data_avg)
# print(fold_change.index)

# Transform Series to dataframe and name columns
fold_change_df = fold_change.to_frame("transitions")
# print(fold_change_df.index)
# transitions = set(["_".join(sample_transition.split("_")[-2:]) for sample_transition in fold_change_df.index])

transitions = []
past = None
for sample_transition in fold_change_df.index:
    current = "_".join(sample_transition.split("_")[-2:])
    if current != past:
        transitions.append(current)
        past = current
print(transitions)

print(transitions)
transition_df_map = {
    transition: fold_change_df[fold_change_df.index.str.contains(transition)]
    for transition in transitions
}
print(list(transition_df_map.keys()))
# print(transition_df_map["fw_stim"].head(10))
# %%
# Merge
# df_key = name, df-list = dataframe
transition_df_list = list(transition_df_map.items())
df_key, all_fold_change_df = transition_df_list.pop(0)
for df_key, right_df in transition_df_list:
    all_fold_change_df = pd.merge(
        all_fold_change_df, right_df, left_index=True, right_index=True, how="outer"
    )

# print(all_fold_change_df.to_string()) #print everything

# Plot fold change
ax = all_fold_change_df.plot.box()  # its a series (diff for data frame)
ax.set_title("")
ax.set_xlabel("Transitions")
ax.set_ylabel("Fold change")
ax.set_xticklabels(list(transition_df_map.keys()))

# ax = fold_change_df.plot.box() #single transition type
plt.show()
# %%
# Mann-Witney U test (non-parametric pair wise test)

# da ist ne Menge significant:()

from scipy.stats import stats
from scipy.stats import mannwhitneyu

# mwu1 = mannwhitneyu(transition_df_map["fw_stim"], transition_df_map["bw_stim"])
# mwu2 = mannwhitneyu(transition_df_map["fw_stim"], transition_df_map["turn_stim"])
# mwu3 = mannwhitneyu(transition_df_map["fw_stim"], transition_df_map["hunch_stim"])
# mwu4 = mannwhitneyu(transition_df_map["fw_stim"], transition_df_map["other_stim"])
# mwu5 = mannwhitneyu(transition_df_map["bw_stim"], transition_df_map["turn_stim"])
# mwu6 = mannwhitneyu(transition_df_map["bw_stim"], transition_df_map["bw_stim"])

# print(mwu1, mwu2, mwu3, mwu4, mwu5, mwu6)
# print(mwu1)
# %%
# Triple transition
# Find first and  second transition between two different! behaviors with intersection:
# first event end = second event start, within a max_delay.

# Open single samples
"""
first_transitions = [
    PostBehaviorTransition('17-08-26L6-cl', 'fw', 'stim', 3)
]
second_transitions = [
    PostBehaviorTransition('17-08-26L6-cl', 'stim', 'fw', 3)
]
"""
# Open all samples

first_transitions = [PostBehaviorTransition(name, "bw", "stim", 3) for name in lm_data]

second_transitions = [
    PostBehaviorTransition(name, "stim", "quiet", 3) for name in lm_data
]


found_transitions = []
for first_bt, second_bt in tqdm(zip(first_transitions, second_transitions)):
    transitions = []
    assert first_bt.sample_id == second_bt.sample_id, "{} does not match {}".format(
        first_bt.sample_id, second_bt.sample_id
    )
    sample_df = sample_data.get(first_bt.sample_id)
    if sample_df is None:
        raise ValueError("No data found for sample {}".format(bt.sample_id))
    if not any(["bw" in column for column in sample_df.columns]):
        continue
    first_transition_duration = None
    second_transition_duration = None
    third_transition_duration = 2
    first_transitions = find_behavior_before(
        first_bt.sample_id,
        sample_df,
        first_bt.event,
        first_bt.post_event,
        first_bt.max_delay,
        first_event_duration=first_transition_duration,
        second_event_duration=second_transition_duration,
    )
    second_transitions = find_behavior_before(
        second_bt.sample_id,
        sample_df,
        second_bt.event,
        second_bt.post_event,
        second_bt.max_delay,
        first_event_duration=second_transition_duration,
        second_event_duration=third_transition_duration,
    )
    # print("{} transitions from {} to {}".format(len(first_transitions), first_bt.event, first_bt.post_event))
    # print("{} transitions from {} to {}".format(len(second_transitions), second_bt.event, second_bt.post_event))

    for ft in first_transitions:
        for st in second_transitions:
            if abs(ft["second_event_start"] - st["first_event_start"]) < 0.00001:
                transitions.append(
                    {
                        "sample_id": ft["sample_id"],
                        "first_event_start": ft["first_event_start"],
                        "first_event_end": ft["first_event_end"],
                        "second_event_start": st["first_event_start"],
                        "second_event_end": st["first_event_end"],
                        "third_event_start": st["second_event_start"],
                    }
                )
    if transitions:
        print("{} transition triples found".format(len(transitions)))
        found_transitions.append(transitions)


print(len(found_transitions))  # number of data sets not the actual stim
print(sum([len(sample_transitions) for sample_transitions in found_transitions]))
print(found_transitions)
# %%
# Use the predefined CellTransConfig to filter by celltype and pattern.
# The results are merged_ordered and interpolated.
cell_Ttrans_configs = []
all_Ttrans_events = []

for sample in tqdm(found_transitions):
    sample_ls_trans = []
    for found_transition in sample:

        # For all behavior
        sample_ls_trans.append(found_transition["second_event_start"])
        cell_Ttrans_configs.append(
            CellTransConfig(
                found_transition["sample_id"],
                "A00c",
                found_transition["second_event_start"],
            )
        )

# Extract for specific time window and align several events.
# Define timepoints pre and post an event (event_df).
# This works for single sample or multiple samples aligned
# Note: In cell_subset_df, time was set to index, because for the previous avg calculation
# Add index and time = column

# Set the window range left and right from the event (in seconds)
left_half_window_size = 100.0
right_half_window_size = 200.0

# trans_df defined in pargraph before
windows = []
n_behavior = 0
# TODO: Split filter columns and  extract windows
for ctc in tqdm(cell_Ttrans_configs):
    sample_df = sample_data.get(ctc.sample_id)
    if sample_df is None:
        raise ValueError("{}: could not find sample data".format(ctc.sample_id))
        continue
    # Extract columns matching our cell type and the optional filter pattern.
    # Pandas' filter() operations works on columns for DataFrames by default.
    cell_subset_df = sample_df.filter(
        regex=ctc.get_filter_regex()
    )  # Get subset of cells
    cell_subset_df.set_index(
        sample_df.time, inplace=True
    )  # Set time to index (essential for min/max...)
    cell_subset_df.reset_index(inplace=True)  # Add index and time = column
    # print(cell_subset_df)

    # Don't apply filter regex, but take all cells from lm_data
    # cell_subset_df =  lm_data.get(ctc.sample_id)#Get subset of cells
    # cell_subset_df.reset_index(inplace = True, drop = True) # Add index and time = column
    # print(cell_subset_df)

    n_behavior += 1
    window_start = ctc.event_time - left_half_window_size
    window_end = ctc.event_time + right_half_window_size

    # Get subset of rows between window_start and window_end
    trans = cell_subset_df[
        (cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)
    ]
    # Normalizing the data to align on beginning of selected
    # behavior (event_df = Zero) by substracting events in window
    # around start of event of interest from start of event interest.
    # Note: using ":" in event.loc[] will select "all rows" in our window.
    # trans.loc[:, 'time'] = trans['time'] - row['time']
    trans.loc[:, "time"] = trans["time"] - ctc.event_time

    # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
    trans.rename(
        lambda x: "{}_{}_{}".format(ctc.sample_id, x, n_behavior),
        axis="columns",
        inplace=True,
    )

    # Rename time collum to time
    trans.rename(columns={trans.columns[0]: "time"}, inplace=True)
    all_Ttrans_events.append(trans)  # Append a list with all event

# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
all_Ttrans_df = all_Ttrans_events.pop(0)
for right_df in all_Ttrans_events:
    all_Ttrans_df = pd.merge_ordered(all_Ttrans_df, right_df, on="time", how="outer")

# Resets the index as time and drops time column
all_Ttrans_df.index = all_Ttrans_df["time"]
del all_Ttrans_df["time"]
# print(all_Ttrans_df)

# Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
int_all_Ttrans_df = all_Ttrans_df.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)
# print(int_all_Ttrans_df)
# %%
# Average and stddev, min, max, sem for post_behavior_transition events
all_Ttrans_avg_df = int_all_Ttrans_df.mean(axis=1)  # Interpolated data used
all_Ttrans_min_df = int_all_Ttrans_df.min(axis=1)
all_Ttrans_max_df = int_all_Ttrans_df.max(axis=1)
# Standard deviation (distribution)
all_Ttrans_std_df = int_all_Ttrans_df.std(axis=1)
# standard error of mean
all_Ttrans_sem_df = int_all_Ttrans_df.sem(axis=1)
# wrong zur haelfte: Want to have avg per celltyp over time point,
# and not avg over all cells per timepoint

# Plotting for multi-events (same_behavioral_transition)
# If a dataframe with NANs is plotted (raw-data = non interpolated), use
# marker = '+', or 'o', since the line in the lineplot only connects
# consecutive data points
def aligned_layout_plot(plot, tick_spacing=10, fov=(-20, 50, 0.0, 1.0), legend=False):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
# sub1 = fig.add_subplot(211)
# all_Ttrans_df.plot(ax=sub1, marker = '*', label = ctc.cell_type)
# aligned_layout_plot(sub1)

sub2 = fig.add_subplot(111)  # 212
all_Ttrans_avg_df.plot(
    ax=sub2, color="r", label=ctc.cell_type
)  # use interpolated df to calculate average...
# all_Ttrans_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
# all_Ttrans_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
all_Ttrans_avg_df.plot.line(
    yerr=all_Ttrans_std_df, ax=sub2, color="lightgrey", alpha=0.1
)
# all_Ttrans_avg_df.plot.line(yerr=all_Ttrans_sem_df, ax=sub2, color = 'grey', alpha = 0.1)
aligned_layout_plot(sub2)
# %%
# Open single sample

# behavior_transitions = [
#    SamePairBehaviorTransition('17-08-26L1-cl', 'turn', 'turn', max_delay=10, max_ignored_quiet_time = math.inf)
# SamePairBehaviorTransition('17-08-24L1-cl', 'fw', 'fw', max_delay=10, max_ignored_quiet_time = math.inf)
# SamePairBehaviorTransition('17-11-06L3-cl', 'bw', 'bw', 3),
# SamePairBehaviorTransition('17-11-29L1-cl', 'bw', 'bw', 3)
# ]


# Open all samples
behavior_transitions = [
    SamePairBehaviorTransition(
        name, "turn", "turn", max_delay=10, max_ignored_quiet_time=math.inf
    )
    for name in lm_data
]


found_transitions = []
for bt in tqdm(behavior_transitions):
    sample_df = sample_data.get(bt.sample_id)

    if not any(["bw" in column for column in sample_df.columns]):
        continue
    if sample_df is None:
        raise ValueError("No data found for sample {}".format(bt.sample_id))
    transitions = find_behavior_next(
        bt.sample_id,
        sample_df,
        bt.pre_event,
        bt.event,
        bt.max_delay,
        bt.max_ignored_quiet_time,
    )

    if transitions:
        found_transitions.append(transitions)

# print(len(transitions))
print(len(found_transitions))
print(sum([len(sample_transitions) for sample_transitions in found_transitions]))
# print(found_transitions)
# %%
# Duration between behavior from SameTransition (min, max, avg)
gap_Strans = []
for sample in tqdm(found_transitions):
    for found_transition in sample:
        gap_Strans.append(
            (found_transition["second_event_start"])
            - (found_transition["first_event_end"])
        )

        if (
            (found_transition["second_event_start"])
            - (found_transition["first_event_end"])
        ) < -5:
            print(
                bt.sample_id,
                found_transition["first_event_end"],
                found_transition["second_event_start"],
            )
# print(gap_Strans)

avg_duration = np.mean(gap_Strans)
max_duration = np.max(gap_Strans)
min_duration = np.min(gap_Strans)

print(avg_duration)
print(max_duration)
print(min_duration)

# Histogram (doesnt show everything)
plt.hist(gap_Strans, bins=100, alpha=0.5)
plt.show()

###there is always the same timestamp and -7.354
# %%
# For same_behavior_transition events get event_time and use cell_trace_config to filter by celltype and pattern.
# The results are merged and interpolated.
cell_Strans_configs = []
all_Strans_events = []

for sample in tqdm(found_transitions):
    sample_ls_trans = []
    for found_transition in sample:
        # print(found_transition["sample_id"], found_transition["second_event_start"])
        # print(found_transition["first_event_end"])
        # Sustitute cell_type with None if you don"t want to sort
        sample_ls_trans.append(found_transition["second_event_start"])
        cell_Strans_configs.append(
            CellTransConfig(
                found_transition["sample_id"],
                "A00c",
                found_transition["second_event_start"],
                "mid",
            )
        )

# Extract for specific time window and align several events.
# Define timepoints pre and post an event (event_df).
# This works for single sample or multiple samples aligned
# Note: In cell_subset_df, time was set to index, because for the previous avg calculation
# Add index and time = column

# Set the window range left and right from the event
left_half_window_size = 20.0  # in seconds
right_half_window_size = 50.0

# trans_df defined in pargraph before
windows = []
n_behavior = 0
# TODO: split  window and filter
for ctc in tqdm(cell_Strans_configs):
    sample_df = sample_data.get(ctc.sample_id)
    if sample_df is None:
        raise ValueError("{}: could not find sample data".format(ctc.sample_id))
        continue

    # Extract columns matching our cell type and the optional filter pattern.
    # Pandas' filter() operations works on columns for DataFrames by default.
    # Apply filter regex
    cell_subset_df = sample_df.filter(
        regex=ctc.get_filter_regex()
    )  # Get subset of cells
    cell_subset_df.set_index(
        sample_df.time, inplace=True
    )  # Set time to index (essential for min/max...)
    cell_subset_df.reset_index(inplace=True)  # Add index and time = column
    # print(cell_subset_df)

    # Don't apply filter regex, but take all cells from lm_data
    # cell_subset_df =  lm_data.get(ctc.sample_id)#Get subset of cells  (includes column with sample_id and exp_id)
    # cell_subset_df.set_index(sample_df.time, inplace=True) #Set time to index (essential for min/max...)
    # cell_subset_df.reset_index(inplace = True, drop = True) # Add index and time = column
    # if 'time' not in cell_subset_df:
    #    print("ignoring " + ctc.sample_id)
    #    continue
    # else:
    #    print("processing " + ctc.sample_id)
    # print(cell_subset_df['time'])

    n_behavior += 1
    window_start = ctc.event_time - left_half_window_size
    window_end = ctc.event_time + right_half_window_size

    # Get subset of rows between window_start and window_end
    trans = cell_subset_df[
        (cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)
    ]
    # print(trans)
    # print(trans['time'])
    # print(trans.loc[:,'time'])
    # print(trans['time']- ctc.event_time)

    # Normalizing the data to align on beginning of selected
    # behavior (event_df = Zero) by substracting events in window
    # around start of event of interest from start of event interest.
    # Note: using ":" in event.loc[] will select "all rows" in our window.
    # trans.loc[:, 'time'] = trans['time'] - row['time']
    # trans.loc[:, 'time'] = trans['time'] - ctc.event_time # Works after applying filer
    trans.loc[:, "time"] = trans["time"] - ctc.event_time

    # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
    trans.rename(
        lambda x: "{}_{}_{}".format(ctc.sample_id, x, n_behavior),
        axis="columns",
        inplace=True,
    )

    # Rename time collum to time
    trans.rename(columns={trans.columns[0]: "time"}, inplace=True)
    all_Strans_events.append(trans)  # Append a list with all event
    # print(trans)
# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
all_Strans_df = all_Strans_events.pop(0)
for right_df in all_Strans_events:
    all_Strans_df = pd.merge_ordered(all_Strans_df, right_df, on="time", how="outer")

# Resets the index as time and drops time column (sollte spaeter kommen)
all_Strans_df.index = all_Strans_df["time"]
del all_Strans_df["time"]
# print(all_Strans_df)

# Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
int_all_Strans_df = all_Strans_df.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)
# print(int_all_Strans_df)
# %%
# Average and stddev, min, max, sem for same_behavior_transition events
all_Strans_avg_df = int_all_Strans_df.mean(axis=1)  # Interpolated data used
all_Strans_min_df = int_all_Strans_df.min(axis=1)
all_Strans_max_df = int_all_Strans_df.max(axis=1)
# Standard deviation (distribution)
all_Strans_std_df = int_all_Strans_df.std(axis=1)
# standard error of mean
all_Strans_sem_df = int_all_Strans_df.sem(axis=1)
# wrong zur haelfte: Want to have avg per celltyp over time point,
# and not avg over all cells per timepoint --> should work with transition-grouper

# Plotting for multi-events (same_behavioral_transition)
# If a dataframe with NANs is plotted (raw-data = non interpolated), use
# marker = '+', or 'o', since the line in the lineplot only connects
# consecutive data points
def aligned_layout_plot(plot, tick_spacing=1, fov=(-10, 10, 0.0, 0.2), legend=False):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
# sub1 = fig.add_subplot(111) #211
# all_trans_df.plot(ax=sub1, marker = '*', label = ctc.cell_type)
# aligned_layout_plot(sub1)

sub2 = fig.add_subplot(111)  # 212
all_Strans_avg_df.plot(
    ax=sub2, color="r", label=ctc.cell_type
)  # use interpolated df to calculate average...
# all_Strans_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
# all_Strans_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
all_Strans_avg_df.plot.line(
    yerr=all_Strans_std_df, ax=sub2, color="lightgrey", alpha=0.5
)
# all_Strans_avg_df.plot.line(yerr=all_Strans_sem_df, ax=sub2, color = 'grey', alpha = 0.1)
aligned_layout_plot(sub2)
# %%
# Using class Transitiongrouper to group cell_type/pattern or both
# This is crucial if statistical methods are applied for sub groups
# Input dataframe from behavioral transitions (Post-, Same-, Triple-Transition)

# grouper = TransitionGrouper(int_all_Ptrans_df)
groupers = [TransitionGrouper(int_pre_data), TransitionGrouper(int_post_data)]


# Function to calculate the average cell/pattern_groups depending on regex
def average_grouping(grouping):
    df = grouping[list(grouping)[0]][1]
    average_df = pd.DataFrame(index=df.index)

    for group_pattern, df in grouping.values():
        average_col = df.mean(axis=1)
        average_df[group_pattern] = average_col

    return average_df


# Function to calculate the std cell/pattern_groups depending on regex
def std_grouping(grouping):
    df = grouping[list(grouping)[0]][1]
    std_df = pd.DataFrame(index=df.index)

    for group_pattern, df in grouping.values():
        std_col = df.std(axis=1)
        std_df[group_pattern] = std_col

    return std_df


# Function to calculate the sem cell/pattern_groups depending on regex
def sem_grouping(grouping):
    df = grouping[list(grouping)[0]][1]
    sem_df = pd.DataFrame(index=df.index)

    for group_pattern, df in grouping.values():
        sem_col = df.sem(axis=1)
        sem_df[group_pattern] = sem_col

    return sem_df


cell_groups = []
cell_averages = []
cell_stds = []
cell_sems = []

pattern_groups = []
pattern_averages = []
pattern_stds = []
pattern_sems = []

cellpattern_groups = []
cellpattern_averages = []
cellpattern_stds = []
cellpattern_sems = []

for grouper in groupers:
    cell_group = grouper.group_cells()
    cell_groups.append(cell_group)
    cell_averages.append(average_grouping(cell_group))
    cell_stds.append(std_grouping(cell_group))
    cell_sems.append(sem_grouping(cell_group))
    # print(cell_averages.head())

    pattern_group = grouper.group_patterns()
    pattern_groups.append(pattern_group)
    pattern_averages.append(average_grouping(pattern_group))
    pattern_stds.append(std_grouping(pattern_group))
    pattern_sems.append(sem_grouping(pattern_group))
    # print(pattern_averages.head())

    cellpattern_group = grouper.group_cellpattern()
    cellpattern_averages.append(average_grouping(cellpattern_group))
    cellpattern_stds.append(std_grouping(cellpattern_group))
    cellpattern_sems.append(sem_grouping(cellpattern_group))
    # print(pattern_averages)

# %%
# Plot average, std-dev from Transition grouper
for i, (cell_average, cell_sem, pattern_average) in enumerate(
    zip(cell_averages, cell_sems, pattern_averages)
):

    def aligned_layout_plot(
        plot, tick_spacing=1, fov=(-10, 10, 0.0, 0.9), legend=False
    ):
        # Set fine x-axis scale
        plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        # Set x and y limits and legend (default = False)
        plot.axis(fov)
        plot.legend().set_visible(legend)

    fig = plt.figure(i)

    sub1 = fig.add_subplot(211)
    pattern_average.plot(ax=sub1)
    aligned_layout_plot(sub1)

    sub2 = fig.add_subplot(212)
    cell_average.plot(ax=sub2, color="r")
    # cell_average.plot(ax=sub2, yerr=cell_std, color = 'lightgrey', alpha = 0.1)
    cell_average.plot(ax=sub2, yerr=cell_sem, color="grey", alpha=0.5)
    aligned_layout_plot(sub2)
# %%
# Duration of each behavior (min, max, avg)

desired_behavior = "hunch"
all_durations = []
sample_ids = []
start_times = []
end_times = []
for sample_id, sample_df in tqdm(sample_data.items()):
    if any(["bw" in column for column in sample_df.columns]):
        for i, row in sample_df.iterrows():
            if row["{}_start".format(desired_behavior)]:
                start_times.append(row["time"])
                sample_ids.append(sample_id)
            if row["{}_end".format(desired_behavior)]:
                end_time = row["time"]
                end_times.append(row["time"])
        assert len(end_times) == len(
            start_times
        ), "{} start times vs {} end times on sample {}".format(
            len(start_times), len(end_times), sample_id
        )
all_durations = [end - start for start, end in zip(start_times, end_times)]
durations_with_sample_id = list(zip(sample_ids, all_durations))

# Test if behavior.csv is correct
for t, sample_id, start_time, end_time in zip(
    all_durations, sample_ids, start_times, end_times
):
    if t > 200:
        print(sample_id, start_time, end_time)

# print(durations_with_sample_id[0:10])

# print(all_durations)
# print(len(all_durations))
avg_duration = np.mean(all_durations)
max_duration = np.max(all_durations)
min_duration = np.min(all_durations)

print(avg_duration)
# print(max_duration)
# print(min_duration)

# Histogram
fig = plt.figure()
plt.hist(all_durations, bins=5, alpha=0.5)
plt.show()

# print(list(filter(lambda x: x> 10, all_durations)))  # not working
# %%
# Neuronal activity of a certain behavior (all behavior from the same type, independent on pre- and post-event)

# For all behavior events (same type) get event_time and use cell_trace_config to filter by celltype and pattern.

# The results are merged and interpolated.

cell = "A00c"
# pattern = 'mid'
desired_behavior = "fw"
align_to_start = True

sample_ids = []
start_times = []
end_times = []
desired_traces = []

# Checkpoint if not all behavior data are used for analysis, but dff and timestamps are available
for sample_id, sample_df in tqdm(sample_data.items()):
    if any(["bw" in column for column in sample_df.columns]):
        for i, row in sample_df.iterrows():
            if row["{}_start".format(desired_behavior)]:
                start_times.append(row["time"])
                sample_ids.append(sample_id)
            if row["{}_end".format(desired_behavior)]:
                end_time = row["time"]
                end_times.append(row["time"])
        assert len(end_times) == len(
            start_times
        ), "{} start times vs {} end times on sample {}".format(
            len(start_times), len(end_times), sample_id
        )

for sample_id, start_time, end_time in zip(sample_ids, start_times, end_times):
    if align_to_start:
        desired_traces.append(CellTransConfig(sample_id, cell, start_time))
    else:
        desired_traces.append(CellTransConfig(sample_id, cell, end_time))

print(len(desired_traces))
# %%
# For testing if samples are different from each other, extract per sample the peak of each behavior wave (eg fw),
# for each A00c, each hemisegment for Basin and each segment for handles. No averaging of cells.

cell = (
    "A00c"
)  # Note: for each cell type we need a specific grouping because of missing identity
pattern = "post"
desired_behavior = "fw"

sample_ids = []
peak_values = []
start_times = []
start_index = []
end_times = []
end_index = []

# Checkpoint if not all behavior data are used for analysis, but dff and timestamps are available
for sample_id, sample_df in tqdm(sample_data.items()):
    if any(["bw" in column for column in sample_df.columns]):
        for i, row in sample_df.iterrows():
            if row["{}_start".format(desired_behavior)]:
                start_times.append(row["time"])
                start_index.append(i)
                sample_ids.append(sample_id)
            if row["{}_end".format(desired_behavior)]:
                end_time = row["time"]
                end_index.append(i)
                end_times.append(row["time"])
        assert len(end_times) == len(
            start_times
        ), "{} start times vs {} end times on sample {}".format(
            len(start_times), len(end_times), sample_id
        )


max_values = {}
for sample_id, start_time, end_time in zip(sample_ids, start_index, end_index):
    max_values.setdefault(sample_id, [])
    regex_pattern = "{}_{}".format(cell, pattern)
    sample_df = sample_data.get(sample_id)
    max_values[sample_id].append(
        sample_df.filter(like=regex_pattern, axis=1).loc[start_time:end_time].max()
    )

print({key: len(values) for key, values in max_values.items()})
print(max_values)

# %%
from vnc_stats_python.hypothesis_tests.permutation_tests import all_pairs_perm_tests

# Pull out max values for each sample, taking care to throw out any nan values
smp_labels = list(max_values.keys())
n_smps = len(smp_labels)
smp_max_vls = [None] * n_smps
for i, s_l in enumerate(smp_labels):
    vls = np.asarray(max_values[s_l])
    vls = vls[np.logical_not(np.isnan(vls))]  # Get rid of nan values
    smp_max_vls[i] = vls


# Now we want to get rid of any samples where we had *no* data (after removing nan values)
good_inds = [len(v) != 0 for v in smp_max_vls]
smp_labels = [s_l for i, sl in enumerate(smp_labels) if good_inds[i]]
smp_max_vls = [vl for i, vl in enumerate(smp_max_vls) if good_inds[i]]

p_vls = all_pairs_perm_tests(smp_max_vls, test_opts={"n_perms": 10000}, update_int=100)
# %%
plt.imshow(p_vls < 0.01)
print(p_vls)
# plt.colorbar()
# %%
# Correlation Matrix: cell/cell before and cell/cell after transition

# Extract for data int_all_Ptrans where second_event_start - second event ends
# I cheat with window size (its actually a good idea:)
# legend would be nice
# header of matrix
# muss average pro cell_id haben sonst plottet er fur jedes sample/transition extra
# average over cells not time should work

# TEST
# print(int_pre_data)
# print(int_post_data)

# Using TransitionGrouper
# Be careful where you generate the dataframe
# TODO after grouper I lose pre and post window
# grouper = [TransitionGrouper(int_pre_data),
#           #TransitionGrouper(int_post_data)
#          ] NOT WORKING + OVERWRITES WITH THE PREVIOUS ONE


# cell_groups = grouper.group_cells()
# cell_averages = average_grouping(cell_groups)
# print(cell_averages)

# pattern_groups = grouper.group_patterns()

# cellpattern_groups = grouper.group_cellpattern()
# cellpattern_averages = average_grouping(cellpattern_groups)
# print(pattern_averages)
print(cellpattern_averages)

# Define time window for pre-transition_start
# pre_data defined before
# Define time window for post-transition_start
# post_data defined before

plt.matshow(cellpattern_averages.corr())
# plt.matshow(int_post_data.corr())
# plt.matshow(int_post_data.corr()-int_pre_data.corr())
cb = plt.colorbar()
cb.ax.tick_params(labelsize=14)
plt.title("Correlation Matrix", fontsize=16)
plt.show()
correlation = np.corrcoef(cellpattern_averages)  # data(cell,time)
print(correlation)
# %%
import seaborn as sns

corr = cellpattern_averages.corr()
ax = sns.heatmap(
    corr,
    vmin=-1,
    vmax=1,
    center=0,
    cmap=sns.diverging_palette(20, 220, n=200),
    square=True,
)

# size(c1, c2) ~ abs(corr(c1, c2)) # should make the visualisation better (not tested yet)
# %%
# Test if sum of left and right turns = turns
# for sample_id, sample_df in tqdm(sample_data.items()):
# print(sample_id, sum(sample_df['turn_start']), sum(sample_df['left turn_start']), sum(sample_df['right turn_start']))
#    if sum(sample_df['turn_start']) != (sum(sample_df['left turn_start']) + sum(sample_df['right turn_start'])):
#        print(sample_id)
# %%
# # Test number of behavior
# alle = []
# for sample_id, sample_df in tqdm(sample_data.items()):
# print(sum(sample_df['turn_start']))
#    alle.append(sum(sample_df['HP_start']))

# print(sum(alle))

