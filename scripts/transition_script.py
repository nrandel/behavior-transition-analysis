# %%
# Imports and Load data

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from tqdm import tqdm

from logging_functions import log_num_transitions
from logging_functions import log_unique_column_values
from logging_functions import log_mann_whitney_u_test

from plotting_functions import plot_transition_gaps_hist
from plotting_functions import plot_transitions
from plotting_functions import plot_fold_change

from transition_functions import extract_transitions
from transition_functions import extract_transition_triples

from functions import merge_dataframe_list
from functions import get_sample_data
from functions import get_fold_change_df
from functions import extract_transition_dfs

from classes.transition_type import TransitionType
from classes.cell_transition_config import CellTransConfig, extract_windows
from classes import PostBehaviorTransition

# Logging info:
# Logging has 4 levels (DEBUG, INFO, WARNING, ERROR)
# appends to log file rather than overwritting it
logging.basicConfig(level=logging.INFO, filename="log.log")

# TODO: Cleanup comments in get_sample_data
sample_data = get_sample_data()

PLOT_EVERY_SAMPLE = True
PLOT_TRANSITION_GAPS = True
USE_CELL_PATTERN_FILTER = True

FILTER_PATTERN = "A1R"

"""
BEHAVIOR_TRANSITIONS = [
    PostBehaviorTransition("17-08-26L2-cl", "turn", "fw", 4.9),
    PostBehaviorTransition("17-08-26L6-cl", "turn", "fw", 4.9),
]

"""
# Open all samples multiple Transitions
BEHAVIOR_TRANSITIONS = (
    [PostBehaviorTransition(name, "fw", "stim", 4.9) for name in sample_data]
    + [PostBehaviorTransition(name, "bw", "stim", 4.9) for name in sample_data]
    + [PostBehaviorTransition(name, "turn", "stim", 4.9) for name in sample_data]
    + [PostBehaviorTransition(name, "hunch", "stim", 4.9) for name in sample_data]
    + [PostBehaviorTransition(name, "other", "stim", 4.9) for name in sample_data]
)

BEHAVIOR_TRANSITIONS = (
    [PostBehaviorTransition("17-08-24L4-cl", "bw", "stim", 1000)]
)

TRANSITION_TYPES = [
    TransitionType(first_event="fw"),
    TransitionType(first_event="bw"),
    TransitionType(first_event="turn"),
    TransitionType(first_event="hunch"),
    TransitionType(first_event="other"),
    TransitionType(first_event="stim"),
]

MANN_WHITNEY_U_TEST_COMPARISONS = (
    ("fw_stim", "bw_stim"),
    ("fw_stim", "turn_stim"),
    ("fw_stim", "hunch_stim"),
    ("fw_stim", "other_stim"),
    ("bw_stim", "turn_stim"),
    ("bw_stim", "bw_stim"),
)


# %%
#######################################################################################################################
# transitions
#######################################################################################################################


found_transitions = extract_transitions(sample_data, BEHAVIOR_TRANSITIONS)

log_num_transitions(found_transitions)


if PLOT_TRANSITION_GAPS:
    plot_transition_gaps_hist(found_transitions)


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
                filter_pattern=FILTER_PATTERN,
                first_event=found_transition["first_event"],
                second_event=found_transition["second_event"],
            )
        )

        # For stimulus as first_event
        # cell_Ptrans_configs.append(CellTransConfig(found_transition["sample_id"], "basin",
        #                                           found_transition["first_event_start"], 'A6',
        #                                           first_event=found_transition["first_event"],
        #                                           second_event=found_transition["second_event"]))

all_Ptrans_events = extract_windows(
    sample_data, cell_Ptrans_configs, cell_pattern_filter=USE_CELL_PATTERN_FILTER
)

all_Ptrans_df = merge_dataframe_list(all_Ptrans_events, on="time", how="outer")

# Resets the index as time and drops time column
all_Ptrans_df.index = all_Ptrans_df["time"]
del all_Ptrans_df["time"]


int_all_Ptrans_df = all_Ptrans_df.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)

log_unique_column_values(
    int_all_Ptrans_df, cell=True, filter_pattern=True, sample_id=True
)


plot_transitions(all_Ptrans_df, TRANSITION_TYPES, use_sem=False)
plot_transitions(all_Ptrans_df, TRANSITION_TYPES, use_sem=True)

plot_transitions(int_all_Ptrans_df, TRANSITION_TYPES, use_sem=False)
plot_transitions(int_all_Ptrans_df, TRANSITION_TYPES, use_sem=True)

fold_change_df = get_fold_change_df(df=all_Ptrans_df, pre=0, post=2)
transitions, transition_dfs = extract_transition_dfs(fold_change_df)
all_fold_change_df = merge_dataframe_list(
    transition_dfs, ordered=False, left_index=True, right_index=True, how="outer"
)
plot_fold_change(df=all_fold_change_df, transitions=transitions)

int_fold_change_df = get_fold_change_df(df=int_all_Ptrans_df, pre=0, post=2)
int_transitions, int_transition_dfs = extract_transition_dfs(int_fold_change_df)
int_all_fold_change_df = merge_dataframe_list(
    int_transition_dfs, ordered=False, left_index=True, right_index=True, how="outer"
)
plot_fold_change(df=int_all_fold_change_df, transitions=int_transitions)

# %%
# Mann-Witney U test (non-parametric pair wise test)

# da ist ne Menge significant:()

log_mann_whitney_u_test(transition_dfs, transitions, MANN_WHITNEY_U_TEST_COMPARISONS)
log_mann_whitney_u_test(
    int_transition_dfs, int_transitions, MANN_WHITNEY_U_TEST_COMPARISONS
)

raise Exception("STOP HERE!")

####################################################################################################
# GROUPINGS
####################################################################################################
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

