#%%
from read_functions import readall_behavior, readall_lm, readall_time

import glob
import os
import pandas as pd
import ast

import pickle

import logging

# Logging info:
# Logging has 4 levels (DEBUG, INFO, WARNING, ERROR)
# appends to log file rather than overwritting it
logging.basicConfig(level=logging.INFO, filename="log.log")

# Implement error messages (Default should be True)
ERROR_ON_MISSING_TIMESTAMPS = False
ERROR_ON_TIME_LIGHT_MISMATCH = False
ERROR_ON_TIME_BEHAVIOR_MISMATCH = False
ERROR_ON_MISSING_BEHAVIORS = False
ERROR_ON_INVALID_BEHAVIOR_RANGE = False

# available behaviors
AVAILABLE_BEHAVIORS = (
    "fw",
    "bw",
    "stim",
    "hunch",
    "turn",
    "other",
    "HP",
    "left turn",
    "right turn",
)

# Files
BEHAVIOR_DIRECTORIES = [
    r"/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_cl_A4/",
    # r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_cl_A9/',
    # r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_ol/',
    # r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_stim_artefact/'
]

LIGHTMICROSOPE_DIRECTORIES = [
    r"/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/Basin_traces/",
    r"/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/Handle-like_Traces",
    r"/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/a00c_traces",
    # r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/candidate_neuron_traces'
]

TIMELAPSE_DIRECTORY = r"/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/timelapse/"

# If True: use time cache (create if it doesnt exist)
# If False: create time cache whether it exists or not
USE_TIME_CACHE = True

#%%
# read in behavior data

logging.warning(
    "Chris's FRAN is zero based!!, whereas old annotations are 1 based. NOT YET HANDLED PROPERLY?"
)

behavior_files = []
for d in BEHAVIOR_DIRECTORIES:
    behavior_files.extend(
        glob.glob(os.path.join(d, "*.csv"))
    )  # join pathname with filename,

behavior_data = readall_behavior(behavior_files)


#%%
# Concatenate behavior_data dataframes. Get behavior counts.
# TODO: this should probaly be put in a funciton since df_behavior is not reused
df_behavior = pd.concat(
    behavior_data.values(), axis=0, ignore_index=True, sort=False
)  # add sorting
logging.info("Behavior counts:\n{}".format(df_behavior[df_behavior == 1].count()))

#%%
# Import and merge fluorescence data: Several LM files for the same sample_id exists, but differ in cell_id).
# List of LM data with two extra columns: sample_id and cell_id
# Open LM files from different directories

# Iterate through LM data and extend files in a list from within and between directory and
# build a list of files from all directories
# (Note: append would 'extend' lists and not single files)
lightmicroscope_files = []
for d in LIGHTMICROSOPE_DIRECTORIES:
    lightmicroscope_files.extend(
        glob.glob(os.path.join(d, "*.csv"))
    )  # join pathname with filename


lm_samples = readall_lm(lightmicroscope_files)

# New dictionary: lm_data{} to build a single dataframe with all cell_ids combined
# for a single sample. Iterate over dict from same sample in one dataframe.
# df.items iterate over pairs and build a list

lm_data = {}

# Iterate over all light samples and merge all found files
# for each sample into a single data frame (per sample)
for sample_id, sample_info in lm_samples.items():
    cells_dataframes = sample_info["data"]
    # check if number of cells >= 1
    if not cells_dataframes:
        raise ValueError("No cells found for sample {}".format(sample_id))
    # first element in the list
    lm_df = None

    # iteration through other df
    for cdf in cells_dataframes:
        if lm_df is None:
            lm_df = cdf
        else:
            if len(lm_df.index) != len(cdf.index):
                raise ValueError(
                    "Data frame frame to merge has not same row count as target",
                    sample_id,
                )
            lm_df = pd.merge(lm_df, cdf, left_index=True, right_index=True)

    lm_df["sample_id"] = sample_id  # add sample_id column
    lm_df["exp_id"] = sample_info["exp_id"]
    lm_data[sample_id] = lm_df
# print(list(lm_data.keys()))

#%%
# Import txt-files from of the absolute time/frame from the Ca-imaging (lm-data).
# All txt-files have to be transposed, which is a memory intensive step. After the
# data are complete, the transposed files should be exported (ToDo). Time-data are
# combined with sample-ID and experiment-ID.

timelapse_files = glob.glob(
    os.path.join(TIMELAPSE_DIRECTORY, "*.txt")
)  # join pathname with filename

# Keep in mind that some of the files has to be changed because of the discrepancy in timestamps
timelapse_cache = "timelapse.cache"

try:
    if not USE_TIME_CACHE:
        raise FileNotFoundError()
    with open(timelapse_cache, "rb") as timelapse_cache_file:
        # TODO
        time_data = pickle.load(timelapse_cache_file)
        # cache_data = timelapse_cache_file.read()
        # time_data = ast.literal_eval(cache_data)
except FileNotFoundError as e:
    print("No cache file found, recomputing")
    # No cache file found, recompute
    time_data = readall_time(timelapse_files)
    # Write cache
    with open(timelapse_cache, "wb") as timelapse_cache_file:
        pickle.dump(time_data, timelapse_cache_file)


#%%

# TODO: This should be in a function

# Time data are merged into light data and checked if number length of lm = timestamp.
# Due to technical conditions, some time.txt-file have too many or not enough time data compared
# to the corresponding LM data. The discrepancy is fixed by either dropping the extra timepoints or
# by taking the average of the difference between each timepoint and extend the dataframe.
# The first 10 timepoints are not included to account for instability of the microscope in
# the beginning due to the moving parts.
# Maximal difference between timepoints fyi.
from functions import combine_lm_time_data

sample_data = lm_data
combine_lm_time_data(
    sample_data, time_data, ERROR_ON_MISSING_TIMESTAMPS, ERROR_ON_TIME_LIGHT_MISMATCH
)

logging.info(
    "Matched {} light data sets with their respective time points".format(
        len(sample_data)
    )
)

#%%
# TODO: This should also be in a function
# Combine behavior data with light data into a single data frame
# per sample ID. To do so, add behavior data to light data frames,
# because the light data is already organizes by frame. To accomodate
# frame ranges without an behavior data, a column named "quiet" is
# added which is True in these cases and False otherwise. Additionally,
# for each behavior column, a behavior start and end column as well as
# an overlap column is added so that parallel and successive behaviors
# of the same type can be differentiated.
from functions import combine_lm_behavior_data

combine_lm_behavior_data(
    sample_data,
    behavior_data,
    AVAILABLE_BEHAVIORS,
    ERROR_ON_MISSING_BEHAVIORS,
    ERROR_ON_INVALID_BEHAVIOR_RANGE,
    ERROR_ON_TIME_BEHAVIOR_MISMATCH,
)

raise Exception("Stop HERE!")


#%%
#########################################################################################################################
# DATA ANALYSIS
#########################################################################################################################
#%%

# TODO: only take all samples into account if event is quiet (maybe make event optional)

# Determine how many cells of the same cell type are in each (abdominal) hemisegment per sample
# Basins:
# cell_trace_configs = [
#    CellTraceConfig('17-09-01L3-cl', 'basin', 'quiet', 'A2'),
#    CellTraceConfig('17-08-26L6-cl', 'basin', 'quiet', 'A2')
# ]

cell_trace_configs = [
    CellTraceConfig(name, "basin", "quiet", "A4L") for name in lm_data
]  # if all samples should be included, you have
# to chose 'quiet' as event


number_basins = []
sample_ids = []
for ctc in cell_trace_configs:
    sample_df = sample_data.get(ctc.sample_id)
    if sample_df is None:
        raise ValueError("{}: could not find sample data".format(ctc.sample_id))
        continue
    # Extract columns matching our cell type and the optional filter pattern.
    cell_subset_df = sample_df.filter(
        regex=ctc.get_filter_regex()
    )  # Get subset of cells
    number = len(cell_subset_df.columns)
    number_basins.append(number)
    sample_ids.append(ctc.sample_id)
avg = np.mean(number_basins)
total = np.sum(number_basins)
print(total)
# print(len(number_basins))
# print(list(zip(sample_ids, number_basins)))
print(avg)

#%%
# TODO: this should probably be a function (extract time windows)
# TODO: check redundency

# Load single samples with specific filter pattern. A certain event is found with TP-start
# and marked with vertical line. No transition events are considered! Data can be merged along
# the event and interpolated. Average happens for all! columns.
cell_trace_configs = [
    CellTraceConfig("17-09-01L3-cl", "A00c", "stim", "mid"),
    CellTraceConfig("17-08-26L6-cl", "A00c", "stim", "mid"),
    CellTraceConfig("17-08-28L3-cl", "A00c", "stim"),
]

# Load all samples with specific filter pattern
"""
cell_trace_configs = [
    CellTraceConfig(name,'A00c', 'fw') for name in lm_data]
"""

all_events = []  # List of events, with raw dff data (no interpolation or other
# processing done at this point). Sample_id is added to the cell name.

for ctc in cell_trace_configs:
    # print(sample_df.keys)
    # break
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
    # print(ctc.sample_id, cell_subset_df)
    # Get rows where current event starts.
    event_df = sample_df[sample_df.loc[:, ctc.get_event_start_col()]]
    # Gives the timestamp for the event_df (start)
    for idx, row in event_df.iterrows():
        # print('TP of {} ='.format(ctc.event_name), row['time'])
        print(
            "TP of {} =".format(ctc.event_name), row["time"], "{}".format(ctc.sample_id)
        )

    # Extract for specific time window and align several events.
    # Define timepoints pre and post an event (event_df).
    # This works for single sample or multiple samples aligned
    # Note: In cell_subset_df, time was set to index, because for the previous avg calculation
    # Add index and time = column

    # Set the window range left and right from the event
    left_half_window_size = 10.0  # in seconds
    right_half_window_size = 50.0

    # Event_df defined in pargraph before
    windows = []
    n_behavior = 0
    for i, row in event_df.iterrows():
        n_behavior += 1
        window_start = row["time"] - left_half_window_size
        window_end = row["time"] + right_half_window_size

        # Get subset of rows between window_start and window_end
        event = cell_subset_df[
            (cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)
        ]
        # Normalizing the data to align on beginning of selected
        # behavior (event_df = Zero) by substracting events in window
        # around start of event of interest from start of event interest.
        # Note: using ":" in event.loc[] will select "all rows" in our window.
        event.loc[:, "time"] = event["time"] - row["time"]

        # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
        event.rename(
            lambda x: "{}_{}_{}".format(ctc.sample_id, x, n_behavior),
            axis="columns",
            inplace=True,
        )

        # Rename time collum to time
        event.rename(columns={event.columns[0]: "time"}, inplace=True)
        all_events.append(event)  # Append a list with all event

        # Round (NR)
        # decimals = 1
        # event['time'] = event['time'].apply(lambda x: round(x, decimals))

print(len(all_events))

# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
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

#%%

"""
    # Single sample analysis
    # For single sample over the whole experimental time
    # Calculate min, max, avg, stddev, sem from cell_subset_df (defined earlier)
    cell_subset_df.set_index(
        sample_df.time, inplace=True
    )  # Set time to index (essential for min/max...)
    del cell_subset_df["time"]  # delete time_column
    cell_subset_df.index.name = None  # delete index name
    cell_avg_df = cell_subset_df.mean(axis=1)
    cell_min_df = cell_subset_df.min(axis=1)
    cell_max_df = cell_subset_df.max(axis=1)
    # Standard deviation (distribution)
    cell_std_df = cell_subset_df.std(axis=1)
    # standard error of mean
    cell_sem_df = cell_subset_df.sem(axis=1)
    # print(ctc.sample_id, cell_avg_df)
"""


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


#%%
# Plotting - single sample for whole time
def layout_plot(plot, tick_spacing=100, fov=(0, 2400, 0, 1.2), legend=False):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


# Get rows where current event is active and draw a vertical
# line to indicate the event in the plot
event_df = sample_df[sample_df.loc[:, ctc.get_event_start_col()] == 1]
fig = plt.figure()
fig.set_facecolor("white")

# Plot all cells from cell_subset_df over entire time (specified in Cell_Trace_Config). <For a single sample>
sub1 = fig.add_subplot(211)  # 211
cell_subset_df.plot(ax=sub1)
ctc.add_event_time_points_to_plot(event_df, sub1)
layout_plot(sub1)

# Avg, min, max, std-dev for multiple cells in single sample over whole time.
# All cells are averaged. For cell and filterpattern see below.
sub2 = fig.add_subplot(212)  # 212
ctc.add_event_time_points_to_plot(event_df, sub2)
cell_avg_df.plot(ax=sub2, color="g", label=ctc.cell_type, linewidth=1)
# cell_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
# cell_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
cell_avg_df.plot.line(yerr=cell_std_df, ax=sub2, color="r", alpha=0.1)
# cell_avg_df.plot.line(yerr=cell_sem_df, ax=sub2, color = 'c', alpha = 0.1)
layout_plot(sub2)

#%%

# PLOTTING THE ALIGNED EVENT, INDEPENDENT OF TRANSITION EVENT
# Plotting for multi-events (all_df) (raw_dff_data) of a single sample
# If a dataframe with NANs is plotted, use
# marker = '+', or 'o', since the line in the lineplot only connects
# consecutive data points
def aligned_layout_plot(
    plot, tick_spacing=0.5, fov=(-20, 50, -0.05, 1.9), legend=False
):
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False)
    plot.axis(fov)
    plot.legend().set_visible(legend)


fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
sub1 = fig.add_subplot(211)
all_df.plot(ax=sub1, marker="*", label=ctc.cell_type)
aligned_layout_plot(sub1)

sub2 = fig.add_subplot(212)
all_cell_avg_df.plot(
    ax=sub2, color="k", label=ctc.cell_type
)  # use interpolated df to calculate average...
# all_cell_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
# all_cell_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
all_cell_avg_df.plot.line(yerr=all_cell_std_df, ax=sub2, color="r", alpha=0.1)
# all_cell_avg_df.plot.line(yerr=all_cell_sem_df, ax=sub2, color = 'c', alpha = 0.1)
aligned_layout_plot(sub2)
