from tqdm import tqdm
import logging
from read_functions import readall_behavior, readall_lm, readall_time

import glob
import os
import pandas as pd
import ast

import pickle


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
    r"/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_cl_A4/",
    # r'/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_cl_A9/',
    # r'/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_ol/',
    # r'/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/behavior_csv_stim_artefact/'
]

LIGHTMICROSOPE_DIRECTORIES = [
    r"/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/Basin_traces/",
    r"/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/Handle-like_Traces",
    r"/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/a00c_traces",
    # r'/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/candidate_neuron_traces'
]

TIMELAPSE_DIRECTORY = r"/Users/Nadine/Documents/Zlatic_lab/close-loop/Notes/timelapse/"

# If True: use time cache (create if it doesnt exist)
# If False: create time cache whether it exists or not
USE_TIME_CACHE = True
TIMELAPSE_CACHE_FILE = "timelapse.cache"


def get_sample_data():
    logging.warning(
        "Chris's FRAN is zero based!!, whereas old annotations are 1 based. NOT YET HANDLED PROPERLY?"
    )

    behavior_files = []
    for d in BEHAVIOR_DIRECTORIES:
        behavior_files.extend(glob.glob(os.path.join(d, "*.csv")))

    lightmicroscope_files = []
    for d in LIGHTMICROSOPE_DIRECTORIES:
        lightmicroscope_files.extend(glob.glob(os.path.join(d, "*.csv")))

    timelapse_files = glob.glob(os.path.join(TIMELAPSE_DIRECTORY, "*.txt"))

    behavior_data = readall_behavior(behavior_files)
    lm_data = readall_lm(lightmicroscope_files)
    time_data = readall_time(
        timelapse_files,
        timelapse_cache=TIMELAPSE_CACHE_FILE,
        use_time_cache=USE_TIME_CACHE,
    )

    sample_data = lm_data
    combine_lm_time_data(
        sample_data,
        time_data,
        ERROR_ON_MISSING_TIMESTAMPS,
        ERROR_ON_TIME_LIGHT_MISMATCH,
    )

    combine_lm_behavior_data(
        sample_data,
        behavior_data,
        AVAILABLE_BEHAVIORS,
        ERROR_ON_MISSING_BEHAVIORS,
        ERROR_ON_INVALID_BEHAVIOR_RANGE,
        ERROR_ON_TIME_BEHAVIOR_MISMATCH,
    )

    return sample_data


def combine_lm_behavior_data(
    lm_data,
    behavior_data,
    available_behaviors,
    error_on_missing_behaviors=True,
    error_on_invalid_behavior_range=True,
    error_on_time_behavior_mismatch=True,
):
    """
    Combine behavior data with light data into a single data frame
    per sample ID. To do so, add behavior data to light data frames,
    because the light data is already organizes by frame. To accomodate
    frame ranges without an behavior data, a column named "quiet" is
    added which is True in these cases and False otherwise. Additionally,
    for each behavior column, a behavior start and end column as well as
    an overlap column is added so that parallel and successive behaviors
    of the same type can be differentiated.
    """
    for sample_id, lm_df in tqdm(lm_data.items(), "Combine Light with Behavior: "):
        behavior_df = behavior_data.get(sample_id)
        if behavior_df is None:
            msg = 'Could not find behavior data for sample "{}"'.format(sample_id)
            if error_on_missing_behaviors:
                raise ValueError(msg)
            logging.warning(msg)
            continue

        # Add extra columns for behavior
        for behavior in available_behaviors:
            lm_df[behavior] = False
            lm_df["{}_start".format(behavior)] = False
            lm_df["{}_end".format(behavior)] = False
            lm_df["{}_overlap".format(behavior)] = False

        # Add 'quiet' column. Set it initially to True and mark frames
        # with actual behavior as quiet = False.
        lm_df["quiet"] = True

        n_light_entries = len(lm_df)

        # Iterate over behavior data and add data to target data frame
        for i, row in behavior_df.iterrows():
            # Start and end are 1-based, make them 0-based
            start = int(row["START"])
            end = int(row["END"])

            if type(row["START"]) == str:
                logging.debug(
                    "sample_id: {}, starts at {}, ends at {}".format(
                        sample_id, start, end
                    )
                )

            if start >= end:
                msg = "{}: start ({}) needs to be strictly smaller than end ({})".format(
                    sample_id, start, end
                )
                if error_on_invalid_behavior_range:
                    raise ValueError(msg)
                logging.warning(msg)
                continue

            # Make sure we capture start/end times that are a fractional number.
            if row["START"] - start > 0 or row["END"] - end > 0:
                raise ValueError(
                    "{}: start and end frame number can't contain fractions".format(
                        sample_id
                    )
                )

            # Ignore behavior entries with an end frame higher than available light data.
            # The behavior data is one-based, which is why a strict larger than test should
            # be correct.
            if end > n_light_entries:
                msg = "Sample: {} - Behavior row with range {}-{} exceeds light time points ({}): {}".format(
                    sample_id, start, end, n_light_entries, row
                )
                if error_on_time_behavior_mismatch:
                    raise ValueError(msg)
                logging.warning(msg)
                continue

            # Find behavior observed in row
            observed_behaviors = []
            for behavior in available_behaviors:
                if row[behavior]:
                    observed_behaviors.append(behavior)

            # We assume that not more than two behaviors are observed at the same time
            if len(observed_behaviors) > 2:
                raise ValueError(
                    "Found multiple behaviors in row {} of sample {}".format(
                        i, sample_id
                    )
                )

            # Add observed behavior information to target data frames in all
            # rows in behavior range.
            for b in observed_behaviors:
                # Iterate over frames valid for current behavior. Every valid
                # frame is mapped into the canonical (light/cell) data frame,
                # which is 0-indexed.
                for j in range(start, end + 1):
                    # Behavior ranges are 1-indexed
                    current_frame = j - 1
                    # If the current behavior has already been observed at this frame,
                    # set overlap to True, because we are about to mark this behavior
                    # again as observed for this frame.
                    if lm_df.at[current_frame, b]:
                        lm_df.at[current_frame, "{}_overlap".format(b)] = True
                    else:
                        lm_df.at[current_frame, b] = True

                    # Mark this row as not quiet, because we observed
                    # a behavior in the current frame.
                    lm_df.at[current_frame, "quiet"] = False

                lm_df.at[start - 1, "{}_start".format(b)] = True
                lm_df.at[end - 1, "{}_end".format(b)] = True

        # Mark quiet ranges with _start, _end and _overlap. By definion,
        # quiet_overlap is always False.
        lm_df["quiet_start"] = False
        lm_df["quiet_end"] = False
        lm_df["quiet_overlap"] = False
        last_sample_idx = n_light_entries - 1
        for i, row in lm_df.iterrows():
            lm_df.at[i, "quiet_start"] = row["quiet"] and (
                i == 0 or not lm_df.at[i - 1, "quiet"]
            )
            lm_df.at[i, "quiet_end"] = row["quiet"] and (
                i == last_sample_idx or not lm_df.at[i + 1, "quiet"]
            )


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def combine_lm_time_data(
    lm_data,
    time_data,
    error_on_missing_timestamps=True,
    error_on_time_light_mismatch=True,
    show_diffs=False,
):
    # Time data are merged into light data and checked if number length of lm = timestamp.
    # Due to technical conditions, some time.txt-file have too many or not enough time data compared
    # to the corresponding LM data. The discrepancy is fixed by either dropping the extra timepoints or
    # by taking the average of the difference between each timepoint and extend the dataframe.
    # The first 10 timepoints are not included to account for instability of the microscope in
    # the beginning due to the moving parts.
    # Maximal difference between timepoints fyi.
    all_diffs = np.array([])

    for sample_id, lm_df in tqdm(lm_data.items(), "Combine Light with Time: "):
        # Add time stamps to data frame of current sample by merging
        # The time data frame for the current sample, which is expected
        # to match the light data (based on index).
        time_df = time_data.get(sample_id)
        if time_df is None:
            msg = "{}: could not find timestamp data for sample".format(sample_id)
            if error_on_missing_timestamps:
                raise ValueError(msg)
            # Ignore, if missing data shouldn't cancel the whole process.
            logging.warning(msg)
            continue

        n_timestamps = len(time_df)
        n_lightdata = len(lm_df)

        # The timestamp and light recordings are done by different systems.
        # This can cause the existence of additional time points/ or missing time points in a
        # dataset, which will be filtered out in the merge operation below.
        if n_lightdata != n_timestamps:
            msg = "{}: time data ({} entries) doesn't match light data ({} entries)".format(
                sample_id, n_timestamps, n_lightdata
            )
            if error_on_time_light_mismatch:
                raise ValueError(msg)
            logging.warning(msg)
            diffs = np.diff(time_df["time"])[10:]  # from 10th row onwards
            all_diffs = np.concatenate((all_diffs, diffs))
            diffs_avg = diffs.mean(axis=0)
            # diff between timedata and lightdata
            missing_data = len(lm_df) - len(time_df)

            # add 'diffs_avg' to fill in missing_timedata
            if missing_data > 0:
                last_valid_index = len(time_df) - 1
                last_timestamp = time_df.iloc[last_valid_index]["time"]
                if pd.isna(last_timestamp):
                    raise ValueError(
                        "Unexpected last valid timestamp for sample {} at index {}".format(
                            sample_id, last_valid_index
                        )
                    )
                for _ in range(missing_data):
                    last_valid_index += 1
                    time_df.loc[last_valid_index] = (
                        time_df.iloc[last_valid_index - 1]["time"] + diffs_avg
                    )
            elif missing_data < 0:
                drop_start = len(time_df) + missing_data
                drop_end = len(time_df)
                time_df.drop(list(range(drop_start, drop_end)))

        # Merge timedata into light data
        # Use an 'inner' join/merge to exclude time points that don't have matching light data.
        new_sample_df = pd.merge(
            lm_df, time_df, left_index=True, right_index=True, how="inner"
        )

        # Store newly created data frame for sample (dictionary)
        lm_data[sample_id] = new_sample_df

    if show_diffs:
        # Max.diffs for timestamps
        # diffs defined earlier
        mx = all_diffs.max()
        logging.info("max-diff", mx)

        plt.hist(all_diffs, bins=10, alpha=0.5)
        plt.show()

    logging.info(
        "Matched {} light data sets with their respective time points".format(
            len(lm_data)
        )
    )


def merge_dataframe_list(list_of_dfs, ordered=True, **kwargs):
    combined_df = list_of_dfs[0]
    for right_df in list_of_dfs[1:]:
        if ordered:
            combined_df = pd.merge_ordered(combined_df, right_df, **kwargs)
        else:
            combined_df = pd.merge(combined_df, right_df, **kwargs)
    return combined_df


def split_window(df, pre=0, post=0):
    pre_data = df[df.index < 0.0]
    post_data = df[df.index > 2.0]
    return pre_data, post_data


def get_fold_change_df(df, pre=0, post=0):
    pre_data, post_data = split_window(df=df, pre=pre, post=post)

    pre_data_avg = pre_data.mean(axis=0)
    post_data_avg = post_data.mean(axis=0)

    # fold changeI: post/pre
    # fold_change = (post_data_avg)/(pre_data_avg)

    # fold changeII: post-pre/pre
    fold_change = ((post_data_avg) - (pre_data_avg)) / (pre_data_avg)

    fold_change_df = fold_change.to_frame("transitions")
    return fold_change_df


def extract_transition_dfs(df):

    transitions = set()
    for sample_transition in df.index:
        current = "_".join(sample_transition.split("_")[-2:])
        transitions.add(current)
    transitions = list(transitions)

    logging.info("Found transitions: {}".format(transitions))

    transition_df_list = [df[df.index.str.contains(t)] for t in transitions]

    return transitions, transition_df_list
