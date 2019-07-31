from tqdm import tqdm
import logging

def combine_lm_behavior_data(
    lm_data,
    behavior_data,
    available_behaviors,
    error_on_missing_behaviors=True,
    error_on_invalid_behavior_range=True,
    error_on_time_behavior_mismatch=True,
):
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
                logging.debug("sample_id: {}, starts at {}, ends at {}".format(sample_id, start, end))

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

    return lm_df

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def combine_lm_time_data(lm_data, time_data, error_on_missing_timestamps=True, error_on_time_light_mismatch=True, show_diffs=False):
    sample_data = {}
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
        sample_data[sample_id] = new_sample_df

    if show_diffs:
        # Max.diffs for timestamps
        # diffs defined earlier
        mx = all_diffs.max()
        logging.info("max-diff", mx)

        plt.hist(all_diffs, bins=10, alpha=0.5)
        plt.show()