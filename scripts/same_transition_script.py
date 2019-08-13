from tqdm import tqdm
import pandas as pd
import logging

from functions import get_sample_data, merge_dataframe_list

from transition_functions import extract_transition

from logging_functions import log_num_transitions

from plotting_functions import plot_transition_gaps_hist, plot_windowed_events

from classes.cell_triple_transition_config import CellTripleTransConfig, extract_windows
from classes.post_behavior_transitions import PostBehaviorTransition


sample_data = get_sample_data()

FILTER_PATTERN = "A1R"

# Open all samples
BEHAVIOR_TRANSITIONS = [
    SamePairBehaviorTransition(
        name, "turn", "turn", max_delay=10, max_ignored_quiet_time=float("inf")
    )
    for name in sample_data
]


FIRST_TRANS_DURATION = None
SECOND_TRANS_DURATION = None

USE_CELL_PATTERN_FILTER = False


# Open all samples

found_transition_triples = extract_transition_triples(
    sample_data,
    TRANSITION_TRIPLES,
    first_trans_duration=FIRST_TRANS_DURATION,
    second_trans_duration=SECOND_TRANS_DURATION,
    third_trans_duration=THIRD_TRANS_DURATION,
)

log_num_transitions(found_transition_triples)

# TODO: Probably doesn't handle  triples properly
plot_transition_gaps_hist(found_transition_triples)
# %%
# Use the predefined CellTransConfig to filter by celltype and pattern.
# The results are merged_ordered and interpolated.
cell_Ttrans_configs = []
all_Ttrans_events = []

# TODO: MAKE NEW SCRIPT FOR TRANSITION TRIPLES
# TODO: MAKE NEW CLASS FOR CellTriplesTransConfig
# TODO: CHANGE COLUMN RENAMING IN extract_windows TO USE ALL 3 EVENTS
for sample in tqdm(found_transition_triples):
    sample_ls_trans = []
    for found_transition in sample:
        logging.info(found_transition)

        # For all behavior
        sample_ls_trans.append(found_transition["second_event_start"])
        cell_Ttrans_configs.append(
            CellTripleTransConfig(
                found_transition["sample_id"],
                "A00c",
                found_transition["second_event_start"],  # event start for alignments
                FILTER_PATTERN,
                found_transition["first_event"],
                found_transition["second_event"],
                found_transition["third_event"],
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

all_Ttrans_windows = extract_windows(
    sample_data,
    cell_Ttrans_configs,
    cell_pattern_filter=USE_CELL_PATTERN_FILTER,
    left_half_window_size=left_half_window_size,
    right_half_window_size=right_half_window_size,
)
all_Ttrans_df = merge_dataframe_list(all_Ttrans_windows, ordered=True, on="time", how="outer")

all_Ttrans_df.index = all_Ttrans_df["time"]
del all_Ttrans_df["time"]

int_all_Ttrans_df = all_Ttrans_df.interpolate(
    method="index", axis=0, limit=None, inplace=False, limit_direction="both"
)

plot_windowed_events(all_Ttrans_df)
