# Put '' [empty string] if you dont want any cell type
from tqdm import tqdm
from functions import AVAILABLE_BEHAVIORS
import logging


class CellTransConfig:
    """
    Define a class with sample_id, cell_type, event_time and filter_pattern (for behavioral_transitions)
    """

    def __init__(
        self,
        sample_id,
        cell_type,
        event_time,
        filter_pattern=None,
        first_event=None,
        second_event=None,
    ):
        self.sample_id = sample_id
        self.cell_type = cell_type
        self.event_time = event_time
        self.filter_pattern = filter_pattern
        self.first_event = first_event
        self.second_event = second_event

    def get_filter_regex(self):
        if self.cell_type is None:
            cell_str = r"[a-zA-Z0-9]+"
        else:
            cell_str = self.cell_type

        filter_regex = "^{}_".format(cell_str)
        if self.filter_pattern:
            filter_regex += ".*{}.*".format(self.filter_pattern)
        return filter_regex


def extract_windows(
    sample_data,
    cell_transition_configs,
    left_half_window_size=18.4,
    right_half_window_size=42.4,
    cell_pattern_filter = True,
):
    # trans_df defined in pargraph before
    all_Ptrans_events = []
    n_behavior_per_sample = {}

    # TODO: Split filter data and extract windows
    for ctc in tqdm(cell_transition_configs, "Extracting windows"):
        sample_df = sample_data.get(ctc.sample_id)
        n_behavior = n_behavior_per_sample.get(ctc.sample_id, 1)
        n_behavior_per_sample.setdefault(ctc.sample_id, 1)
        if sample_df is None:
            raise ValueError("{}: could not find sample data".format(ctc.sample_id))
            continue

        # TODO make regex optional and filter out behavior

        # TODO: make extracting behavior / cell columns easy and accessable anywhere
        all_columns = set(sample_df.columns)
        behavior_columns = set(
            column
            for column in all_columns
            if any(behavior in column for behavior in AVAILABLE_BEHAVIORS + ("quiet",))
        )
        cell_df = sample_df.filter(items=all_columns - behavior_columns)
        if cell_pattern_filter:
            cell_subset_df = cell_df.filter(regex=ctc.get_filter_regex())
        else:
            cell_subset_df = cell_df
        cell_subset_df.set_index(sample_df.time, inplace=True)
        cell_subset_df.reset_index(inplace=True)

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
    return all_Ptrans_events
