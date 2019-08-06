import matplotlib.pyplot as plt

# Define a class with sample_id, cell_type, event_name and filter_pattern

# TODO: This class is very similar to the CellTransitionConfig, maybe combine them


class CellTraceConfig:
    """
    Takes a sample_id, cell_type, event_name, and optional filter_pattern
    """

    def __init__(self, sample_id, cell_type, event_name, filter_pattern=None):
        self.sample_id = sample_id
        self.cell_type = cell_type
        self.event_name = event_name
        self.filter_pattern = filter_pattern

    def get_filter_regex(self):
        """
        produces a filter_regex that matches strings starting with cell_type followed
        by filter_pattern if available
        """
        filter_regex = "^{}_".format(self.cell_type)
        if self.filter_pattern:
            filter_regex += ".*{}.*".format(self.filter_pattern)
        return filter_regex

    def get_event_start_col(self):
        return "{}_start".format(self.event_name)

    def add_event_time_points_to_plot(self, source_df, plot):
        """
        takes source dataframe and plot, adds event_name rows to plot
        """
        for idx, row in source_df.iterrows():
            plot.annotate(self.event_name, xy=(row["time"], 1))
            plt.axvline(row["time"], color="k", linestyle="-")


# CELL TRACE CONFIG FUNCTIONS


def extract_windows(
    sample_data, cell_trace_configs, left_half_window_size, right_half_window_size
):
    # extract cell columns
    filtered_data = cell_trace_config_filter(sample_data, cell_trace_configs)
    # extract event start rows
    event_start_dfs = get_event_start_rows(sample_data, cell_trace_configs)

    all_events = []
    for sample_id, cell_subset_df in filtered_data.items():
        event_df = event_start_dfs[sample_id]

        # Extract for specific time window and align several events.
        # Define timepoints pre and post an event (event_df).
        # This works for single sample or multiple samples aligned
        # Note: In cell_subset_df, time was set to index, because for the previous avg calculation
        # Add index and time = column

        # Event_df defined in pargraph before
        windows = (
            []
        )  # TODO: use list per sample_id (simplify column names and sample_id tracking)
        n_behavior = 0  # the nth behavior of this type for this sample
        for i, row in event_df.iterrows():
            n_behavior += 1
            window_start = row["time"] - left_half_window_size
            window_end = row["time"] + right_half_window_size

            # Get subset of rows between window_start and window_end
            event = cell_subset_df[
                (cell_subset_df.time >= window_start)
                & (cell_subset_df.time <= window_end)
            ]
            # Normalizing the data to align on beginning of selected
            # behavior (event_df = Zero) by substracting events in window
            # around start of event of interest from start of event interest.
            # Note: using ":" in event.loc[] will select "all rows" in our window.
            event.loc[:, "time"] = event["time"] - row["time"]

            # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
            event.rename(
                lambda x: "{}_{}_{}".format(sample_id, x, n_behavior),
                axis="columns",
                inplace=True,
            )

            # Rename time collum to time
            event.rename(columns={event.columns[0]: "time"}, inplace=True)
            all_events.append(event)  # Append a list with all event

            # Round (NR)
            # decimals = 1
            # event['time'] = event['time'].apply(lambda x: round(x, decimals))

    return all_events


def cell_trace_config_filter(sample_data, ctcs):
    """
    Extract columns matching our cell type and the optional filter pattern.
    """
    filtered_data = {}
    for ctc in ctcs:
        sample_df = sample_data.get(ctc.sample_id)
        if sample_df is None:
            raise ValueError("{}: could not find sample data".format(ctc.sample_id))
        cell_subset_df = sample_df.filter(regex=ctc.get_filter_regex())
        # set index = time, then reset index and make time a colum
        # should be the same as cell_subset_df["time"] = sample_df["time"]
        cell_subset_df.set_index(sample_df.time, inplace=True)
        cell_subset_df.reset_index(inplace=True)
        filtered_data[ctc.sample_id] = cell_subset_df
    return filtered_data


def get_event_start_rows(sample_data, ctcs):
    """
    filters forws of sample_data dfs by event start
    """
    event_start_dfs = {}
    for ctc in ctcs:
        sample_df = sample_data.get(ctc.sample_id)
        event_df = sample_df[sample_df.loc[:, ctc.get_event_start_col()]]
        for idx, row in event_df.iterrows():
            print(
                "TP of {} = {}, sample_id: {}".format(
                    ctc.event_name, row["time"], "{}".format(ctc.sample_id)
                )
            )
        event_start_dfs[ctc.sample_id] = event_df
    return event_start_dfs
