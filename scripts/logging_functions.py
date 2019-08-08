from classes.cell_trace_config import CellTraceConfig, cell_trace_config_filter
import logging
import numpy as np
from scipy.stats import mannwhitneyu


def count_cell_filterpatterns(sample_data, cell, filter_pattern):
    cell_trace_configs = [
        CellTraceConfig(sample_id, cell, "quiet", filter_pattern)
        for sample_id in sample_data
    ]  # if all samples should be included, you have
    # to chose 'quiet' as event

    # ALL THIS DOES IS FILTER ALL SAMPLES BY CELL_TRACE_CONFIGS
    # FILTER OPERATION
    number_basins = []
    sample_ids = []
    filtered_data = cell_trace_config_filter(sample_data, cell_trace_configs)
    for sample_id, filtered_df in filtered_data.items():
        number = len(filtered_df.columns)
        number_basins.append(number)
        sample_ids.append(sample_id)
    avg = np.mean(number_basins)
    total = np.sum(number_basins)
    logging.info(total)
    logging.info(avg)


def log_num_transitions(transitions):
    logging.info("Found transitions for {} samples".format(len(transitions)))
    logging.info(
        "Found total of {} transitions".format(
            sum([len(sample_transitions) for sample_transitions in transitions])
        )
    )


def log_overlapping_transition(transition, seconds_overlap=0):
    # Test minus valus
    if (
        (transition["second_event_start"]) - (transition["first_event_end"])
    ) < seconds_overlap:
        logging.info(
            'Overlapping Transition on "{}": first_event_end: {},  second_event_start: {}, overlap: {}'.format(
                transition["sample_id"],
                transition["first_event_end"],
                transition["second_event_start"],
                transition["first_event_end"] - transition["second_event_start"],
            )
        )


def log_unique_column_values(
    df,
    sample_id=False,
    cell=False,
    filter_pattern=False,
    n_obs=False,
    first_event=False,
    second_event=False,
):
    values = (sample_id, cell, filter_pattern, n_obs, first_event, second_event)
    string_values = (
        "sample_id",
        "cell",
        "filter_pattern",
        "n_obs",
        "first_event",
        "second_event",
    )
    indicies = tuple(i for i in range(len(values)) if values[i])

    def extract_parts(column: str, indicies):
        parts = column.split("_")
        return tuple(parts[i] for i in indicies)

    parts = set([extract_parts(column, indicies) for column in df.columns])
    logging.info(
        "Counted {} unique columns for {}".format(
            len(parts), tuple(string_values[i] for i in indicies)
        )
    )
    logging.info("parts: {}".format(parts))


def log_mann_whitney_u_test(transition_dfs, transitions, comparisons):
    transition_df_map = {t: t_df for t, t_df in zip(transitions, transition_dfs)}
    mwu_results = [
        mannwhitneyu(transition_df_map[a], transition_df_map[b]) for a, b in comparisons
    ]
    info = "\n".join(
        [
            "Mann-Whitney u test comparing {} to {}: {}".format(a, b, result)
            for (a, b), result in zip(comparisons, mwu_results)
        ]
    )

    logging.info(info)
