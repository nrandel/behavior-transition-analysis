from classes.cell_trace_config import CellTraceConfig, cell_trace_config_filter
import logging
import numpy as np


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
                transition["first_event_end"] - transition["second_event_start"]
            )
        )

