import logging
from functions import get_sample_data
from classes import CellTraceConfig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker

from logging_functions import count_cell_filterpatterns
from plotting_functions import plot_all_events
from classes.cell_trace_config import extract_windows, cell_trace_config_filter
from plotting_functions import plot_windowed_events

PLOT_EVERY_SAMPLE = True

# Logging info:
# Logging has 4 levels (DEBUG, INFO, WARNING, ERROR)
# appends to log file rather than overwritting it
logging.basicConfig(level=logging.INFO, filename="log.log")

# TODO: Cleanup comments in get_sample_data
sample_data = get_sample_data()


# %%
########################################################################################################################
# DATA ANALYSIS
########################################################################################################################

cell_trace_configs = [
    CellTraceConfig("17-09-01L3-cl", "A00c", "stim", "mid"),
    CellTraceConfig("17-08-26L6-cl", "A00c", "stim", "mid"),
    CellTraceConfig("17-08-28L3-cl", "A00c", "stim"),
]

# Load all samples with specific filter pattern
"""
cell_trace_configs = [
    CellTraceConfig(name,'A00c', 'fw') for name in sample_data]
"""


# %%

# TODO: only take all samples into account if event is quiet (maybe make event optional)

count_cell_filterpatterns(sample_data, "basin", "A4L")

all_events = extract_windows(
    sample_data, cell_trace_configs, left_half_window_size=10, right_half_window_size=50
)

logging.info("found {} events".format(len(all_events)))
filtered_data = cell_trace_config_filter(sample_data, cell_trace_configs)

if PLOT_EVERY_SAMPLE:
    plot_all_events(sample_data, cell_trace_configs)

plot_windowed_events(all_events)
