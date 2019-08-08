import pickle

from plotting_functions import plot_transitions, plot_fold_change
from logging_functions import log_mann_whitney_u_test
from functions import get_fold_change_df, extract_transition_dfs, merge_dataframe_list

from classes.transition_type import TransitionType

import logging

logging.basicConfig(level=logging.INFO)

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

all_Ptrans_df = pickle.load(open("all_ptrans_df.obj", "rb"))

fold_change_df = get_fold_change_df(df=all_Ptrans_df, pre=0, post=2)
transitions, transition_dfs = extract_transition_dfs(fold_change_df)
all_fold_change_df = merge_dataframe_list(
    transition_dfs, ordered=False, left_index=True, right_index=True, how="outer"
)


log_mann_whitney_u_test(transition_dfs, transitions, MANN_WHITNEY_U_TEST_COMPARISONS)
