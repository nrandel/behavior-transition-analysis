from tqdm import tqdm
from functions import AVAILABLE_BEHAVIORS
import logging


# TODO: COMBINE `find_behavior_before` with find_same_behavior_transitions`
# TODO: THESE SHOULD BE TESTED


def find_behavior_before(
    sample_id,
    sample_df,
    first_event,
    second_event,
    max_delay=0,
    first_event_duration=None,
    second_event_duration=None,
):
    """For the data frame of a single sample <df>, find all behaviors
    of type <first_event> that is followed by the event <second_event>,
    separated by <max_delay> time. The end of <second_event> is expected
    to happen strictly after the end of <first_event>. The start time
    of <second_event> however can overlap with the end time of <first_event>.
    In this case, the time difference is negative, and still smaller than
    <max_delay>. The start time of <second_event> can be before, at or after the
    end of <first_event>.
    """
    results = []
    first_event_start_col = "{}_start".format(first_event)
    first_event_end_col = "{}_end".format(first_event)
    second_event_start_col = "{}_start".format(second_event)
    second_event_end_col = "{}_end".format(second_event)
    second_event_overlap_col = "{}_overlap".format(second_event)

    first_event_start_time = None
    first_event_end_time = None
    second_event_start_time = None
    second_event_end_time = None

    for i, row in sample_df.iterrows():
        if row[first_event_start_col]:
            first_event_start_time = row["time"]
            logging.warning(
                "{} starts at {}".format(first_event, first_event_start_time)
            )
        if row[first_event_end_col] and first_event_start_time is not None:
            first_event_end_time = row["time"]
            logging.warning(
                "{} ends at {}".format(first_event, first_event_end_time)
            )
        # Look for start of second behavior and remember its time.
        if row[second_event_start_col] and first_event_start_time is not None:  # and not row[second_event_overlap_col]:
            second_event_start_time = row["time"]
            logging.warning("{} starts at {}".format(second_event, row["time"]))
        for column in sample_df.columns:
            if (
                first_event_start_time is not None
                and column.endswith("_start")
                and "quiet" not in column
                and column != first_event_start_col
                and column != second_event_start_col
                and first_event not in column
                and second_event not in column
            ):
                if row[column]:
                    logging.warning(
                        "{} ended at {}, but then found {} at {}".format(
                            first_event, first_event_end_time, column, row["time"]
                        )
                    )
                    first_event_start_time = None
                    first_event_end_time = None
                    second_event_start_time = None
                    second_event_end_time = None

        # As long as we haven't collected all needed time points,
        # keep on searching.
        if None in (
            first_event_start_time,
            first_event_end_time,
            second_event_start_time,
        ):
            continue

        # Define rules for event_start_time and event_end_time
        if first_event_start_time > second_event_start_time:
            logging.warning("first event start time is not less than second event start time")
            continue
        if first_event_start_time > first_event_end_time:
            logging.warning("first event start time is not less than first event end time")
            continue

        # Test if first_event_start_time = second_event_start_time
        if abs(first_event_start_time - second_event_start_time) < 0.00001:
            logging.info(
                "{}: start time (first) event {} and start time of (second) event {} are the same: {}".format(
                    sample_id, first_event, second_event, first_event_start_time
                )
            )

        if second_event_end_time is None:
            for j, row in sample_df.loc[i:, :].iterrows():
                if row[second_event_end_col]:
                    second_event_end_time = row["time"]
                    break
        if second_event_end_time is None:
            logging.warning("end time not found for second event")

        # Test time between first event end and second event start. If it
        # is smaller than <max_delay>, store start of second event as result.
        # The first event end time being greater than the second event start
        # time, is explicitly allowed.
        # implement event duration (for quiet)
        if (second_event_start_time - first_event_end_time) <= max_delay:
            if (
                first_event_duration is not None
                and first_event_end_time - first_event_start_time < first_event_duration
            ):
                logging.warning("first event duration too short")
                continue
            if (
                second_event_duration is not None
                and second_event_end_time - second_event_start_time
                < second_event_duration
            ):
                logging.warning("second event duration too short")
                continue
            logging.warning(
                "found transition {}, first event {}-{}, second event: {}-{}".format(
                    "-".join([first_event, second_event]),
                    first_event_start_time,
                    first_event_end_time,
                    second_event_start_time,
                    second_event_end_time,
                )
            )

            results.append(
                {
                    "sample_id": sample_id,
                    "first_event_start": first_event_start_time,
                    "first_event_end": first_event_end_time,
                    "second_event_start": second_event_start_time,
                    "second_event_end": second_event_end_time,
                    "first_event": first_event,
                    "second_event": second_event,
                }
            )
        else:
            logging.warning("delay too long")

        # Reset behavior tracking variables to find new pattern match.
        first_event_start_time = None
        first_event_end_time = None
        second_event_start_time = None
        second_event_end_time = None

    return results


def find_same_behavior_transitions(
    sample_id,
    sample_df,
    first_event,
    second_event,
    max_delay=0,
    max_ignored_quiet_time=float("inf"),
    first_event_duration=None,
    second_event_duration=None,
):
    """For the data frame of a single sample <df>, find all behaviors
    of type <first_event> that will be followed by the same event <second_event>,
    separated by <max_delay> time. The start of <first_event> is expected
    to happen strictly before the start of <second_event>. The end time
    of <first_event> however can overlap with the start time of <second_event>.
    In this case, the time difference is negative, and still smaller than
    <max_delay>. The end time of <first_event> can be before, at or after the
    end of <second_event>.

    If <first_event> and <second_event> are the same type of behavior,
    overlaps have to be taken into account differently to match start and end times
    to the correct event. During iteration for one loop, we have to exclude the
    fact that the first_event == second_event.

    Assumption: No overlapp
    Not in data sets identified
    """

    assert first_event == second_event, "This function only works for same behaviors"

    first_event_duration = (
        first_event_duration if first_event_duration is not None else 0
    )
    second_event_duration = (
        second_event_duration if second_event_duration is not None else 0
    )

    results = []
    first_event_start_col = "{}_start".format(first_event)
    first_event_end_col = "{}_end".format(first_event)
    first_event_overlap_col = "{}_overlap".format(first_event)
    second_event_start_col = "{}_start".format(second_event)
    second_event_end_col = "{}_end".format(second_event)
    second_event_overlap_col = "{}_overlap".format(second_event)

    # Construct a list without quiet and passed in behaviors
    other_behavior_cols = get_unrelated_nonquiet_behaviors((first_event, second_event))
    other_behavior_start_cols = ["{}_start".format(c) for c in other_behavior_cols]

    first_event_start_time = None
    first_event_end_time = None
    quiet_event_start_time = None
    quiet_event_end_time = None
    second_event_start_time = None
    second_event_end_time = None

    # Check for overlap between the same behaviors (print index, where 'True') and use
    # it as a check that there is not this error in the behavior data
    # print(sample_id, sample_df.index[sample_df['bw_overlap']].tolist())

    # Note: The overlap statement was removed. This part has to be
    # checked if overlapping events are found in the data
    for i, row in sample_df.iterrows():
        # Look for start of first behavior and remember its time.
        if row[first_event_start_col] and first_event_start_time is None:
            first_event_start_time = row["time"]
        if row[first_event_end_col] and first_event_end_time is None:
            first_event_end_time = row["time"]
        if row[second_event_start_col] and first_event_start_time is not None:
            second_event_start_time = row["time"]
        if row[second_event_end_col] and first_event_end_time is not None:
            second_event_end_time = row["time"]

        if first_event_start_time is not None:
            if row["quiet_start"]:
                quiet_event_start_time = row["time"]
            elif quiet_event_start_time is not None and row["quiet_end"]:
                quiet_event_end_time = row["time"]

            reset = False

            # If we found a quiet start and end time *after* the start of
            # the first event has been found, reset pattern search if the
            # quiet period is larger than max_ignored_quiet_time.
            if quiet_event_start_time is not None and quiet_event_end_time is not None:
                quiet_duration = quiet_event_end_time - quiet_event_start_time
                if quiet_duration > max_ignored_quiet_time:
                    reset = True

            # If any non-quiet and non-passed-in behavior is start
            if any(row[c] for c in other_behavior_start_cols):
                reset = True

            if reset:
                first_event_start_time = None
                first_event_end_time = None
                quiet_event_start_time = None
                quiet_event_end_time = None
                second_event_start_time = None
                second_event_end_time = None

        # As long as we haven't collected all needed time points,
        # keep on searching.
        if None in (
            first_event_start_time,
            first_event_end_time,
            second_event_start_time,
            second_event_end_time,
        ):
            continue

        if first_event_start_time == second_event_start_time:
            continue
        if first_event_end_time == second_event_end_time:
            continue
        if first_event_start_time > first_event_end_time:
            continue
        if first_event_start_time > second_event_start_time:
            continue
        if first_event_end_time - first_event_start_time < first_event_duration:
            continue
        if second_event_end_time - second_event_start_time < second_event_duration:
            continue

        # Test time between first event end and second event start. If it
        # is smaller than <max_delay>, store start of second event as result.
        # The first event end time being greater than the second event start
        # time, is explicitly allowed. During iteration the first_event == second_event.
        if second_event_start_time == first_event_end_time:  # NR
            continue
        if (second_event_start_time - first_event_end_time) <= max_delay:
            results.append(
                {
                    "sample_id": sample_id,
                    "first_event_start": first_event_start_time,
                    "first_event_end": first_event_end_time,
                    "second_event_start": second_event_start_time,
                    "second_event_end": second_event_end_time,
                    "first_event": first_event,
                    "second_event": second_event,
                }
            )

        # Reset behavior tracking variables to find new pattern match.
        # THIS HANDLES SAME TRANSITION BEHAVIOR
        first_event_start_time = second_event_start_time
        first_event_end_time = second_event_end_time
        second_event_start_time = None
        second_event_end_time = None
        quiet_event_start_time = None
        quiet_event_end_time = None

    return results


def get_unrelated_nonquiet_behaviors(excluding=[]):
    if "turn" in excluding:
        turn_behaviors = ["left turn", "right turn"]
        turn_behaviors.extend(excluding)
        excluding = turn_behaviors
    return [ab for ab in AVAILABLE_BEHAVIORS if ab not in excluding]


def extract_transitions(
    sample_data,
    behavior_transitions,
    first_event_duration=None,
    second_event_duration=None,
):
    found_transitions = []
    for bt in tqdm(behavior_transitions, "Extracting Transitions: "):
        sample_df = sample_data.get(bt.sample_id)
        # skip samples without behavior data
        if not any(["bw" in column for column in sample_df.columns]):
            continue
        if sample_df is None:
            raise ValueError("No data found for sample {}".format(bt.sample_id))
        transitions = find_behavior_before(
            bt.sample_id,
            sample_df,
            bt.event,
            bt.post_event,
            bt.max_delay,
            first_event_duration=first_event_duration,
            second_event_duration=second_event_duration,
        )  # For 'quiet' change *_event_duration. Defaul = None.

        if transitions:
            found_transitions.append(transitions)
    return found_transitions


def extract_same_transitions(
    sample_data,
    behavior_transitions,
    first_event_duration=None,
    second_event_duration=None,
):
    found_transitions = []
    for bt in tqdm(behavior_transitions, "Extracting Same Transitions: "):
        sample_df = sample_data.get(bt.sample_id)
        if not any(["bw" in column for column in sample_df.columns]):
            continue
        if sample_df is None:
            raise ValueError("No data found for sample {}".format(bt.sample_id))
        transitions = find_same_behavior_transitions(
            bt.sample_id,
            sample_df,
            bt.pre_event,
            bt.event,
            bt.max_delay,
            bt.max_ignored_quiet_time,
            first_event_duration=first_event_duration,
            second_event_duration=second_event_duration,
        )

        if transitions:
            found_transitions.append(transitions)
    return found_transitions


def extract_transition_triples(
    sample_data,
    behavior_transition_triples,
    first_trans_duration=None,
    second_trans_duration=None,
    third_trans_duration=None,
):
    found_transitions = []
    for first_bt, second_bt in tqdm(behavior_transition_triples):
        transitions = []
        # This function should assume working input transition tripples, so the assert is redundent
        assert first_bt.sample_id == second_bt.sample_id, "{} does not match {}".format(
            first_bt.sample_id, second_bt.sample_id
        )
        sample_df = sample_data.get(first_bt.sample_id)
        if sample_df is None:
            raise ValueError("No data found for sample {}".format(first_bt.sample_id))
        if not any(["bw" in column for column in sample_df.columns]):
            continue
        first_transitions = find_behavior_before(
            first_bt.sample_id,
            sample_df,
            first_bt.event,
            first_bt.post_event,
            first_bt.max_delay,
            first_event_duration=first_trans_duration,
            second_event_duration=second_trans_duration,
        )
        second_transitions = find_behavior_before(
            second_bt.sample_id,
            sample_df,
            second_bt.event,
            second_bt.post_event,
            second_bt.max_delay,
            first_event_duration=second_trans_duration,
            second_event_duration=third_trans_duration,
        )
        # print("{} transitions from {} to {}".format(len(first_transitions), first_bt.event, first_bt.post_event))
        # print("{} transitions from {} to {}".format(len(second_transitions), second_bt.event, second_bt.post_event))

        for ft in first_transitions:
            for st in second_transitions:
                if abs(ft["second_event_start"] - st["first_event_start"]) < 0.00001:
                    transitions.append(
                        {
                            "sample_id": ft["sample_id"],
                            "first_event_start": ft["first_event_start"],
                            "first_event_end": ft["first_event_end"],
                            "first_event": ft["first_event"],
                            "second_event_start": st["first_event_start"],
                            "second_event_end": st["first_event_end"],
                            "second_event": st["first_event"],
                            "third_event_start": st["second_event_start"],
                            "third_event_end": st["second_event_end"],
                            "third_event": st["second_event"],
                        }
                    )
        if transitions:
            print("{} transition triples found".format(len(transitions)))
            found_transitions.append(transitions)

    return found_transitions
