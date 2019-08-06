from tqdm import tqdm
from functions import AVAILABLE_BEHAVIORS


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
        # Look for start of second behavior and remember its time.
        if row[second_event_start_col] and not row[second_event_overlap_col]:
            # print("{} starts at {}".format(second_event, row["time"]))
            second_event_start_time = row["time"]
        if row[first_event_end_col]:
            # print("{} ends at {}".format(first_event, row["time"]))
            first_event_end_time = row["time"]
        if row[first_event_start_col]:
            # print("{} starts at {}".format(first_event, row["time"]))
            first_event_start_time = row["time"]
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
                    # print("{} ended at {}, but then found {} at {}".format(first_event, first_event_end_time, column, row["time"]))
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
            continue
        if first_event_start_time > first_event_end_time:
            continue

        # Test if first_event_start_time = second_event_start_time
        if abs(first_event_start_time - second_event_start_time) < 0.00001:
            print(
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
            print("warning: end time not found for second event")

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
                continue
            if (
                second_event_duration is not None
                and second_event_end_time - second_event_start_time
                < second_event_duration
            ):
                continue

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
        first_event_start_time = None
        first_event_end_time = None
        second_event_start_time = None
        second_event_end_time = None

    return results


def find_behavior_next(
    sample_id,
    sample_df,
    first_event,
    second_event,
    max_delay=0,
    max_ignored_quiet_time=float("inf"),
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

    print("finding same behaviors only")

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

        if first_event_start_time:
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
                }
            )

        # Reset behavior tracking variables to find new pattern match.
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


def extract_transitions(sample_data, behavior_transitions):
    found_transitions = []
    for bt in tqdm(behavior_transitions):
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
            first_event_duration=None,
            second_event_duration=None,
        )  # For 'quiet' change *_event_duration. Defaul = None.

        if transitions:
            found_transitions.append(transitions)
    return found_transitions
