class PostBehaviorTransition:
    # TODO: maybe use max_ignored_quite_time
    def __init__(self, sample_id, event, post_event, max_delay=0):
        self.sample_id = sample_id
        self.post_event = post_event
        self.event = event
        self.max_delay = max_delay

class SamePairBehaviorTransition:
    
    def __init__(self, sample_id, pre_event, event, max_delay=0, max_ignored_quiet_time=float("inf")):
        self.sample_id = sample_id
        self.pre_event = pre_event
        self.event = event
        self.max_delay = max_delay
        self.max_ignored_quiet_time = max_ignored_quiet_time