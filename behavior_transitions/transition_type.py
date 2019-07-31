# Specific after Post-transitions for multiple transition kinds, used for plotting. For multiple transition 
# events, group after transition (first, or second event) <most useful> with option to group
# after celltype, filterpattern, sample_id, observations.
class TransitionType:
    def __init__(self, sample_id=".*", cell=".*", filter_pattern=".*", n_obs=".*", first_event=".*", second_event=".*"):
        self.sample_id = sample_id
        self.cell = cell
        self.filter_pattern = filter_pattern
        self.n_obs = n_obs
        self.first_event = first_event
        self.second_event = second_event
        
        self.pattern = "{}_{}_{}_{}_{}_{}"
        
    def get_filter_regex(self, use_all=False, use_cell=False, use_sample=False, use_filter_pattern=False, use_n_obs=False, use_first_event=False, use_second_event=False):
        filter_regex = self.pattern.format(self.sample_id if use_sample or use_all else ".*",
                                          self.cell if use_cell or use_all else ".*",
                                          self.filter_pattern if use_filter_pattern or use_all else ".*",
                                          self.n_obs if use_n_obs or use_all else ".*",
                                          self.first_event if use_first_event or use_all else ".*",
                                          self.second_event if use_second_event or use_all else ".*")
        return filter_regex    