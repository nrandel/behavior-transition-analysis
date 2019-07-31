# Put '' [empty string] if you dont want any cell type

class CellTransConfig:
    """
    Define a class with sample_id, cell_type, event_time and filter_pattern (for behavioral_transitions)
    """

    def __init__(self, sample_id, cell_type, event_time, filter_pattern=None, first_event=None, second_event=None):
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
            
        filter_regex = '^{}_'.format(cell_str)
        if self.filter_pattern:
            filter_regex += '.*{}.*'.format(self.filter_pattern)
        return filter_regex