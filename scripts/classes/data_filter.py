# For example to average not only over all A00cs but all A00c_midL.

class DataFilter():
    """
    Define a class for filtering after behavioral_transitions for either only cell_type or filter_pattern or both.
    """
    def __init__(self, cell='.*', pattern='.*'):
        # cell and pattern are both optional since default '.*' is provided
        self.cell = cell
        self.pattern = pattern
        
    def get_cell_filter_regex(self):
        filter_regex = '.*_{}_.*_.*'.format(self.cell)
        return filter_regex
    
    def get_pattern_filter_regex(self):
        filter_regex = '.*_.*_{}_.*'.format(self.pattern)
        return filter_regex
    
    def get_cellpattern_filter_regex(self):
        filter_regex = '.*_{}_{}_.*'.format(self.cell, self.pattern)
        return filter_regex
        
    def __str__(self):
        return "{}_{}".format(self.cell, self.pattern)  