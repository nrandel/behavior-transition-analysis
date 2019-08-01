# Define a class with sample_id, cell_type, event_name and filter_pattern

# TODO: This class is very similar to the CellTransitionConfig, maybe combine them

class CellTraceConfig:
    """
    Takes a sample_id, cell_type, event_name, and optional filter_pattern
    
    """
    
    def __init__(self, sample_id, cell_type, event_name, filter_pattern=None):
        self.sample_id = sample_id
        self.cell_type = cell_type
        self.event_name = event_name
        self.filter_pattern = filter_pattern
        
    def get_filter_regex(self):
        """
        produces a filter_regex that matches strings starting with cell_type followed
        by filter_pattern if available
        """
        filter_regex = '^{}_'.format(self.cell_type)
        if self.filter_pattern:
            filter_regex += '.*{}.*'.format(self.filter_pattern)
        return filter_regex
    
    def get_event_start_col(self):
        return '{}_start'.format(self.event_name)

    def add_event_time_points_to_plot(self, source_df, plot):
        """
        takes source dataframe and plot, adds event_name rows to plot
        """
        for idx, row in source_df.iterrows():
            plot.annotate(self.event_name, xy=(row['time'], 1))
            plt.axvline(row['time'], color='k', linestyle='-')  