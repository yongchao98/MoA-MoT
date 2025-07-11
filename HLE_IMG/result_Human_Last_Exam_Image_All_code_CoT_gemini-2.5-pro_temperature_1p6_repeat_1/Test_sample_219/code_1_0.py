import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from a graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1: Identify the timestamps of the synchronized bursts from the graph.
    # These are visually estimated from the densest vertical bands of activity.
    burst_timestamps = [5, 17, 36, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the intervals between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")
    
    # Step 3: Calculate the average of the inter-burst intervals.
    average_ibi = np.mean(inter_burst_intervals)
    
    # Print the calculation steps for the average
    intervals_str = ' + '.join(map(str, inter_burst_intervals))
    print(f"Calculation for average IBI: ({intervals_str}) / {len(inter_burst_intervals)}")
    print(f"Unrounded average IBI: {average_ibi:.4f} s")

    # Step 4: Round the average to two significant figures.
    # The format specifier 'g' is used for general format, which works well for significant figures.
    rounded_average_ibi = float(f"{average_ibi:.2g}")
    print(f"\nAverage inter-burst interval rounded to two significant figures: {rounded_average_ibi} s")

calculate_average_ibi()