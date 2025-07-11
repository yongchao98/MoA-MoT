import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from a neuronal activity graph,
    calculates the inter-burst intervals (IBIs), and computes their average.
    """
    # Step 1: Identify the timestamps of synchronized bursts from the graph.
    # These are the approximate center times (in seconds) of the dense vertical bands of activity.
    burst_timestamps = [2.5, 17, 35, 48]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (IBIs).
    # The IBI is the time difference between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    
    print("\nCalculating Inter-Burst Intervals (IBIs):")
    for i in range(len(inter_burst_intervals)):
        print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {inter_burst_intervals[i]} s")

    # Step 3: Calculate the average IBI.
    average_ibi = np.mean(inter_burst_intervals)
    
    # Step 4: Print the final calculation and the result rounded to two significant figures.
    intervals_sum_str = " + ".join(map(str, inter_burst_intervals))
    print(f"\nAverage Inter-Burst Interval = ({intervals_sum_str}) / {len(inter_burst_intervals)}")
    print(f"Average Inter-Burst Interval = {np.sum(inter_burst_intervals)} / {len(inter_burst_intervals)} = {average_ibi:.4f} s")
    
    # Rounding to two significant figures.
    # Using format specifier 'g' for significant figures.
    final_answer = f"{average_ibi:.2g}"
    print(f"\nThe average inter-burst interval rounded to two significant figures is: {final_answer} s")

calculate_average_ibi()