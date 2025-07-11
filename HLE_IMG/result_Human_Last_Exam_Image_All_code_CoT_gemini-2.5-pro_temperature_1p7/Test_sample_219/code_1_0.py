import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from the provided graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1: Identify the timestamps of the major synchronized bursts from the graph.
    # These are the dense vertical bands where many electrodes fire at once.
    # From visual inspection, the bursts occur at approximately these times (in seconds).
    burst_timestamps = [3, 14, 19, 25, 35, 46]
    print(f"Identified burst timestamps (s): {burst_timestamps}\n")

    # Step 2: Calculate the Inter-Burst Intervals (IBIs)
    # The IBI is the time difference between consecutive bursts.
    intervals = np.diff(burst_timestamps)
    print("Calculating the inter-burst intervals (IBIs):")
    for i in range(len(intervals)):
        t1 = burst_timestamps[i]
        t2 = burst_timestamps[i+1]
        print(f"  Interval {i+1}: {t2} s - {t1} s = {intervals[i]} s")
    print(f"\nCalculated IBIs (s): {list(intervals)}\n")

    # Step 3: Calculate the average IBI
    average_interval = np.mean(intervals)
    
    # Step 4: Display the final calculation and result
    # We construct the string for the final equation to show all numbers.
    interval_sum_str = " + ".join(map(str, intervals))
    num_intervals = len(intervals)
    total_sum = sum(intervals)

    print("Final average IBI calculation:")
    print(f"  Average = ({interval_sum_str}) / {num_intervals}")
    print(f"  Average = {total_sum} / {num_intervals}")
    
    # Format the final result to two significant figures.
    formatted_average = "{:.2g}".format(average_interval)
    print(f"\nThe average inter-burst interval is {average_interval:.2f} s.")
    print(f"Rounded to two significant figures, the average inter-burst interval is {formatted_average} s.")

calculate_average_ibi()