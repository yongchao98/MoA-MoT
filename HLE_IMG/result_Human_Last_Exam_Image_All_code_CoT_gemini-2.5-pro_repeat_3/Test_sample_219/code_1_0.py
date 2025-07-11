import numpy as np

def calculate_average_ibi():
    """
    Identifies synchronized bursts from a neuronal activity graph,
    calculates the inter-burst intervals, and computes the average.
    """
    # Step 1: Identify the timestamps of the synchronized bursts.
    # By visual inspection of the graph, the dense vertical bands (bursts)
    # are centered at approximately the following times in seconds.
    burst_timestamps = [3, 17, 35, 48]
    print(f"Identified synchronized burst timestamps: {burst_timestamps} s")

    # Step 2: Calculate the inter-burst intervals (IBIs).
    # The IBI is the difference in time between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    print("Calculating the inter-burst intervals:")
    for i in range(len(inter_burst_intervals)):
        print(f"  Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {inter_burst_intervals[i]} s")

    # Step 3: Calculate the average inter-burst interval.
    average_ibi = np.mean(inter_burst_intervals)
    
    # Create the string for the equation
    interval_sum_str = " + ".join(map(str, inter_burst_intervals))
    calculation_str = f"({interval_sum_str}) / {len(inter_burst_intervals)}"
    
    print("\nCalculating the average inter-burst interval:")
    print(f"  Average = {calculation_str} = {average_ibi:.4f} s")

    # Step 4: Round the result to two significant figures.
    # The result 15.0 already has two significant figures.
    # We will format it to ensure it is displayed correctly.
    formatted_average_ibi = f"{average_ibi:.2g}"
    
    print(f"\nThe average inter-burst interval rounded to two significant figures is {formatted_average_ibi} s.")

calculate_average_ibi()