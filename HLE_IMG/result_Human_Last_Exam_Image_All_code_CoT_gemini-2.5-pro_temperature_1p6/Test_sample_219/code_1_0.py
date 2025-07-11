import numpy as np

def calculate_average_ibi():
    """
    Identifies synchronized burst timestamps from a graph,
    calculates the inter-burst intervals (IBIs),
    and computes the average IBI.
    """
    # Step 1 & 2: Identify and list the timestamps of synchronized bursts.
    # Visually inspecting the graph, dense vertical bands of firing activity
    # are prominent at approximately these times in seconds.
    burst_timestamps = [3, 17, 36, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 3: Calculate the inter-burst intervals.
    # The interval is the difference in time between consecutive bursts.
    intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {intervals.tolist()}")

    # Step 4: Calculate the average inter-burst interval.
    sum_of_intervals = np.sum(intervals)
    number_of_intervals = len(intervals)
    average_interval = np.mean(intervals)

    # Step 5: Format the final answer to two significant figures.
    # We use format specifiers to handle the rounding.
    num_sig_figs = 2
    formatted_average = f"{average_interval:.{num_sig_figs}g}"

    # Print the detailed calculation steps
    print("\n--- Calculation of Average Inter-Burst Interval ---")
    interval_sum_str = " + ".join(map(str, intervals))
    print(f"Average = ({interval_sum_str}) / {number_of_intervals}")
    print(f"Average = {sum_of_intervals} / {number_of_intervals}")
    print(f"Average = {average_interval:.4f} s")
    print(f"\nResult rounded to two significant figures: {formatted_average} s")

# Execute the function to get the result
calculate_average_ibi()