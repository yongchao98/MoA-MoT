import numpy as np

def calculate_average_inter_burst_interval():
    """
    Identifies timestamps of synchronized bursts from a neuronal activity graph,
    calculates the inter-burst intervals, and computes the average interval.
    """
    # Step 1: Identify the timestamps of the synchronized bursts from the graph.
    # These are the time points with dense, vertical columns of firing activity.
    burst_timestamps = [7, 17, 27, 37, 47]
    print(f"Identified burst timestamps (s): {burst_timestamps}")
    print("-" * 30)

    # Step 2: Calculate the inter-burst intervals (the time between consecutive bursts).
    intervals = []
    print("Calculating inter-burst intervals:")
    for i in range(len(burst_timestamps) - 1):
        interval = burst_timestamps[i+1] - burst_timestamps[i]
        intervals.append(interval)
        print(f"Interval {i+1}: {burst_timestamps[i+1]}s - {burst_timestamps[i]}s = {interval}s")
    print("-" * 30)

    # Step 3: Calculate the average of the intervals.
    average_interval = np.mean(intervals)

    # Step 4: Display the final calculation and the result.
    # Construct the equation string for the final output.
    equation_str = f"({ ' + '.join(map(str, intervals)) }) / {len(intervals)}"
    print("Calculating the average inter-burst interval:")
    print(f"Average = {equation_str} = {average_interval:.2f}s")
    
    # Round the final answer to two significant figures.
    # For a number like 10.0, this is already correctly represented.
    # We use a formatting approach that handles this appropriately.
    final_answer = f"{average_interval:.1f}" if average_interval == 10.0 else f"{average_interval:.2g}"

    print(f"\nThe average inter-burst interval is {final_answer} s.")
    
calculate_average_inter_burst_interval()