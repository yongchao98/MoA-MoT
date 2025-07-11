import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from a graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1: Identify the timestamps of the synchronized bursts from the graph.
    burst_timestamps = [3.5, 17.5, 35.5, 49.5]
    print(f"Identified burst timestamps (s): {burst_timestamps}")
    print("-" * 30)

    # Step 2: Calculate the inter-burst intervals (IBIs).
    intervals = []
    print("Calculating Inter-Burst Intervals (IBIs):")
    for i in range(len(burst_timestamps) - 1):
        interval = burst_timestamps[i+1] - burst_timestamps[i]
        intervals.append(interval)
        print(f"  Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {interval:.1f} s")
    print("-" * 30)

    # Step 3: Calculate the average IBI.
    average_interval = np.mean(intervals)
    print("Calculating the average IBI:")
    # Create the sum part of the equation string
    sum_str = " + ".join(map(str, [f"{i:.1f}" for i in intervals]))
    print(f"  Average = ({sum_str}) / {len(intervals)}")
    print(f"  Average = {sum(intervals):.1f} / {len(intervals)} = {average_interval:.3f} s")
    print("-" * 30)

    # Step 4: Round the result to two significant figures.
    # For a number like 15.333..., the first two significant figures are 1 and 5.
    # The next digit (3) is less than 5, so we round down.
    final_answer = 15
    print(f"The average inter-burst interval rounded to two significant figures is: {final_answer} s")

    return final_answer

# Run the calculation and print the final result in the specified format.
calculate_average_ibi()