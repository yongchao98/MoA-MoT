import numpy as np

def calculate_burst_interval():
    """
    Identifies burst timestamps from a visual analysis of the provided graph,
    calculates the inter-burst intervals, and computes the average interval.
    """
    # Step 1: Identify timestamps of synchronized bursts from the graph.
    # These are visually estimated from the centers of the dense vertical bands.
    burst_timestamps = [3.5, 12.5, 24.5, 35.5, 47.5]
    print(f"Step 1: Identified synchronized burst timestamps (in seconds): {burst_timestamps}\n")

    # Step 2: Calculate the inter-burst intervals (the time between consecutive bursts).
    intervals = np.diff(burst_timestamps)
    print("Step 2: Calculate the inter-burst intervals.")
    for i, interval in enumerate(intervals):
        print(f"  Interval {i+1}: {burst_timestamps[i+1]}s - {burst_timestamps[i]}s = {interval:.1f}s")
    print("\nCalculated intervals:", [round(i, 1) for i in intervals], "s\n")

    # Step 3: Calculate the average inter-burst interval.
    average_interval = np.mean(intervals)
    print("Step 3: Calculate the average inter-burst interval.")
    
    # Constructing and printing the equation for the average calculation.
    interval_strs = [f"{i:.1f}" for i in intervals]
    equation_str = f"({ ' + '.join(interval_strs) }) / {len(intervals)}"
    print(f"  Equation: Average = {equation_str}")
    print(f"  Calculation: Average = {np.sum(intervals):.1f} / {len(intervals)} = {average_interval:.3f}s\n")

    # Step 4: Round the result to two significant figures.
    # The format specifier '.2g' rounds to 2 significant figures.
    final_answer = float(f"{average_interval:.2g}")
    print(f"Step 4: Rounding the average to two significant figures gives: {final_answer}s")

calculate_burst_interval()