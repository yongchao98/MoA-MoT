import numpy as np

def calculate_average_interval():
    """
    Identifies timestamps of synchronized bursts from a graph, calculates the 
    inter-burst intervals, and computes the average interval.
    """
    # Step 1 & 2: Identify and list the approximate timestamps of the major synchronized bursts from the graph.
    # Bursts are visible as dense vertical columns of firing activity.
    burst_timestamps = [2, 11, 18, 27, 36, 47]
    print(f"Identified burst timestamps (s): {burst_timestamps}\n")

    # Step 3: Calculate the inter-burst intervals (time difference between consecutive bursts).
    intervals = np.diff(burst_timestamps)
    print("Calculating inter-burst intervals:")
    for i in range(len(intervals)):
        print(f"  Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {intervals[i]} s")
    print(f"\nCalculated inter-burst intervals (s): {list(intervals)}\n")
    
    # Step 4: Calculate the average of these intervals.
    sum_of_intervals = np.sum(intervals)
    number_of_intervals = len(intervals)
    average_interval = np.mean(intervals)

    # Step 5: Display the full calculation and the final answer.
    # The equation string is built to show all numbers used in the calculation.
    equation_str = f"({ ' + '.join(map(str, intervals)) }) / {number_of_intervals}"
    print(f"Average inter-burst interval calculation: {equation_str} = {average_interval:.2f} s")

    # Round the final answer to two significant figures. 
    # For a number like 9.0, this is equivalent to one decimal place.
    final_answer = float(f"{average_interval:.1f}")
    print(f"\nThe average inter-burst interval rounded to two significant figures is {final_answer} s.")

calculate_average_interval()