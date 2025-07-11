import numpy as np

def calculate_average_ibi():
    """
    Identifies burst timestamps from the graph, calculates the inter-burst intervals (IBIs),
    and computes the average IBI.
    """
    # Step 1: Visually identify the approximate center timestamp for each synchronized burst.
    # The bursts are the dense vertical bands of activity.
    burst_timestamps = [3.0, 16.5, 25.5, 36.0, 49.0]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the time difference between consecutive bursts (inter-burst intervals).
    # Using np.diff is a concise way to do this.
    intervals = np.diff(burst_timestamps)
    # The np.diff function returns a numpy array, we can convert it to a list for printing if desired
    intervals_list = list(intervals)
    print(f"Calculated inter-burst intervals (s): {intervals_list}")

    # Step 3: Calculate the average of the intervals.
    average_interval = np.mean(intervals)

    # Step 4: Display the full calculation equation and the result.
    # We construct the string for the equation to show how the average was calculated.
    sum_str = " + ".join(map(str, intervals_list))
    num_intervals = len(intervals_list)
    
    print("\n--- Calculation Details ---")
    print(f"Average = ({sum_str}) / {num_intervals}")
    print(f"Average = {sum(intervals_list)} / {num_intervals}")
    print(f"Result = {average_interval:.3f} s")

    # Format the final answer to two significant figures as requested.
    # The 'g' format specifier is used for significant figures.
    formatted_average = f"{average_interval:.2g}"
    print(f"\nThe average inter-burst interval to two significant figures is {formatted_average} s.")

# Run the function to perform the calculation and print the results.
calculate_average_ibi()