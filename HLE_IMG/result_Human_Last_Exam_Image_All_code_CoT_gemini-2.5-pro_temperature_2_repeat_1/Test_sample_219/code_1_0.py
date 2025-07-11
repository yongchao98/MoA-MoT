import numpy as np

def calculate_average_ibi():
    """
    Identifies burst timestamps from the graph, calculates the inter-burst intervals,
    and computes the average interval.
    """
    # Step 1: List the identified timestamps of the synchronized bursts (in seconds).
    burst_timestamps = [4, 17, 27, 37, 48]
    print(f"Identified synchronized burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (IBIs).
    # This is the difference in time between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {inter_burst_intervals.tolist()}")

    # Step 3: Calculate the average inter-burst interval.
    average_ibi = np.mean(inter_burst_intervals)

    # Step 4: Display the final calculation and the result.
    # The final answer must be formatted to two significant figures.
    # The calculated average is 11.0, which already has two significant figures when written as 11.
    
    # Building the equation string for clear output
    interval_sum_str = " + ".join(map(str, inter_burst_intervals))
    num_intervals = len(inter_burst_intervals)
    
    print("\nCalculating the average inter-burst interval:")
    # Print the full equation with the numbers used in the calculation
    print(f"Average = ({interval_sum_str}) / {num_intervals}")
    print(f"Average = {sum(inter_burst_intervals)} / {num_intervals}")
    
    # The result 11 has two significant figures.
    final_answer = round(average_ibi)
    print(f"The average inter-burst interval is {final_answer} s.")

calculate_average_ibi()

# The final numerical answer.
final_answer_value = 11