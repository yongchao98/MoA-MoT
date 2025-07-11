import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from the provided graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1 & 2: Identify and list the timestamps of the synchronized bursts (in seconds)
    # from visual inspection of the graph.
    burst_timestamps = [3, 16, 27, 36, 48]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 3: Calculate the inter-burst intervals (IBIs)
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")

    # Step 4: Calculate the average IBI
    average_ibi = np.mean(inter_burst_intervals)

    # Create the equation string for display
    interval_sum_str = " + ".join(map(str, inter_burst_intervals))
    equation = f"({interval_sum_str}) / {len(inter_burst_intervals)}"

    print("\nCalculating the average inter-burst interval:")
    print(f"Average IBI = {equation} = {average_ibi:.4f} s")

    # Step 5: Round the result to two significant figures
    # For a number like 11.25, rounding to 2 sig figs results in 11.
    final_answer = float(f"{average_ibi:.2g}")
    print(f"\nThe average inter-burst interval rounded to two significant figures is: {final_answer} s")

    return final_answer

# Execute the function and store the final answer
result = calculate_average_ibi()

# The final answer in the required format
# print(f"\n<<<{result}>>>")