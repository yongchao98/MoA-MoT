import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from the graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1 & 2: Identify and list the approximate timestamps of synchronized bursts (in seconds)
    # These are estimated from the center of the dense vertical firing patterns in the graph.
    burst_timestamps = [2.5, 11.5, 18.5, 26.5, 36.5, 48.5]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 3: Calculate the inter-burst intervals (IBIs)
    # The IBI is the difference in time between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")

    # Step 4: Calculate the average IBI
    average_ibi = np.mean(inter_burst_intervals)
    
    # Step 5: Print the final calculation and result
    # Create the equation string
    intervals_sum_str = " + ".join(map(str, inter_burst_intervals))
    equation = f"({intervals_sum_str}) / {len(inter_burst_intervals)}"
    
    print(f"\nCalculation for the average inter-burst interval:")
    print(f"Average = {equation} = {average_ibi:.2f} s")
    
    # Format the final answer to two significant figures
    final_answer = f"{average_ibi:.1f}"
    print(f"\nThe average inter-burst interval to two significant figures is {final_answer} s.")

calculate_average_ibi()