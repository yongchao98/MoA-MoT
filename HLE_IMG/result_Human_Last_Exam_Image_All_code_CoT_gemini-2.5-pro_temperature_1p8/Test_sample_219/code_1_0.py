import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from a neuronal activity graph,
    calculates the inter-burst intervals, and computes their average.
    """
    # Step 1: Identify timestamps of synchronized bursts from the graph
    burst_timestamps = [2.5, 16, 27, 37, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (IBIs)
    inter_burst_intervals = np.diff(burst_timestamps)
    
    # Step 3: Calculate the average inter-burst interval
    average_ibi = np.mean(inter_burst_intervals)
    
    # Step 4: Round the result to two significant figures
    # To round to n significant figures, we can use a format string.
    # The format '{value:.{sf}g}' formats the 'value' to 'sf' significant figures.
    rounded_average_ibi = float(f"{average_ibi:.2g}")

    # Print the detailed calculation
    intervals_str = " + ".join(map(str, inter_burst_intervals))
    print("\nCalculation of the average inter-burst interval:")
    print(f"({intervals_str}) / {len(inter_burst_intervals)} = {average_ibi:.4f} s")
    
    print(f"\nThe average inter-burst interval is {average_ibi:.4f} s.")
    print(f"Rounding to two significant figures, the result is {rounded_average_ibi} s.")

calculate_average_ibi()