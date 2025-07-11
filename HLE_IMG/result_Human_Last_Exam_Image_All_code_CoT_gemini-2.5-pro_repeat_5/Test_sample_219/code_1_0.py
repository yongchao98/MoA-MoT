import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from the provided graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1 & 2: Identify and list the timestamps of the synchronized bursts (in seconds)
    burst_timestamps = [3, 12, 18, 27, 37, 48]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 3: Calculate the inter-burst intervals
    # np.diff calculates the difference between consecutive elements in an array
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")

    # Step 4: Calculate the average inter-burst interval
    average_ibi = np.mean(inter_burst_intervals)
    
    # Step 5 & 6: Format the output and print the full calculation
    # Create the equation string
    intervals_str = " + ".join(map(str, inter_burst_intervals))
    num_intervals = len(inter_burst_intervals)
    
    # Format the average to two significant figures
    # The format specifier '.1f' works for this specific result (9.0)
    # A more general approach would be needed for other numbers, but this is sufficient here.
    average_ibi_formatted = f"{average_ibi:.1f}"

    print("\nCalculation of the average inter-burst interval:")
    print(f"({intervals_str}) / {num_intervals} = {average_ibi_formatted} s")
    
    # Print the final numerical answer separately for clarity
    print(f"\nThe average inter-burst interval is {average_ibi_formatted} s.")

calculate_average_ibi()