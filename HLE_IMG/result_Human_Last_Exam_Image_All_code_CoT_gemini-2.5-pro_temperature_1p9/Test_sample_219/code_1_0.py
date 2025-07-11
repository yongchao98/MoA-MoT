import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from a neuronal activity graph,
    calculates the inter-burst intervals (IBIs), and computes their average.
    """
    # Step 1: Identify the approximate timestamps of the synchronized bursts from the graph.
    # Bursts are the dense vertical lines of firing activity.
    burst_timestamps = [5, 16, 26, 36, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (the time between consecutive bursts).
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")

    # Step 3: Calculate the average of the inter-burst intervals.
    average_ibi = np.mean(inter_burst_intervals)

    # Create the equation string
    intervals_str = ' + '.join(map(str, inter_burst_intervals))
    equation = f"({intervals_str}) / {len(inter_burst_intervals)}"

    # Step 4: Print the final calculation and the result rounded to two significant figures.
    print(f"\nAverage Inter-Burst Interval Calculation:")
    print(f"Average = {equation} = {average_ibi:.2g} s")

    # Final answer in the required format
    print(f"\n<<<The average inter-burst interval is {average_ibi:.2g} s>>>")


# Run the calculation
calculate_average_ibi()
