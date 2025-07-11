import numpy as np

def calculate_average_ibi():
    """
    This function calculates the average inter-burst interval from a list of burst timestamps.
    """
    # Step 1: Identify the timestamps of synchronized bursts from the graph.
    # The timestamps are estimated by finding the center of each dense vertical band of activity.
    burst_timestamps = [2, 13, 22, 34, 45, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (IBIs).
    # IBI is the difference in time between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    
    print("\nCalculating Inter-Burst Intervals (IBIs):")
    # Using a for loop to print the calculation for each interval
    equation_parts = []
    for i in range(len(inter_burst_intervals)):
      interval = burst_timestamps[i+1] - burst_timestamps[i]
      print(f"  Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {interval} s")
      equation_parts.append(str(interval))
    
    # Step 3: Calculate the average IBI.
    average_ibi = np.mean(inter_burst_intervals)

    # Step 4: Display the final calculation and round to two significant figures.
    sum_of_intervals = sum(inter_burst_intervals)
    num_intervals = len(inter_burst_intervals)
    
    print("\nCalculating the average IBI:")
    print(f"  Average IBI = ({' + '.join(equation_parts)}) / {num_intervals}")
    print(f"  Average IBI = {sum_of_intervals} / {num_intervals}")
    print(f"  Average IBI = {average_ibi:.3f} s")

    # Round to two significant figures
    # To round to 'n' significant figures, you can use a format string.
    # For a number like 9.4, '{:.1f}' gives one decimal place, which is two sig figs.
    # For a number like 10.5, '{:.2f}' would be needed. A general solution is more complex.
    # For 9.4, formatting to 1 decimal place is equivalent to 2 significant figures.
    final_answer = float(f"{average_ibi:.1f}")

    print(f"\nThe average inter-burst interval rounded to two significant figures is: {final_answer} s")

calculate_average_ibi()