import numpy as np

def calculate_ibi():
    """
    Identifies timestamps of synchronized bursts from the graph,
    calculates the inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1 & 2: Identify and list the estimated timestamps of synchronized bursts in seconds.
    burst_timestamps = [4, 17, 36, 48]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 3: Calculate the inter-burst intervals.
    # An IBI is the difference between consecutive burst times.
    inter_burst_intervals = []
    for i in range(len(burst_timestamps) - 1):
        interval = burst_timestamps[i+1] - burst_timestamps[i]
        inter_burst_intervals.append(interval)

    print(f"Calculated Inter-Burst Intervals (s): {inter_burst_intervals}")

    # Step 4: Calculate the average IBI.
    average_ibi = np.mean(inter_burst_intervals)

    # Display the equation used for the average calculation
    intervals_sum_str = " + ".join(map(str, inter_burst_intervals))
    print(f"Average IBI calculation: ({intervals_sum_str}) / {len(inter_burst_intervals)}")
    print(f"Result before rounding: {average_ibi:.4f} s")
    
    # Step 5: Format the final result to two significant figures.
    # Using the 'g' format specifier handles significant figures correctly.
    average_ibi_formatted = '{:.2g}'.format(average_ibi)

    print(f"\nThe average inter-burst interval is {average_ibi_formatted} s (rounded to two significant figures).")
    
    return average_ibi_formatted

# Run the calculation and print the final answer in the specified format.
final_answer = calculate_ibi()
# The final answer is commented out to avoid execution in the user's environment,
# but it shows the value that should be extracted.
# print(f"<<<{final_answer}>>>")