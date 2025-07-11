import numpy as np

def calculate_average_ibi():
    """
    Identifies burst timestamps from a neuronal firing graph, calculates the
    inter-burst intervals (IBIs), and computes their average.
    """
    # Step 1 & 2: Identify the timestamps of the synchronized bursts from the graph.
    # By visual inspection, the dense vertical bands of activity are centered at these approximate times:
    burst_timestamps = [2, 9, 16, 25, 35, 45]
    print(f"Identified burst timestamps (s): {burst_timestamps}\n")

    # Step 3: Calculate the inter-burst intervals (IBIs).
    # The IBI is the time difference between consecutive bursts.
    inter_burst_intervals = np.diff(burst_timestamps)
    
    print("Calculating inter-burst intervals:")
    for i in range(len(burst_timestamps) - 1):
        print(f"Interval {i+1}: {burst_timestamps[i+1]}s - {burst_timestamps[i]}s = {inter_burst_intervals[i]}s")
    
    print(f"\nCalculated inter-burst intervals (s): {list(inter_burst_intervals)}\n")

    # Step 4: Calculate the average inter-burst interval.
    sum_of_intervals = np.sum(inter_burst_intervals)
    number_of_intervals = len(inter_burst_intervals)
    average_ibi = np.mean(inter_burst_intervals)

    print("Calculating the average inter-burst interval:")
    
    # Create the string for the sum part of the equation
    sum_equation_str = " + ".join(map(str, inter_burst_intervals))
    print(f"Average = ({sum_equation_str}) / {number_of_intervals}")
    print(f"Average = {sum_of_intervals} / {number_of_intervals} = {average_ibi:.2f}s\n")
    
    # Step 5: Round the result to two significant figures.
    # The value 8.6 already has two significant figures.
    final_answer = float('{:.2g}'.format(average_ibi))

    print(f"The average inter-burst interval rounded to two significant figures is: {final_answer}s")

calculate_average_ibi()