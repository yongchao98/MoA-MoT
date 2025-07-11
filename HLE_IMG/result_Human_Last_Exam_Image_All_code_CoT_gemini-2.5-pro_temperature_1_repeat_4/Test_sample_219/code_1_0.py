import numpy as np

def calculate_average_ibi():
    """
    Identifies burst timestamps from a neuronal activity graph, calculates the
    inter-burst intervals (IBIs), and computes the average IBI.
    """
    # Step 1 & 2: Identify and estimate the timestamps of the synchronized bursts.
    # These are the dense vertical bands of firing activity across many electrodes.
    burst_timestamps = [2.5, 9.0, 17.5, 25.5, 36.5, 47.5]
    print(f"Identified burst timestamps (s): {burst_timestamps}\n")

    # Step 3: Calculate the inter-burst intervals.
    inter_burst_intervals = np.diff(burst_timestamps)
    print("Calculating the inter-burst intervals (time between consecutive bursts):")
    for i in range(len(inter_burst_intervals)):
        print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {inter_burst_intervals[i]:.1f} s")
    print(f"\nCalculated inter-burst intervals (s): {list(inter_burst_intervals)}\n")

    # Step 4: Calculate the average inter-burst interval.
    average_ibi = np.mean(inter_burst_intervals)
    
    # Step 5: Display the final calculation and round the result.
    intervals_sum_str = " + ".join(map(str, [f"{i:.1f}" for i in inter_burst_intervals]))
    num_intervals = len(inter_burst_intervals)
    sum_intervals = sum(inter_burst_intervals)

    print("Calculating the average inter-burst interval:")
    print(f"Average = ({intervals_sum_str}) / {num_intervals}")
    print(f"Average = {sum_intervals:.1f} / {num_intervals}")
    
    # The final answer is rounded to two significant figures.
    final_answer = f"{average_ibi:.1f}"
    print(f"Average = {final_answer} s")

    return final_answer

if __name__ == '__main__':
    # Run the calculation and store the final answer.
    result = calculate_average_ibi()
    # The final answer format is specified by the problem.
    # The instruction says "directly return the answer with the format <<<answer content>>>"
    # But since the prompt also says "Instead, use 'print' function for the output when relevant."
    # and the final format is just a number, I will print the result as requested,
    # and then add the special format at the end.
    # Let's assume the special format is the final expected output.
    # Final answer is 9.0 which has two significant figures.
