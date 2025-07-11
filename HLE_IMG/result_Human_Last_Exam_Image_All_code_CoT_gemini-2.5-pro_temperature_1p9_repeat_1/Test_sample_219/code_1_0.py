import numpy as np

def calculate_average_ibi():
    """
    This function calculates the average inter-burst interval from a list of burst timestamps.
    """
    # Step 1: Identify the timestamps of synchronized bursts from the graph.
    burst_timestamps = [4, 16, 27, 37, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (IBIs).
    inter_burst_intervals = []
    print("\nCalculating Inter-Burst Intervals (IBIs):")
    for i in range(len(burst_timestamps) - 1):
        interval = burst_timestamps[i+1] - burst_timestamps[i]
        inter_burst_intervals.append(interval)
        print(f"Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {interval} s")

    print(f"\nCalculated IBIs (s): {inter_burst_intervals}")

    # Step 3: Calculate the average IBI.
    average_ibi = np.mean(inter_burst_intervals)
    
    # Create the equation string for the average calculation
    interval_sum_str = " + ".join(map(str, inter_burst_intervals))
    average_equation = f"Average IBI = ({interval_sum_str}) / {len(inter_burst_intervals)} = {average_ibi} s"
    print("\nCalculating the average IBI:")
    print(average_equation)
    
    # Step 4: Round the result to two significant figures.
    # To round to two significant figures, we can use a format specifier.
    # 11.25 has four significant figures. Rounding to two gives 11.
    final_answer = float(f"{average_ibi:.2g}")
    print(f"\nAverage IBI rounded to two significant figures: {final_answer} s")

    # This part is just to fulfill the final answer format requirement.
    # The actual answer is printed above.
    return final_answer

if __name__ == "__main__":
    result = calculate_average_ibi()
    # The final answer format will be printed outside the standard output logic
    # print(f"<<<{result}>>>")