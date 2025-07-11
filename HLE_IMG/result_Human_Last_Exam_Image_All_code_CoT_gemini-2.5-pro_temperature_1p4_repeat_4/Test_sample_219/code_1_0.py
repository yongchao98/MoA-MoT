import numpy as np

def calculate_average_ibi():
    """
    Identifies burst timestamps from the graph, calculates the inter-burst intervals (IBIs),
    and computes the average IBI, rounding the result to two significant figures.
    """
    # Step 1: Identify the timestamps of the synchronized bursts from the graph.
    burst_timestamps = [3, 15, 26, 37, 48]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (the time between consecutive bursts).
    intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {intervals.tolist()}")

    # Step 3: Calculate the average inter-burst interval.
    average_interval = np.mean(intervals)
    
    # Create the equation string for the final output
    # e.g., "Average Interval = (12 + 11 + 11 + 11) / 4 = 11.25"
    interval_sum_str = " + ".join(map(str, intervals))
    equation = f"Average Interval = ({interval_sum_str}) / {len(intervals)} = {average_interval}"
    print(equation)

    # Step 4: Round the final answer to two significant figures.
    # To format to 2 significant figures, we can use a format specifier.
    # For a number like 11.25, two significant figures is 11.
    final_answer = float(f"{average_interval:.2g}")
    print(f"\nThe average inter-burst interval rounded to two significant figures is: {final_answer} s")

    return final_answer

if __name__ == "__main__":
    result = calculate_average_ibi()
    # The final answer is wrapped in <<<>>> as requested.
    # Note: the print statements within the function provide the step-by-step thinking process.
    # The final print here is for the final answer format.
    print(f"\n<<<{result}>>>")
