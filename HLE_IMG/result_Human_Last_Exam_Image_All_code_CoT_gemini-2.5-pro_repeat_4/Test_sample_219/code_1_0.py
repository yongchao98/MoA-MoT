import numpy as np

def calculate_average_ibi():
    """
    This function identifies the timestamps of synchronized bursts from the graph,
    calculates the inter-burst intervals, and computes their average.
    """
    # Step 1: Identify the approximate center timestamp for each synchronized burst.
    # These are estimated visually from the dense vertical columns of activity in the graph.
    burst_timestamps = [3, 17, 26, 36, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}\n")

    # Step 2: Calculate the intervals between consecutive bursts.
    intervals = []
    print("Calculating inter-burst intervals:")
    for i in range(len(burst_timestamps) - 1):
        interval = burst_timestamps[i+1] - burst_timestamps[i]
        intervals.append(interval)
        # The prompt requires showing the numbers in the final equation,
        # but showing the individual subtractions is also helpful for clarity.
        print(f"  Interval {i+1}: {burst_timestamps[i+1]} s - {burst_timestamps[i]} s = {interval} s")

    # Step 3: Calculate the average of the intervals.
    average_interval = np.mean(intervals)
    
    # Create the equation string as requested
    interval_strings = [str(i) for i in intervals]
    equation = f"({ ' + '.join(interval_strings) }) / {len(intervals)}"

    print(f"\nAverage interval calculation: {equation} = {average_interval:.3f} s")

    # Step 4: Round the final answer to two significant figures.
    # For a number like 11.5, two significant figures is 12.
    average_rounded = float(f"{average_interval:.2g}")
    print(f"\nThe average inter-burst interval, rounded to two significant figures, is: {average_rounded} s")

if __name__ == "__main__":
    calculate_average_ibi()