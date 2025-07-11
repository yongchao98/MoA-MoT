import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts from a raster plot,
    calculates the inter-burst intervals (IBI), and computes the average IBI.
    """
    
    # Step 1: Identify the timestamps of the synchronized bursts from the graph.
    # These are the approximate times (in seconds) where dense vertical bands of activity are observed.
    burst_timestamps = [3, 18, 36, 49]
    print(f"Step 1: The identified timestamps of the synchronized bursts are {burst_timestamps[0]} s, {burst_timestamps[1]} s, {burst_timestamps[2]} s, and {burst_timestamps[3]} s.")
    print("-" * 30)

    # Step 2: Calculate the inter-burst intervals (the time between consecutive bursts).
    inter_burst_intervals = np.diff(burst_timestamps)
    print("Step 2: Calculate the inter-burst intervals (IBIs).")
    for i in range(len(inter_burst_intervals)):
        t1 = burst_timestamps[i]
        t2 = burst_timestamps[i+1]
        interval = inter_burst_intervals[i]
        print(f"  Interval {i+1}: {t2} s - {t1} s = {interval} s")
    print("-" * 30)

    # Step 3: Calculate the average inter-burst interval.
    average_ibi = np.mean(inter_burst_intervals)
    print("Step 3: Calculate the average IBI.")
    # Display the full calculation equation
    interval_sum_str = " + ".join(map(str, inter_burst_intervals))
    print(f"  Average IBI = ({interval_sum_str}) / {len(inter_burst_intervals)}")
    print(f"  Average IBI = {sum(inter_burst_intervals)} / {len(inter_burst_intervals)} = {average_ibi:.4f} s")
    print("-" * 30)
    
    # Step 4: Round the average IBI to two significant figures.
    # For a number like 15.33..., the first two significant digits are 1 and 5.
    final_answer = float('{:.2g}'.format(average_ibi))
    print("Step 4: Round the result to two significant figures.")
    print(f"The average inter-burst interval is {final_answer} s.")
    
    return final_answer

if __name__ == "__main__":
    calculate_average_ibi()
