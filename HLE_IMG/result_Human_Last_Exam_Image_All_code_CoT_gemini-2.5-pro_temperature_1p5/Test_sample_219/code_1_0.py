import numpy as np

def calculate_average_ibi():
    """
    Identifies timestamps of synchronized bursts, calculates the inter-burst intervals (IBIs),
    and computes the average IBI.
    """
    # Step 1: Identify timestamps of synchronized bursts from the graph
    burst_timestamps = [5, 17, 28, 37, 49]
    print(f"Identified burst timestamps (s): {burst_timestamps}")

    # Step 2: Calculate the inter-burst intervals (IBIs)
    inter_burst_intervals = np.diff(burst_timestamps)
    print(f"Calculated inter-burst intervals (s): {list(inter_burst_intervals)}")

    # Step 3: Calculate the average IBI
    average_ibi = np.mean(inter_burst_intervals)

    # Step 4: Print the calculation and the final result
    # Constructing the equation string
    intervals_str = ' + '.join(map(str, inter_burst_intervals))
    equation = f"({intervals_str}) / {len(inter_burst_intervals)}"
    
    print(f"\nCalculation for average inter-burst interval:")
    # Using the print function to display each number in the equation
    print(f"Average = ({inter_burst_intervals[0]} + {inter_burst_intervals[1]} + {inter_burst_intervals[2]} + {inter_burst_intervals[3]}) / 4 = {average_ibi:.2f} s")
    
    # Rounding to two significant figures for the final answer.
    # For a number like 11.00, two significant figures is 11.
    final_answer = float(f"{average_ibi:.1g}") if average_ibi >= 10 else float(f"{average_ibi:.2g}")

    print(f"\nThe average inter-burst interval is {final_answer} s (rounded to two significant figures).")
    
    return final_answer

if __name__ == '__main__':
    result = calculate_average_ibi()
    # The final answer format is handled by the platform,
    # but we can print it here for clarity.
    # print(f"<<<{result}>>>")

calculate_average_ibi()