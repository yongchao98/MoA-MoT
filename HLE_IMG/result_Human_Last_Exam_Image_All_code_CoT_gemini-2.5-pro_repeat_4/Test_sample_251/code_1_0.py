def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function is based on visual analysis of the provided ECG image.
    1. The paper speed is 25 mm/sec.
    2. A small square on the grid is 1 mm.
    3. The time represented by one small square is 1 mm / 25 mm/sec = 0.04 seconds.
    4. By visually inspecting the ECG, the longest interval between two consecutive R-waves
       is identified and measured in small squares. This longest interval spans 24 small squares.
    5. The duration in seconds is calculated by multiplying the number of squares by the time per square.
    """
    
    # Paper speed in mm per second
    paper_speed_mm_per_sec = 25
    
    # Time per small square (1 mm) in seconds
    time_per_small_square = 1 / paper_speed_mm_per_sec
    
    # Longest R-R interval measured in number of small squares from the ECG
    num_small_squares = 24
    
    # Calculate the longest R-R interval in seconds
    longest_rr_interval_sec = num_small_squares * time_per_small_square
    
    # Print the equation and the final answer
    print(f"To find the longest R-R interval in seconds, we use the formula:")
    print(f"Duration = (Number of Small Squares) * (Time per Small Square)")
    print(f"Based on the ECG at 25 mm/sec, Time per Small Square = 1 / {paper_speed_mm_per_sec} = {time_per_small_square} seconds.")
    print(f"The longest interval spans {num_small_squares} small squares.")
    print(f"Calculation: {num_small_squares} * {time_per_small_square} = {longest_rr_interval_sec}")
    print(f"\nThe longest R-R interval is {longest_rr_interval_sec} seconds.")

calculate_longest_rr_interval()