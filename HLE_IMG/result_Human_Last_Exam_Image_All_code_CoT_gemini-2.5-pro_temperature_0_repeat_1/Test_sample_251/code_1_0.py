def calculate_rr_interval():
    """
    Calculates the longest R-R interval in seconds from an ECG strip.
    """
    # ECG paper speed in mm/sec
    paper_speed_mm_per_sec = 25

    # Time for one small 1mm square in seconds
    time_per_small_square_sec = 1 / paper_speed_mm_per_sec

    # Visually measured number of small squares for the longest R-R interval
    num_small_squares = 41

    # Calculate the longest R-R interval in seconds
    longest_rr_interval_sec = num_small_squares * time_per_small_square_sec

    # Print the calculation steps and the final result
    print(f"The ECG paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"The time for one small square is 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec} seconds.")
    print(f"The longest R-R interval spans {num_small_squares} small squares.")
    print("\nCalculation:")
    print(f"{num_small_squares} small squares * {time_per_small_square_sec:.2f} seconds/square = {longest_rr_interval_sec:.2f} seconds")

calculate_rr_interval()