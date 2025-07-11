def calculate_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    The calculation is based on the standard ECG paper speed and the
    measurement of the interval in small squares.
    """
    # Standard paper speed in mm/sec
    paper_speed_mm_per_sec = 25

    # Time duration of one small square (1 mm) in seconds
    time_per_small_square_sec = 1 / paper_speed_mm_per_sec

    # The longest R-R interval was measured to be approximately 43 small squares
    num_small_squares = 43

    # Calculate the longest R-R interval in seconds
    longest_rr_interval_sec = num_small_squares * time_per_small_square_sec

    # Print the explanation and the result
    print("To find the longest R-R interval:")
    print(f"1. The paper speed is {paper_speed_mm_per_sec} mm/sec, so each 1mm small square is 1/{paper_speed_mm_per_sec} = {time_per_small_square_sec:.2f} seconds.")
    print(f"2. The longest interval on the strip spans approximately {num_small_squares} small squares.")
    print(f"3. The calculation is: {num_small_squares} squares * {time_per_small_square_sec:.2f} sec/square = {longest_rr_interval_sec:.2f} seconds.")
    print(f"\nThe longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")

calculate_rr_interval()