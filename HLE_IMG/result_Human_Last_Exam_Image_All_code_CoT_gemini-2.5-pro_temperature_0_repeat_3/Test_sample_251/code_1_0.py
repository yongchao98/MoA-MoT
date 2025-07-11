def calculate_rr_interval():
    """
    Calculates the longest R-R interval in seconds from an ECG strip.

    The function uses the standard ECG paper speed and the measured interval
    in small squares to compute the time duration.
    """
    # ECG paper speed is 25 mm/sec.
    # One small square on the grid is 1 mm.
    # Time per small square = 1 mm / 25 mm/sec = 0.04 seconds.
    time_per_small_square = 0.04  # seconds

    # By visual inspection, the longest R-R interval is found in the fourth strip.
    # It spans 6 large squares and 2 small squares.
    # 1 large square = 5 small squares.
    num_large_squares = 6
    num_small_squares_extra = 2
    total_small_squares = (num_large_squares * 5) + num_small_squares_extra

    # Calculate the longest R-R interval in seconds.
    longest_rr_interval_sec = total_small_squares * time_per_small_square

    print(f"The longest R-R interval spans {total_small_squares} small squares.")
    print(f"Each small square represents {time_per_small_square} seconds.")
    print(f"Calculation: {total_small_squares} squares * {time_per_small_square} sec/square = {longest_rr_interval_sec:.2f} seconds.")
    print(f"The longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")

calculate_rr_interval()