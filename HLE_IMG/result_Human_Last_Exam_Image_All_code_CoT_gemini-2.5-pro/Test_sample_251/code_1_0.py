def calculate_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function uses the standard ECG paper speed and grid measurements
    to find the duration of the longest interval identified on the strip.
    """

    # ECG parameters
    paper_speed_mm_per_sec = 25
    mm_per_large_square = 5

    # Time calculation
    time_per_mm = 1 / paper_speed_mm_per_sec  # seconds
    time_per_large_square = time_per_mm * mm_per_large_square # seconds

    # Measurement from the ECG strip
    # The longest interval is visually identified and measured to be 7 large squares.
    longest_interval_in_large_squares = 7

    # Calculate the final duration in seconds
    longest_rr_interval_sec = longest_interval_in_large_squares * time_per_large_square

    print(f"ECG Paper Speed: {paper_speed_mm_per_sec} mm/sec")
    print(f"Time per large square: {time_per_large_square:.2f} seconds")
    print(f"Longest measured interval in large squares: {longest_interval_in_large_squares}")
    print("\nCalculation:")
    print(f"{longest_interval_in_large_squares} large squares * {time_per_large_square:.2f} sec/large square = {longest_rr_interval_sec:.2f} seconds")

    print(f"\nThe longest R-R interval is: {longest_rr_interval_sec:.2f} seconds")

calculate_rr_interval()