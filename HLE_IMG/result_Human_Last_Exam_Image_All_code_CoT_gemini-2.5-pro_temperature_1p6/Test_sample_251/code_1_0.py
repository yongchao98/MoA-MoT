def calculate_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function uses standard ECG paper measurements and visual analysis
    of the provided image to determine the duration of the longest R-R interval.
    """

    # ECG paper speed in mm/sec, as stated on the strip
    paper_speed_mm_per_sec = 25

    # Each small square on the grid is 1 mm.
    # Time for one small square (sec) = 1 mm / (paper_speed_mm_per_sec)
    time_per_small_square_sec = 1 / paper_speed_mm_per_sec

    # From visual inspection of the ECG, the longest R-R interval spans
    # approximately 36 small squares. This is observed in the fourth strip,
    # where a pause spans 7 large squares (7 * 5 = 35 small squares) plus
    # one additional small square.
    longest_interval_in_small_squares = 36

    # Calculate the total duration of the interval in seconds.
    longest_rr_interval_sec = longest_interval_in_small_squares * time_per_small_square_sec

    print("ECG Analysis for Longest R-R Interval")
    print("---------------------------------------")
    print(f"Paper Speed: {paper_speed_mm_per_sec} mm/sec")
    print(f"Time per small square: 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec:.2f} seconds")
    print(f"Longest measured interval in small squares: {longest_interval_in_small_squares}")
    print("\nCalculation:")
    print(f"{longest_interval_in_small_squares} small squares * {time_per_small_square_sec:.2f} sec/square = {longest_rr_interval_sec:.2f} seconds")
    print("\nFinal Answer:")
    print(f"The longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")

calculate_rr_interval()