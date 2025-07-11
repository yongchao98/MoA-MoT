def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function uses standard ECG paper measurements to calculate the duration
    of the longest R-R interval identified visually on the provided strip.
    """
    # Standard ECG paper speed in mm/sec
    paper_speed_mm_per_sec = 25

    # Dimensions of ECG grid squares
    # 1 large square = 5 small squares (5 mm)
    # 1 small square = 1 mm

    # Time calculation
    # Time for 1 large square (5 mm) = 5 mm / 25 mm/sec
    time_per_large_square_sec = 5 / paper_speed_mm_per_sec

    # From visual inspection of the ECG, the longest R-R interval spans 7 large squares.
    num_large_squares = 7

    # Calculate the longest R-R interval in seconds
    longest_rr_interval_sec = num_large_squares * time_per_large_square_sec

    # Print the explanation and the final result
    print("To find the longest R-R interval:")
    print(f"1. The paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"2. The time for one large square (5 mm) is 5 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_large_square_sec} seconds.")
    print(f"3. By inspecting the ECG, the longest interval between two R-waves spans {num_large_squares} large squares.")
    print(f"4. The duration of the longest R-R interval is calculated as: {num_large_squares} large squares * {time_per_large_square_sec} seconds/large square = {longest_rr_interval_sec} seconds.")

calculate_longest_rr_interval()