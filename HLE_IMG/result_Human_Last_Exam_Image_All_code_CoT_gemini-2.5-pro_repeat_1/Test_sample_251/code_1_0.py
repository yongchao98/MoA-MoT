def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval in seconds from an ECG strip.
    """
    # 1. ECG parameters
    paper_speed_mm_per_sec = 25  # Standard paper speed in mm/sec

    # 2. Time conversion
    # Each small square on the ECG grid is 1 mm.
    time_per_small_square_sec = 1 / paper_speed_mm_per_sec

    # 3. Measurement from the ECG strip
    # By visual inspection, the longest R-R interval is found on the fourth line.
    # It spans 9 large squares and 4 small squares.
    # Each large square is 5 small squares.
    num_large_squares = 9
    num_additional_small_squares = 4
    total_small_squares = (num_large_squares * 5) + num_additional_small_squares

    # 4. Final calculation
    longest_rr_interval_sec = total_small_squares * time_per_small_square_sec

    # Print the explanation and result
    print("Step 1: The paper speed is 25 mm/sec.")
    print(f"Step 2: The time for one small square (1 mm) is 1 / {paper_speed_mm_per_sec} = {time_per_small_square_sec} seconds.")
    print(f"Step 3: The longest R-R interval on the strip spans {total_small_squares} small squares.")
    print(f"Step 4: The calculation is: {total_small_squares} squares * {time_per_small_square_sec} sec/square.")
    print(f"The longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")

calculate_longest_rr_interval()