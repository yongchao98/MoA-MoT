def calculate_rr_interval():
    """
    Calculates the longest R-R interval in seconds from an ECG strip.
    """
    # Step 1: Define the constants from the ECG.
    paper_speed_mm_per_sec = 25
    # One small square on the grid is 1 mm.
    seconds_per_small_square = 1 / paper_speed_mm_per_sec

    # Step 2: Measure the longest R-R interval in small squares from the image.
    # The longest interval is observed in the fourth strip.
    # It spans 22 large squares and 1 additional small square.
    # Each large square is 5 small squares.
    large_squares = 22
    small_squares_per_large_square = 5
    additional_small_squares = 1
    total_small_squares = (large_squares * small_squares_per_large_square) + additional_small_squares

    # Step 3: Calculate the interval in seconds.
    longest_rr_interval_seconds = total_small_squares * seconds_per_small_square

    # Step 4: Print the explanation and the result.
    print("To find the longest R-R interval:")
    print(f"1. The paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"2. Each small 1 mm square on the horizontal axis is 1/{paper_speed_mm_per_sec} = {seconds_per_small_square:.2f} seconds.")
    print("3. The longest interval on the strip is measured to be 22 large squares and 1 small square.")
    print(f"4. The total number of small squares is ({large_squares} * {small_squares_per_large_square}) + {additional_small_squares} = {total_small_squares} small squares.")
    print(f"5. The duration in seconds is {total_small_squares} squares * {seconds_per_small_square:.2f} sec/square.")
    print(f"\nFinal Equation: {total_small_squares} * {seconds_per_small_square:.2f} = {longest_rr_interval_seconds:.2f}")
    print(f"The longest R-R interval is {longest_rr_interval_seconds:.2f} seconds.")

calculate_rr_interval()