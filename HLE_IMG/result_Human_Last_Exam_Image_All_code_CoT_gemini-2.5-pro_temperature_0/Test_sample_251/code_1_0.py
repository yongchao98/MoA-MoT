def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval in seconds from an ECG strip.

    The function uses the standard ECG paper speed and the measured interval
    in small squares to compute the duration.
    """
    # Step 1: Define the given parameters from the ECG.
    paper_speed_mm_per_sec = 25  # mm/sec
    small_square_mm = 1          # mm

    # Step 2: Calculate the time duration of one small square.
    time_per_small_square = small_square_mm / paper_speed_mm_per_sec

    # Step 3: Measure the longest R-R interval in small squares from the image.
    # By visual inspection, the longest pause spans 7 large squares and 2 small squares.
    # One large square = 5 small squares.
    num_large_squares = 7
    num_extra_small_squares = 2
    total_small_squares = (num_large_squares * 5) + num_extra_small_squares

    # Step 4: Calculate the longest R-R interval in seconds.
    longest_rr_interval_seconds = total_small_squares * time_per_small_square

    # Print the explanation and the final calculation.
    print("To find the longest R-R interval:")
    print(f"1. The paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"2. The time for one small square (1 mm) is 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square:.2f} seconds.")
    print("3. The longest interval on the strip is measured to be 7 large squares and 2 small squares.")
    print(f"4. In terms of small squares, this is ({num_large_squares} * 5) + {num_extra_small_squares} = {total_small_squares} small squares.")
    print(f"5. The duration in seconds is {total_small_squares} squares * {time_per_small_square:.2f} sec/square.")
    print(f"Final Equation: {total_small_squares} * {time_per_small_square:.2f} = {longest_rr_interval_seconds:.2f}")
    print(f"The longest R-R interval is {longest_rr_interval_seconds:.2f} seconds.")

calculate_longest_rr_interval()