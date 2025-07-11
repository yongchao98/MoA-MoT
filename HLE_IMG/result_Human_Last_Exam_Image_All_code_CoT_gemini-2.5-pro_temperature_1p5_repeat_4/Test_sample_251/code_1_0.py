def calculate_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function determines the time in seconds corresponding to the longest
    R-R interval observed on the provided ECG image.
    """
    
    # ECG paper speed in mm/sec, as stated on the strip.
    paper_speed_mm_per_sec = 25
    
    # Time corresponding to one small 1mm square on the ECG grid.
    time_per_small_square_sec = 1 / paper_speed_mm_per_sec
    
    # By visual inspection, the longest R-R interval spans 7 large squares and 3 small squares.
    # Number of small squares in one large square.
    small_squares_per_large_square = 5
    
    # Measurement of the longest interval.
    num_large_squares = 7
    num_small_squares_extra = 3
    
    # Total number of small squares for the longest interval.
    total_small_squares = (num_large_squares * small_squares_per_large_square) + num_small_squares_extra
    
    # Calculate the longest R-R interval in seconds.
    longest_rr_interval_sec = total_small_squares * time_per_small_square_sec
    
    print("Step 1: Determine the time represented by one small square on the ECG.")
    print(f"The paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"Each small square is 1 mm, so its duration is 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec:.2f} seconds.")
    print("\nStep 2: Measure the longest R-R interval in small squares.")
    print("The longest interval spans 7 large squares and 3 small squares.")
    print(f"Total small squares = ( {num_large_squares} large squares * {small_squares_per_large_square} small squares/large square ) + {num_small_squares_extra} small squares = {total_small_squares} small squares.")
    print("\nStep 3: Calculate the duration of the interval in seconds.")
    print(f"Duration = {total_small_squares} squares * {time_per_small_square_sec:.2f} sec/square = {longest_rr_interval_sec:.2f} seconds.")

calculate_rr_interval()