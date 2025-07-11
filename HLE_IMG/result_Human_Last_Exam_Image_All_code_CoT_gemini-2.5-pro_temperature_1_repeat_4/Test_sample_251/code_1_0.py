def calculate_rr_interval():
    """
    This function calculates the longest R-R interval from an ECG strip.

    The standard ECG paper speed is 25 mm/sec.
    The grid consists of large squares, each 5 mm wide.
    Therefore, the time for one large square is 5 mm / 25 mm/sec = 0.2 seconds.

    By inspecting the ECG strip, the longest interval between two R-waves
    is found to be 8 large squares.
    """
    
    # ECG parameters
    time_per_large_square = 0.2  # seconds
    
    # Measurement from the ECG strip
    num_large_squares = 8
    
    # Calculation
    longest_rr_interval = num_large_squares * time_per_large_square
    
    # Output the result
    print(f"The longest R-R interval is determined by measuring the number of large squares between the two R-waves of the longest pause.")
    print(f"Number of large squares = {num_large_squares}")
    print(f"Time per large square = {time_per_large_square} seconds")
    print(f"Calculation: {num_large_squares} large squares * {time_per_large_square} seconds/large square = {longest_rr_interval:.1f} seconds")

calculate_rr_interval()