def calculate_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function is based on the standard parameters of an ECG recording.
    The paper speed is 25 mm/sec.
    The time for one small square (1 mm) is calculated from the paper speed.
    The longest R-R interval is measured by counting the number of small squares
    between the two R-waves that are furthest apart.
    """
    
    # ECG paper speed in mm/sec
    paper_speed_mm_per_sec = 25
    
    # Duration of one small square (1 mm) in seconds
    small_square_duration_sec = 1 / paper_speed_mm_per_sec
    
    # By visual inspection, the longest R-R interval is measured.
    # It spans approximately 4 large squares and 1 small square.
    # Number of small squares in a large square = 5
    # Total small squares = (4 * 5) + 1 = 21
    num_small_squares = 21
    
    # Calculate the longest R-R interval in seconds
    rr_interval_seconds = num_small_squares * small_square_duration_sec
    
    # Print the explanation and the calculation
    print("To find the longest R-R interval in seconds:")
    print(f"1. The paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"2. The duration of one small square (1 mm) is 1 / {paper_speed_mm_per_sec} = {small_square_duration_sec} seconds.")
    print(f"3. The longest R-R interval on the strip measures {num_small_squares} small squares.")
    print(f"4. The calculation is: {num_small_squares} small squares * {small_square_duration_sec} seconds/small square = {rr_interval_seconds:.2f} seconds.")
    print(f"\nThe final equation is: {num_small_squares} * {small_square_duration_sec} = {rr_interval_seconds:.2f}")

calculate_rr_interval()