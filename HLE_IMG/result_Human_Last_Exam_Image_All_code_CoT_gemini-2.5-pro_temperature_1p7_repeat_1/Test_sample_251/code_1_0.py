import math

def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function determines the longest R-R interval in seconds based on visual
    analysis of the provided ECG image.

    Steps:
    1. The paper speed is given as 25 mm/sec. Each small square on the grid is 1 mm.
       Therefore, the time for one small square is 1 sec / 25 mm = 0.04 seconds.
    2. By visual inspection of the ECG strip, the longest pause is located on the
       fourth tracing.
    3. We count the number of small squares between the R-wave peaks for this interval.
       The interval spans approximately 34 small squares.
    4. The total time for the longest R-R interval is calculated by multiplying
       the number of squares by the time per square.
    """
    
    paper_speed_mm_per_sec = 25
    
    # Calculate the time represented by one small square (1 mm)
    time_per_small_square = 1 / paper_speed_mm_per_sec
    
    # Count the number of small squares for the longest R-R interval from the image
    longest_interval_in_squares = 34
    
    # Calculate the longest R-R interval in seconds
    longest_rr_interval_seconds = longest_interval_in_squares * time_per_small_square
    
    # Print the explanation and the calculation
    print("Step 1: Determine the time duration of one small square on the ECG grid.")
    print(f"Given the paper speed of {paper_speed_mm_per_sec} mm/sec, each 1 mm small square represents 1 / {paper_speed_mm_per_sec} = {time_per_small_square:.2f} seconds.")
    print("\nStep 2: Measure the longest R-R interval on the ECG strip in small squares.")
    print(f"By visual inspection, the longest interval between two R-waves spans {longest_interval_in_squares} small squares.")
    print("\nStep 3: Calculate the duration of the longest R-R interval in seconds.")
    print("The final calculation is:")
    print(f"{longest_interval_in_squares} squares * {time_per_small_square:.2f} seconds/square = {longest_rr_interval_seconds:.2f} seconds")

calculate_longest_rr_interval()