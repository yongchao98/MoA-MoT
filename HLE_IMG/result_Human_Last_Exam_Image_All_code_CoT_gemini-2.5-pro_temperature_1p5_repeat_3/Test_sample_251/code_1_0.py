import math

def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval from an ECG strip.

    This function determines the longest R-R interval based on standard ECG parameters
    and visual analysis of the provided image.
    """

    # Step 1: Define ECG parameters from the image.
    paper_speed_mm_per_sec = 25  # The paper speed is 25 mm/sec.
    large_square_width_mm = 5    # The width of one large square is 5 mm.

    # Step 2: Calculate the time duration of one large square.
    time_per_large_square = large_square_width_mm / paper_speed_mm_per_sec

    # Step 3: Visually identify and measure the longest R-R interval from the image.
    # By inspecting the ECG strip, the longest interval is found on the fourth line.
    # Counting the squares between the peaks of the two relevant R-waves reveals
    # that the interval covers 7 large squares.
    num_large_squares = 7

    # Step 4: Calculate the total duration of the longest R-R interval.
    longest_rr_interval_sec = num_large_squares * time_per_large_square

    # Step 5: Print the detailed explanation and the final calculation.
    print("To determine the longest R-R interval:")
    print(f"1. The ECG paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"2. Each large square is {large_square_width_mm} mm, which corresponds to {time_per_large_square} seconds ({large_square_width_mm} mm / {paper_speed_mm_per_sec} mm/sec).")
    print(f"3. The longest interval found on the strip spans {num_large_squares} large squares.")
    print("\nFinal Calculation:")
    print(f"{num_large_squares} large squares * {time_per_large_square} sec/large square = {longest_rr_interval_sec} seconds")
    print(f"\nThe longest R-R interval is {longest_rr_interval_sec} seconds.")

# Execute the function
calculate_longest_rr_interval()