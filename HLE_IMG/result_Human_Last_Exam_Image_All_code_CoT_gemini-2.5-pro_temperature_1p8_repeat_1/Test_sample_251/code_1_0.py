import math

def calculate_longest_rr_interval():
    """
    Calculates the longest R-R interval from the provided ECG strip information.

    The steps are:
    1. Identify the paper speed.
    2. Measure the longest R-R interval in terms of grid squares from the image.
    3. Convert the grid measurement to time in seconds.
    """

    # 1. ECG Parameters
    paper_speed_mm_per_sec = 25  # mm/sec, as stated on the ECG
    mm_per_small_square = 1
    small_squares_per_large_square = 5

    # 2. Measurement from ECG
    # By visually inspecting the long pause in the fourth strip, we measure the distance.
    # It spans approximately 8 large squares and 3 extra small squares.
    num_large_squares = 8
    num_extra_small_squares = 3

    # 3. Calculation
    # Calculate the time represented by one small square
    time_per_small_square = mm_per_small_square / paper_speed_mm_per_sec

    # Calculate the total number of small squares for the interval
    total_small_squares = (num_large_squares * small_squares_per_large_square) + num_extra_small_squares

    # Calculate the final R-R interval duration in seconds
    longest_rr_interval_sec = total_small_squares * time_per_small_square

    # Print the explanation and the result
    print(f"The paper speed is {paper_speed_mm_per_sec} mm/sec.")
    print(f"One small 1mm square represents 1/{paper_speed_mm_per_sec} = {time_per_small_square:.2f} seconds.")
    print(f"The longest pause on the ECG measures approximately {num_large_squares} large squares and {num_extra_small_squares} small squares.")
    print(f"The total distance in small squares is ({num_large_squares} * {small_squares_per_large_square}) + {num_extra_small_squares} = {total_small_squares} squares.")
    print(f"The calculation for the longest R-R interval is: {total_small_squares} squares * {time_per_small_square:.2f} sec/square.")
    print(f"The longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")

calculate_longest_rr_interval()