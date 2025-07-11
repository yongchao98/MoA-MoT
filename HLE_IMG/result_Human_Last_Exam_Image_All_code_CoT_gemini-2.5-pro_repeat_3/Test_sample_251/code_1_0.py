# Plan:
# 1. Define the time represented by one small square on the ECG paper.
#    The paper speed is 25 mm/sec, and a small square is 1 mm wide.
#    Time per small square = 1 mm / 25 mm/sec = 0.04 seconds.
# 2. Measure the longest R-R interval in terms of small squares from the image.
#    By visual inspection, the longest pause spans 7 large squares and 1 small square.
#    Since 1 large square = 5 small squares, the total is (7 * 5) + 1 = 36 small squares.
# 3. Calculate the interval in seconds by multiplying the number of squares by the time per square.

# 1. Define constants
time_per_small_square_sec = 0.04  # seconds

# 2. Measurement from the ECG
num_large_squares = 7
num_extra_small_squares = 1
small_squares_per_large_square = 5
total_small_squares = (num_large_squares * small_squares_per_large_square) + num_extra_small_squares

# 3. Calculation
longest_rr_interval_sec = total_small_squares * time_per_small_square_sec

# Output the calculation and the result
print(f"The longest R-R interval spans {total_small_squares} small squares.")
print(f"Each small square represents {time_per_small_square_sec} seconds.")
print(f"Calculation: {total_small_squares} * {time_per_small_square_sec} = {longest_rr_interval_sec:.2f} seconds.")
print(f"\nThe longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")