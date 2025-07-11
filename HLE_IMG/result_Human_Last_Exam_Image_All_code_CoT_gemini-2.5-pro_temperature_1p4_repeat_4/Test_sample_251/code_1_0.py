# Step 1: Define the constants based on the ECG strip information.
paper_speed_mm_per_sec = 25  # The paper speed is 25 mm/sec.
time_per_small_square_sec = 1 / paper_speed_mm_per_sec  # Each small square is 1mm.

# Step 2: Measure the longest R-R interval in terms of grid squares from the image.
# The longest interval is in the fourth strip, spanning about 8 large squares and 2 small squares.
num_large_squares = 8
num_small_squares_in_large_square = 5
num_extra_small_squares = 2

# Step 3: Calculate the total number of small squares for the interval.
total_small_squares = (num_large_squares * num_small_squares_in_large_square) + num_extra_small_squares

# Step 4: Calculate the longest R-R interval in seconds.
longest_rr_interval_sec = total_small_squares * time_per_small_square_sec

# Step 5: Print the results and the calculation steps.
print("ECG Analysis for Longest R-R Interval")
print("-" * 40)
print(f"Paper Speed: {paper_speed_mm_per_sec} mm/sec")
print(f"Time per small square (1mm): {time_per_small_square_sec:.2f} seconds")
print("\nMeasurement of the longest interval:")
print(f"Number of large squares: {num_large_squares}")
print(f"Number of extra small squares: {num_extra_small_squares}")
print("\nCalculation:")
print(f"Total small squares = ({num_large_squares} large squares * {num_small_squares_in_large_square} small squares/large square) + {num_extra_small_squares} extra small squares = {total_small_squares} small squares")
print(f"Longest R-R interval (sec) = {total_small_squares} small squares * {time_per_small_square_sec:.2f} sec/square")
print(f"\nThe longest R-R interval is: {longest_rr_interval_sec:.2f} seconds.")