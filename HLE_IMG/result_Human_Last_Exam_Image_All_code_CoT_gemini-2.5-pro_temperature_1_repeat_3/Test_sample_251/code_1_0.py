# Step 1: Define the parameters from the ECG strip.
paper_speed_mm_per_sec = 25
small_square_mm = 1
large_square_mm = 5

# Step 2: Calculate the time duration for one small square.
time_per_small_square_sec = small_square_mm / paper_speed_mm_per_sec

# Step 3: Visually identify and measure the longest R-R interval in grid squares.
# The longest interval is found in the fourth row.
# It spans 8 large squares and 2 small squares.
num_large_squares = 8
num_small_squares_extra = 2

# Step 4: Convert the measurement to total small squares.
small_squares_per_large_square = large_square_mm / small_square_mm
total_small_squares = (num_large_squares * small_squares_per_large_square) + num_small_squares_extra

# Step 5: Calculate the final R-R interval in seconds.
longest_rr_interval_sec = total_small_squares * time_per_small_square_sec

# Print the calculation steps and the final answer.
print("ECG Analysis for Longest R-R Interval")
print(f"Paper Speed: {paper_speed_mm_per_sec} mm/sec")
print(f"Time per small square: {time_per_small_square_sec:.2f} seconds")
print("\nMeasurement of the longest interval:")
print(f"Number of large squares: {num_large_squares}")
print(f"Number of extra small squares: {num_small_squares_extra}")
print(f"\nCalculation:")
print(f"Total small squares = ({num_large_squares} * {int(small_squares_per_large_square)}) + {num_small_squares_extra} = {int(total_small_squares)}")
print(f"Longest R-R Interval (sec) = {int(total_small_squares)} small squares * {time_per_small_square_sec:.2f} sec/square")
print(f"Result: {longest_rr_interval_sec:.2f} seconds")
