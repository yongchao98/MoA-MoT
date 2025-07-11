# Plan:
# 1. Define the time values for small and large squares on a standard ECG (25 mm/sec).
# 2. Identify and count the number of large and small squares in the longest R-R interval from the image.
# 3. Calculate the total time in seconds.
# 4. Print the calculation steps and the final result.

# Time conversion factors for an ECG at 25 mm/sec
time_per_small_square_sec = 0.04
time_per_large_square_sec = 0.20

# Measurement from the ECG strip for the longest R-R interval
num_large_squares = 8
num_small_squares = 3

# Calculate the total duration in seconds
duration_from_large_squares = num_large_squares * time_per_large_square_sec
duration_from_small_squares = num_small_squares * time_per_small_square_sec
total_rr_interval_sec = duration_from_large_squares + duration_from_small_squares

# Output the result
print("Calculation for the longest R-R interval:")
print(f"The longest interval spans {num_large_squares} large squares and {num_small_squares} small squares.")
print(f"Time for large squares: {num_large_squares} * {time_per_large_square_sec} s = {duration_from_large_squares:.2f} s")
print(f"Time for small squares: {num_small_squares} * {time_per_small_square_sec} s = {duration_from_small_squares:.2f} s")
print(f"Total R-R interval = {duration_from_large_squares:.2f} s + {duration_from_small_squares:.2f} s = {total_rr_interval_sec:.2f} s")
print("\nThe longest R-R interval is {:.2f} seconds.".format(total_rr_interval_sec))