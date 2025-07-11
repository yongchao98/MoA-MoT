# ECG paper speed information
# 1 large square = 5 mm
# Paper speed = 25 mm/sec
# Time per large square = 5 mm / 25 mm/sec = 0.2 seconds

# Measure the longest R-R interval in large squares from the ECG strip
longest_rr_interval_large_squares = 7

# Define the time duration of one large square in seconds
time_per_large_square_sec = 0.2

# Calculate the longest R-R interval in seconds
longest_rr_interval_sec = longest_rr_interval_large_squares * time_per_large_square_sec

# Print the calculation and the final result
print("To find the longest R-R interval in seconds:")
print("1. Identify the longest interval on the ECG strip and count the number of large squares. We count 7 large squares.")
print("2. Know that at a paper speed of 25 mm/sec, each large square represents 0.2 seconds.")
print("3. Multiply the number of squares by the time per square to get the total duration.")
print(f"The calculation is: {longest_rr_interval_large_squares} large squares * {time_per_large_square_sec} seconds/large square = {longest_rr_interval_sec:.1f} seconds.")