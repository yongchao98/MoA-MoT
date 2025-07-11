# ECG parameters
paper_speed_mm_per_sec = 25

# Grid measurements
mm_per_small_square = 1
time_per_mm = 1 / paper_speed_mm_per_sec # in seconds

# Time per small square
time_per_small_square = mm_per_small_square * time_per_mm

# Visually identify and count the number of small squares in the longest R-R interval.
# The longest interval is found in the fourth strip and measures approximately 42 small squares.
num_small_squares_longest_interval = 42

# Calculate the longest R-R interval in seconds
longest_rr_interval_sec = num_small_squares_longest_interval * time_per_small_square

print("Calculation Steps:")
print(f"Paper speed: {paper_speed_mm_per_sec} mm/sec")
print(f"Time per small square (1 mm): 1 sec / {paper_speed_mm_per_sec} mm = {time_per_small_square:.2f} seconds")
print(f"Longest R-R interval measurement: {num_small_squares_longest_interval} small squares")
print(f"Longest R-R interval in seconds = {num_small_squares_longest_interval} squares * {time_per_small_square:.2f} sec/square = {longest_rr_interval_sec:.2f} seconds")