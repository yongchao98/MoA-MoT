# Step 1: Define the given parameters from the ECG strip.
paper_speed_mm_per_sec = 25

# Step 2: Calculate the time represented by one small square (1 mm).
time_per_small_square_sec = 1 / paper_speed_mm_per_sec

# Step 3: Measure the longest R-R interval found on the strip in terms of small squares.
# By visual inspection of the pause in the fourth strip, the distance is measured.
longest_interval_in_small_squares = 42

# Step 4: Calculate the longest R-R interval in seconds.
longest_rr_interval_sec = longest_interval_in_small_squares * time_per_small_square_sec

# Step 5: Print the result and the calculation.
print("ECG Analysis for Longest R-R Interval")
print(f"Paper Speed: {paper_speed_mm_per_sec} mm/sec")
print(f"Time per small square (1mm): 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec} seconds")
print(f"Longest interval measurement: {longest_interval_in_small_squares} small squares")
print("\nCalculating the longest R-R interval in seconds:")
print(f"{longest_interval_in_small_squares} squares * {time_per_small_square_sec} sec/square = {longest_rr_interval_sec:.2f} seconds")