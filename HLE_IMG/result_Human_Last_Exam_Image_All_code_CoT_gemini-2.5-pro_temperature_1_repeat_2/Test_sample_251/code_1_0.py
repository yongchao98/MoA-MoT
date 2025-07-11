# Step 1: Define the constants based on the ECG strip information.
paper_speed_mm_per_sec = 25  # The paper speed is 25 mm/sec.
small_square_mm = 1          # Each small square on the grid is 1 mm wide.

# Step 2: Calculate the time duration of one small square.
time_per_small_square_sec = small_square_mm / paper_speed_mm_per_sec

# Step 3: Count the number of small squares for the longest R-R interval from the ECG.
# Upon visual inspection, the longest interval is found in the fourth strip.
# It spans approximately 39 small squares.
num_small_squares_for_longest_rr = 39

# Step 4: Calculate the longest R-R interval in seconds.
longest_rr_interval_sec = num_small_squares_for_longest_rr * time_per_small_square_sec

# Print the calculation and the final result.
print(f"The paper speed is {paper_speed_mm_per_sec} mm/sec.")
print(f"Each small square represents {time_per_small_square_sec} seconds.")
print(f"The longest R-R interval spans {num_small_squares_for_longest_rr} small squares.")
print(f"Calculation: {num_small_squares_for_longest_rr} squares * {time_per_small_square_sec} sec/square = {longest_rr_interval_sec} seconds")
print(f"The longest R-R interval is {longest_rr_interval_sec} seconds.")