# Step 1: Define the parameters from the ECG strip.
paper_speed_mm_per_sec = 25  # The paper speed is 25 mm/sec.
small_square_mm = 1          # Each small square on the grid is 1 mm.

# Step 2: Calculate the time duration of a single small square.
time_per_small_square_sec = small_square_mm / paper_speed_mm_per_sec

# Step 3: Measure the longest R-R interval in terms of small squares.
# By visual inspection of the fourth rhythm strip, the longest interval
# between two consecutive R-waves spans 40 small squares.
num_small_squares_for_longest_rr = 40

# Step 4: Calculate the longest R-R interval in seconds.
longest_rr_interval_sec = num_small_squares_for_longest_rr * time_per_small_square_sec

# Output the calculation and the result.
print(f"The longest R-R interval is determined by counting the number of small squares and multiplying by the time per square.")
print(f"Number of small squares for the longest interval: {num_small_squares_for_longest_rr}")
print(f"Time per small square: {time_per_small_square_sec} seconds")
print(f"Calculation: {num_small_squares_for_longest_rr} squares * {time_per_small_square_sec} sec/square")
print(f"The longest R-R interval is: {longest_rr_interval_sec:.2f} seconds")