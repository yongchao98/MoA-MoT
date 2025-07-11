# 1. Define the parameters from the ECG strip.
paper_speed_mm_per_sec = 25
# Each small square on the ECG grid is 1 mm.
time_per_small_square_sec = 1 / paper_speed_mm_per_sec

# 2. Visually measure the longest R-R interval in terms of small squares.
# The longest interval is found in the second strip and spans approximately 38 small squares.
longest_interval_in_small_squares = 38

# 3. Calculate the longest R-R interval in seconds.
longest_rr_interval_sec = longest_interval_in_small_squares * time_per_small_square_sec

# 4. Print the calculation and the result.
print(f"The paper speed is {paper_speed_mm_per_sec} mm/sec.")
print(f"The time for one small square is 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec} seconds.")
print(f"The longest R-R interval spans {longest_interval_in_small_squares} small squares.")
print(f"Calculation: {longest_interval_in_small_squares} squares * {time_per_small_square_sec} sec/square")
print(f"The longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")