# ECG parameters
paper_speed_mm_per_sec = 25

# Time calculation for one small square (1mm)
time_per_small_square_sec = 1 / paper_speed_mm_per_sec

# Measurement of the longest R-R interval from the ECG strip
# The longest interval spans 7 large squares and 2 small squares.
# Since 1 large square = 5 small squares, the total is:
num_large_squares = 7
num_small_squares_extra = 2
total_small_squares = (num_large_squares * 5) + num_small_squares_extra

# Calculate the longest R-R interval in seconds
longest_rr_interval_sec = total_small_squares * time_per_small_square_sec

# Output the calculation
print(f"The longest R-R interval spans {total_small_squares} small squares.")
print(f"The paper speed is {paper_speed_mm_per_sec} mm/sec, so each small square represents {time_per_small_square_sec:.2f} seconds.")
print(f"Calculation: {total_small_squares} squares * {time_per_small_square_sec:.2f} sec/square = {longest_rr_interval_sec:.2f} seconds.")
print(f"\nThe longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")
