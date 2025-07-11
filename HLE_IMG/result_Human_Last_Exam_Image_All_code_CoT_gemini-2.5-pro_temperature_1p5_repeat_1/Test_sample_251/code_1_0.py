# The ECG paper speed is 25 mm/sec.
paper_speed_mm_per_sec = 25

# One small square on the ECG grid is 1 mm wide.
# We calculate the time it takes for the paper to advance by 1 mm.
time_per_small_square_sec = 1 / paper_speed_mm_per_sec

# By inspecting the ECG, we find the longest R-R interval.
# This interval spans 4 large squares. Since each large square is 5 small squares wide,
# the total distance is 4 * 5 = 20 small squares.
longest_rr_interval_small_squares = 20

# Calculate the duration of the longest R-R interval in seconds.
longest_rr_interval_sec = longest_rr_interval_small_squares * time_per_small_square_sec

print("To find the longest R-R interval in seconds, we use the formula:")
print("Duration (s) = Number of small squares * Time per small square (s)")
print(f"Number of small squares in the longest interval: {longest_rr_interval_small_squares}")
print(f"Time per small square: 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec} seconds")
print(f"Calculation: {longest_rr_interval_small_squares} * {time_per_small_square_sec} = {longest_rr_interval_sec}")
print(f"\nThe longest R-R interval is {longest_rr_interval_sec} seconds.")