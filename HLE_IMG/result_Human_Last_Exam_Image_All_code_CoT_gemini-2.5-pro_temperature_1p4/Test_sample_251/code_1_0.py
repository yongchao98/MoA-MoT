# The ECG paper speed is 25 mm/sec.
# Each small square on the grid is 1 mm.
speed_mm_per_sec = 25

# First, calculate the time duration of one small square (1 mm).
time_per_small_square_sec = 1 / speed_mm_per_sec

# Next, identify the longest R-R interval on the strip by counting the number of small squares.
# Visual inspection shows a long pause in the fourth row.
# Let's count the number of small squares for this longest interval.
# The interval spans 6 large squares and 2 small squares.
# Since one large square = 5 small squares:
num_large_squares = 6
num_small_squares_extra = 2
num_small_squares_total = (num_large_squares * 5) + num_small_squares_extra

# Finally, calculate the duration of the longest R-R interval in seconds.
longest_rr_interval_sec = num_small_squares_total * time_per_small_square_sec

print(f"The ECG paper speed is {speed_mm_per_sec} mm/sec.")
print(f"The time for one small square is 1 / {speed_mm_per_sec} = {time_per_small_square_sec:.2f} seconds.")
print(f"The longest R-R interval spans {num_small_squares_total} small squares.")
print(f"Calculation: {num_small_squares_total} small squares * {time_per_small_square_sec:.2f} seconds/square = {longest_rr_interval_sec:.2f} seconds.")
print(f"The longest R-R interval is {longest_rr_interval_sec:.2f} seconds.")