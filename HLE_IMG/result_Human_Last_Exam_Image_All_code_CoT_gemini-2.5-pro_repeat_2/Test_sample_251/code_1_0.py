# Plan:
# 1. Define the time represented by one large square on the ECG grid.
#    - Paper speed = 25 mm/sec
#    - Large square width = 5 mm
#    - Time per large square = 5 mm / 25 mm/sec = 0.2 seconds
# 2. Visually determine the longest R-R interval in terms of large squares from the ECG image.
#    - By inspection, the longest interval spans approximately 8 large squares.
# 3. Calculate the total time in seconds for this interval.

# Number of large squares for the longest R-R interval
num_large_squares = 8

# Time in seconds for one large square
time_per_large_square = 0.2

# Calculate the longest R-R interval in seconds
longest_rr_interval = num_large_squares * time_per_large_square

# Print the calculation and the result
print(f"The longest R-R interval is calculated by multiplying the number of large squares by the time per large square.")
print(f"Calculation: {num_large_squares} large squares * {time_per_large_square} seconds/large square = {longest_rr_interval:.1f} seconds")