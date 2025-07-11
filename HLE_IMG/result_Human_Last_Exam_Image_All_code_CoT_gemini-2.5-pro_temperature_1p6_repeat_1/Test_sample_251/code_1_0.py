# The standard speed for an ECG rhythm strip is 25 mm/sec.
# The grid is made of 1mm small squares and 5mm large squares.
# Time per small square = 1 mm / 25 mm/sec = 0.04 seconds.
# Time per large square = 5 mm / 25 mm/sec = 0.2 seconds.

# By observing the ECG strip, the longest interval between two R-waves (R-R interval)
# is found on the fourth line of the trace.
# We can measure this interval by counting the number of large squares between the two R-wave peaks.
num_large_squares = 8

# The time for one large square is 0.2 seconds.
time_per_large_square_sec = 0.2

# Calculate the longest R-R interval in seconds.
longest_rr_interval_sec = num_large_squares * time_per_large_square_sec

print(f"To find the longest R-R interval, we measure the distance in large squares and multiply by the time per large square.")
print(f"Longest interval in large squares: {num_large_squares}")
print(f"Time per large square (sec): {time_per_large_square_sec}")
print(f"Calculation: {num_large_squares} * {time_per_large_square_sec} = {longest_rr_interval_sec}")
print(f"The longest R-R interval is {longest_rr_interval_sec} seconds.")
