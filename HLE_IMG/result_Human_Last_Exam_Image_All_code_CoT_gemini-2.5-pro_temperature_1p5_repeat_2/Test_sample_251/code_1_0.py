# Plan:
# 1. Define the ECG paper speed in mm/sec.
# 2. Define the width of a small square in mm.
# 3. Calculate the time duration of one small square.
# 4. Count the number of small squares for the longest R-R interval from the ECG image.
# 5. Calculate the total time for the longest R-R interval in seconds.
# 6. Print the calculation and the result.

# 1. ECG paper speed
paper_speed_mm_per_sec = 25

# 2. Small square width (standard for ECG paper)
small_square_mm = 1

# 3. Time per small square
time_per_small_square_s = small_square_mm / paper_speed_mm_per_sec

# 4. Count the number of small squares for the longest R-R interval.
# Upon visual inspection of the ECG strip, the longest pause is on the fourth
# line of the rhythm strip. Counting the small squares between the peaks of the
# R-waves for this interval gives approximately 56 squares.
num_small_squares = 56

# 5. Calculate the longest R-R interval in seconds
longest_rr_interval_s = num_small_squares * time_per_small_square_s

# 6. Print the equation and the final answer
print(f"The longest R-R interval is calculated by multiplying the number of small squares by the time per square.")
print(f"Equation: {num_small_squares} small squares * {time_per_small_square_s} s/square = {longest_rr_interval_s:.2f} seconds")