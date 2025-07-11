# 1. Define the given parameters.
paper_speed_mm_per_sec = 25

# Each small square on the grid is 1 mm.
time_per_small_square_sec = 1 / paper_speed_mm_per_sec

# 2. Measure the longest R-R interval from the ECG strip by counting small squares.
# By inspection, the longest interval is found on the fourth strip and measures approximately 37 small squares.
longest_interval_in_small_squares = 37

# 3. Calculate the duration of the longest R-R interval in seconds.
longest_interval_in_seconds = longest_interval_in_small_squares * time_per_small_square_sec

# 4. Print the calculation steps and the final answer.
print("Step 1: Determine the time represented by one small square on the ECG.")
print(f"Paper speed = {paper_speed_mm_per_sec} mm/sec")
print(f"Time per small square = 1 mm / {paper_speed_mm_per_sec} mm/sec = {time_per_small_square_sec} seconds")
print("\nStep 2: Measure the longest R-R interval in number of small squares.")
print(f"The longest interval visually identified on the strip is {longest_interval_in_small_squares} small squares.")
print("\nStep 3: Calculate the interval duration in seconds.")
print("Duration = Number of small squares * Time per small square")
print(f"Duration = {longest_interval_in_small_squares} * {time_per_small_square_sec}")
print(f"\nThe longest R-R interval is {longest_interval_in_seconds:.2f} seconds.")
