# Step 1: Identify the test with the highest fault coverage.
# The fault coverage of a standard March test increases with its complexity (the 'x' in 'xN').
# Comparing the complexities: MSCAN (4N), MATS (~5N), March X (6N), March Y (8N),
# March C- (10N), March C (11N), March CL (12N), and March RAW (13N).
# March RAW has the highest complexity of 13N, thus it is chosen.

# Step 2: Define the parameters for the calculation.
complexity_factor = 13  # From March RAW (13N)
ram_size_bits = 1_000_000  # 1Mbit RAM
cycle_time_ns = 5  # 5 nanoseconds per cycle

# Step 3: Calculate the total test time.
# Total time (ns) = complexity_factor * ram_size_bits * cycle_time_ns
total_time_ns = complexity_factor * ram_size_bits * cycle_time_ns

# Step 4: Convert time from nanoseconds to milliseconds (1 ms = 1,000,000 ns).
total_time_ms = total_time_ns / 1_000_000

# Step 5: Print the explanation and the final equation with the result.
print("The RAM test with the highest fault coverage from the list is March RAW (13N).")
print("The time taken for this test on a 1Mbit RAM with a 5ns cycle time is calculated as follows:")
print(f"Test Duration = {complexity_factor} * {ram_size_bits} * {cycle_time_ns} ns")
print(f"Result = {int(total_time_ms)} ms")