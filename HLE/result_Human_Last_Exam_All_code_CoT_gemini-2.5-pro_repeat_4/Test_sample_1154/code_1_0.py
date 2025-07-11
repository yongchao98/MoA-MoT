# Define the parameters for the calculation
test_name = "March RAW"
complexity_factor = 13  # Based on the 13N complexity of March RAW
num_bits = 1000000      # 1 Mbit RAM
cycle_time_ns = 5       # Time for one read/write cycle in nanoseconds

# Calculate the total number of operations
total_operations = complexity_factor * num_bits

# Calculate the total test time in nanoseconds
total_time_ns = total_operations * cycle_time_ns

# Convert the total time from nanoseconds to milliseconds
# 1 millisecond = 1,000,000 nanoseconds
total_time_ms = total_time_ns / 1000000

# Print the breakdown of the calculation
print(f"Test with highest fault coverage: {test_name} (Complexity: {complexity_factor}N)")
print("--- Calculation ---")
# The user requested to output each number in the final equation
print(f"Total Time (ms) = ({complexity_factor} * {num_bits:,} * {cycle_time_ns}) / 1,000,000")
print(f"Total Time (ms) = {total_time_ns:,} ns / 1,000,000")
print(f"The total test duration is: {total_time_ms} ms")
