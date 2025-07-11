# Step 1: Define the parameters for the calculation.
# Chosen test with the highest fault coverage.
test_name = "March RAW"
# The complexity of the March RAW test is 12N, meaning 12 operations per bit.
complexity = 12
# The size of the RAM in bits (N).
ram_size_bits = 1000000
# The time for one read/write cycle in nanoseconds (tc).
cycle_time_ns = 5

# Step 2: Calculate the total number of operations.
total_operations = complexity * ram_size_bits

# Step 3: Calculate the total test time in nanoseconds.
total_time_ns = total_operations * cycle_time_ns

# Step 4: Convert the total time from nanoseconds to milliseconds.
# 1 millisecond = 1,000,000 nanoseconds.
total_time_ms = total_time_ns / 1000000

# Print the breakdown of the calculation.
print(f"Chosen Test: {test_name}")
print(f"Test Complexity: {complexity}N")
print(f"RAM Size (N): {ram_size_bits:,} bits")
print(f"Cycle Time (tc): {cycle_time_ns} ns")
print("\nCalculating the total test duration:")
# The prompt requires outputting each number in the final equation.
print(f"Total Time (ns) = Complexity * RAM Size * Cycle Time")
print(f"Total Time (ns) = {complexity} * {ram_size_bits:,} * {cycle_time_ns} = {total_time_ns:,} ns")
print("\nConverting nanoseconds to milliseconds:")
print(f"Total Time (ms) = Total Time (ns) / 1,000,000")
print(f"Total Time (ms) = {total_time_ns:,} / 1,000,000 = {total_time_ms}")
print(f"\nThe duration of the {test_name} test is {total_time_ms} milliseconds.")
