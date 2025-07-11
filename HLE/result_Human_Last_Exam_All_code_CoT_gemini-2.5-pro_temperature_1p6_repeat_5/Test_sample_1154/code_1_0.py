# Parameters for the calculation
# The number of bits in the RAM.
num_bits = 1_000_000
# The complexity of the March RAW test (operations per bit). This is chosen as it has the highest fault coverage.
complexity = 14
# The time for a single read/write operation in nanoseconds.
cycle_time_ns = 5
# The conversion factor from nanoseconds to milliseconds.
ns_to_ms_conversion = 1_000_000

# Calculate the total test duration in milliseconds
total_time_ms = (complexity * num_bits * cycle_time_ns) / ns_to_ms_conversion

# Display the test choice and the equation used for the calculation,
# showing each number involved as requested.
print("Test Chosen: March RAW (highest fault coverage, complexity 14N)")
print("Equation: (Complexity * Number of Bits * Cycle Time_ns) / ns_per_ms = Total Time_ms")
print(f"Calculation: ({complexity} * {num_bits} * {cycle_time_ns}) / {ns_to_ms_conversion} = {total_time_ms} ms")