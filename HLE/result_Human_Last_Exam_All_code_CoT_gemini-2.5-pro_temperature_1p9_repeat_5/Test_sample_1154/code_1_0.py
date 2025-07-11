# 1. Define the parameters for the calculation.
# The test with the highest fault coverage is March RAW, which has a complexity of 14N.
test_name = "March RAW"
complexity = 14  # Number of operations per bit for March RAW test.

# Given values for the RAM.
num_bits = 1000000  # Size of the RAM (1 Mbit).
cycle_time_ns = 5  # Time for one read/write operation in nanoseconds.

# Conversion factor from nanoseconds to milliseconds.
ns_per_ms = 1000000

# 2. Calculate the total test time.
# Total time = (operations per bit) * (number of bits) * (time per operation).
total_time_ns = complexity * num_bits * cycle_time_ns

# Convert the total time from nanoseconds to milliseconds.
total_time_ms = total_time_ns / ns_per_ms

# 3. Print the equation and the final answer.
# As requested, here are the numbers used in the final equation.
print(f"Test selected: {test_name} (Complexity: {complexity}N)")
print("Equation for total time in milliseconds:")
print(f"({complexity} * {num_bits} * {cycle_time_ns}) / {ns_per_ms}")
print("Result (ms):")
print(total_time_ms)