# Parameters for the RAM test time calculation.

# N represents the size of the RAM in bits.
N = 1000000

# The complexity factor for the March CL test, which is 12N.
# This test is chosen for its high fault coverage.
complexity_factor = 12

# tc_ns is the time for a single read/write cycle in nanoseconds.
tc_ns = 5

# To convert the final result from nanoseconds to milliseconds,
# we divide by this conversion factor.
ns_to_ms_conversion = 1000000

# Calculate the total test time in milliseconds.
# The formula is (Complexity Factor * RAM Size * Cycle Time) / Conversion Factor
total_time_ms = (complexity_factor * N * tc_ns) / ns_to_ms_conversion

# As requested, here are the numbers in the final equation.
print("Equation used for the calculation (time in ms):")
print(f"({complexity_factor} * {N} * {tc_ns}) / {ns_to_ms_conversion}")

# Print the final result.
print("\nTest duration in milliseconds:")
print(total_time_ms)