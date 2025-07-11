# The chosen test is one with a 12N complexity, such as March CL or March RAW, 
# as it offers the highest fault coverage among the options.
# N = 1,000,000 bits (1Mbit RAM).
# tc = 5 ns (time for one read/write operation).
# The formula for total test duration in milliseconds is:
# (Complexity_Factor * N * tc_ns) / 1,000,000

# The final equation includes the following numbers:
# Complexity Factor: 12
# Number of bits (N): 1000000
# Cycle time (tc): 5
# Conversion factor ns to ms: 1000000

# The calculation is performed with these numbers directly.
# Integer division (//) is used because the result is a whole number.
test_duration_ms = (12 * 1000000 * 5) // 1000000

# Print the final result, which is the duration of the test in milliseconds.
print(test_duration_ms)