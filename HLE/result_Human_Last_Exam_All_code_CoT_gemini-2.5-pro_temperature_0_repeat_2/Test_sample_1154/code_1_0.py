# Define the parameters for the calculation
ram_size_bits = 1000000
cycle_time_ns = 5
# March CL is the test with the highest fault coverage, with a complexity of 12N.
march_cl_complexity = 12

# Calculate the total time in nanoseconds
total_time_ns = march_cl_complexity * ram_size_bits * cycle_time_ns

# Convert the total time from nanoseconds to milliseconds
# 1 millisecond = 1,000,000 nanoseconds
total_time_ms = total_time_ns / 1000000

# Print the final equation with all the numbers and the result
print("The chosen test is March CL (12N operations).")
print("The calculation for the test duration is:")
print(f"({march_cl_complexity} * {ram_size_bits} * {cycle_time_ns}) / 1000000")
print(f"The total test duration is: {total_time_ms} ms")