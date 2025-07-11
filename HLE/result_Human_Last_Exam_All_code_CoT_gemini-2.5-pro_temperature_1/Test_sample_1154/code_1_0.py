# Step 1: Define the parameters for the calculation.
# The chosen test is March CL, which has the highest fault coverage.
# Its complexity is typically 12N, meaning 12 operations per bit.
complexity_factor = 12

# The size of the RAM in bits.
ram_size_bits = 1000000

# The time for one read/write cycle in nanoseconds.
cycle_time_ns = 5

# Conversion factor from nanoseconds to milliseconds (1 ms = 1,000,000 ns).
ns_to_ms_conversion = 1000000

# Step 2: Calculate the total test duration.
# The equation is: (complexity_factor * ram_size_bits * cycle_time_ns) / ns_to_ms_conversion
total_time_ms = (complexity_factor * ram_size_bits * cycle_time_ns) / ns_to_ms_conversion

# Step 3: Print the final result in milliseconds.
# The numbers in the final equation are represented by the variables above.
print(int(total_time_ms))