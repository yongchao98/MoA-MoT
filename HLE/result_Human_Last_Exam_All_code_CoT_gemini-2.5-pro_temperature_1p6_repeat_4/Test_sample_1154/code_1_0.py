# The chosen test is March RAW, which has a complexity of 13N.
# N = number of bits in the RAM
# One read/write cycle takes tc.

# Number of operations per bit for March RAW test
operations_per_bit = 13

# Number of bits in the RAM
num_bits = 1000000

# Time for one read/write cycle in nanoseconds
cycle_time_ns = 5

# Conversion factor from nanoseconds to milliseconds
ns_to_ms_conversion = 1000000

# Calculate the total test duration in milliseconds
# Total time = (operations_per_bit * num_bits * cycle_time_ns) / ns_to_ms_conversion
# This simplifies to: operations_per_bit * cycle_time_ns
total_time_ms = (operations_per_bit * num_bits * cycle_time_ns) / ns_to_ms_conversion

print(f"The equation to calculate the test duration for a {num_bits:,} bit RAM is:")
print(f"({operations_per_bit} operations/bit * {num_bits:,} bits * {cycle_time_ns} ns/operation) / {ns_to_ms_conversion:,} ns/ms")
print("\nThe final duration of the test is:")
print(int(total_time_ms))