# Number of bits in the RAM
N = 1000000

# Complexity of the March RAW test (14N)
complexity_factor = 14

# Time per read/write cycle in nanoseconds
tc_ns = 5

# Calculate the total number of operations
total_operations = complexity_factor * N

# Calculate the total test time in nanoseconds
total_time_ns = total_operations * tc_ns

# Convert nanoseconds to milliseconds (1 ms = 1,000,000 ns)
total_time_ms = total_time_ns / 1000000

print("The test with the highest fault coverage is March RAW (14N).")
print("Calculating the test duration for a 1Mbit RAM:")
# The user wants the equation printed, so we format the numbers with commas for readability
print(f"Equation: {complexity_factor:,} (operations/bit) * {N:,} (bits) * {tc_ns} (ns/operation) = {total_time_ns:,.0f} ns")
print(f"Converting to milliseconds: {total_time_ns:,.0f} ns / 1,000,000 = {total_time_ms:.0f} ms")
print("\nThe duration of the test is:")
print(f"{total_time_ms:.0f}")
