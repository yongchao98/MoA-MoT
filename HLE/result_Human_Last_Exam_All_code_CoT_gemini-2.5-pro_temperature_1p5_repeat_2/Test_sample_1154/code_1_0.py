# The RAM test with the highest fault coverage from the list is March RAW.
# Its complexity is 14N, where N is the number of memory bits.

# --- Variables ---
# Complexity of the March RAW test (operations per bit)
complexity = 14

# Number of bits in the RAM (N)
ram_size_bits = 1_000_000

# Time for one read/write cycle in nanoseconds (tc)
cycle_time_ns = 5

# --- Calculation ---
# Total number of operations
total_operations = complexity * ram_size_bits

# Total test time in nanoseconds
total_time_ns = total_operations * cycle_time_ns

# Convert total time from nanoseconds to milliseconds
# 1 ms = 1,000,000 ns
ns_to_ms_conversion_factor = 1_000_000
total_time_ms = total_time_ns / ns_to_ms_conversion_factor

# --- Output ---
print("Chosen Test: March RAW (14N complexity)")
print(f"RAM Size (N): {ram_size_bits:,} bits")
print(f"Cycle Time (tc): {cycle_time_ns} ns")
print("\nCalculating the total test duration:")
# As requested, printing the equation with the numbers
print(f"Equation: ({complexity} ops/bit * {ram_size_bits:,} bits * {cycle_time_ns} ns/op) / {ns_to_ms_conversion_factor:,} ns/ms")
print(f"Result: {total_time_ms} ms")
<<<70>>>