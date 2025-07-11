# This script calculates the time to run a March RAW test on a 1Mbit RAM.

# --- Parameters ---
# RAM size in bits
ram_size_bits = 1_000_000

# Time for one read or write cycle in nanoseconds
time_per_cycle_ns = 5

# The test with the highest fault coverage is March RAW.
# Its complexity is 14N, meaning 14 operations per bit.
complexity_factor = 14

# --- Calculation ---
# 1. Calculate the total number of operations
total_operations = complexity_factor * ram_size_bits

# 2. Calculate the total test duration in nanoseconds
total_time_ns = total_operations * time_per_cycle_ns

# 3. Convert the duration from nanoseconds to milliseconds
#    (1 millisecond = 1,000,000 nanoseconds)
conversion_factor_ns_to_ms = 1_000_000
total_time_ms = total_time_ns / conversion_factor_ns_to_ms

# --- Output ---
# Print the full equation for the calculation as requested
print("The test with the highest fault coverage is March RAW (14N).")
print("Calculation of the test duration:")
print(f"({complexity_factor} operations/bit * {ram_size_bits:,} bits * {time_per_cycle_ns} ns/operation) / {conversion_factor_ns_to_ms:,} ns/ms = {total_time_ms} ms")