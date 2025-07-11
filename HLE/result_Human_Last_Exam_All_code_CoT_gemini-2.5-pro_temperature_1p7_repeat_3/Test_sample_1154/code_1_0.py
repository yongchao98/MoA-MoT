import math

# Step 1: Define the parameters for the calculation.
# The list of tests is: MSCAN, MATS, March X, March Y, March RAW, March C, March C-, March CL.
# Among these, March CL has the highest complexity (14N), providing the best fault coverage.
complexity_factor = 14
# Number of bits in the RAM (N)
num_bits = 1_000_000
# Time per read/write cycle in nanoseconds (tc)
time_per_cycle_ns = 5

# Step 2: Calculate the total number of operations.
total_operations = complexity_factor * num_bits

# Step 3: Calculate the total time in nanoseconds.
total_time_ns = total_operations * time_per_cycle_ns

# Step 4: Convert the total time from nanoseconds to milliseconds.
# 1 millisecond = 1,000,000 nanoseconds
total_time_ms = total_time_ns / 1_000_000

# Step 5: Print the equation with all the numbers, as requested.
# The formatters :, serve to add thousand separators for readability.
print("Test with highest fault coverage: March CL (14N complexity)")
print("Calculation for test duration:")
print(f"{complexity_factor} * {num_bits:,} * {time_per_cycle_ns} ns = {total_time_ns:,.0f} ns")
print(f"{total_time_ns:,.0f} ns is equal to {total_time_ms:.0f} ms.")
print("\nThe final answer is (in milliseconds):")
print(int(total_time_ms))