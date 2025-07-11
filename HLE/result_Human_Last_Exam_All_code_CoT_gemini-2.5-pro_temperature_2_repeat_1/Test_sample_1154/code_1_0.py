# Step 1: Define the parameters for the calculation.
# We choose March CL as it generally has the highest fault coverage among the listed options.
# Its algorithmic complexity is 12N, meaning 12 operations per memory bit.

# Number of bits in the RAM
num_bits = 1_000_000

# Time for one read/write cycle in nanoseconds
time_per_cycle_ns = 5

# Algorithmic complexity factor for March CL
complexity_factor = 12

# Step 2: Calculate the total number of operations.
total_operations = complexity_factor * num_bits

# Step 3: Calculate the total time for the test in nanoseconds.
total_time_ns = total_operations * time_per_cycle_ns

# Step 4: Convert the total time from nanoseconds to milliseconds.
# 1 millisecond = 1,000,000 nanoseconds.
total_time_ms = total_time_ns / 1_000_000

# Step 5: Print the explanation and the final result, showing the numbers used.
print("Selected Test: March CL (Complexity: 12N)")
print("Reason: Offers high and comprehensive fault coverage compared to the other listed tests.")
print("\n--- Calculation ---")
print(f"Number of bits (N): {num_bits:,}")
print(f"Complexity Factor: {complexity_factor}")
print(f"Time per cycle (tc): {time_per_cycle_ns} ns")
print("\nEquation for total time:")
print("Total Time = (Complexity Factor * Number of Bits * Time per Cycle)")
print(f"Total Time = {complexity_factor} * {num_bits:,} * {time_per_cycle_ns} ns")
print(f"Total Time = {total_time_ns:,} ns")
print(f"\nTo convert to milliseconds (ms), we divide by 1,000,000:")
print(f"Total Time = {total_time_ns:,} ns / 1,000,000 = {int(total_time_ms)} ms")